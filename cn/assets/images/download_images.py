#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
图片下载脚本 - 自组织与模式形成课程笔记

功能：
1. 扫描 docs/cn/ 目录下所有 Markdown 文件
2. 提取文件中的图片 URL（支持 mdnice、imgur 等常见图床）
3. 下载图片到本地 assets/images/ 目录
4. 支持增量下载（跳过已下载的图片）
5. 生成下载日志，记录已处理的文件

使用方法：
    python download_images.py

作者：Zhihang Liu
日期：2026-02
"""

import os
import re
import sys
import json
import hashlib
import requests
from pathlib import Path
from datetime import datetime
from urllib.parse import urlparse, unquote
from typing import Dict, List, Set, Tuple

# 强制刷新输出
def plog(msg: str):
    print(msg, flush=True)

# ==================== 配置 ====================

# 脚本所在目录
SCRIPT_DIR = Path(__file__).parent.absolute()

# Markdown 文件目录（docs/cn/assets/images -> docs/cn）
DOCS_DIR = SCRIPT_DIR.parent.parent

# 图片保存目录
IMAGES_DIR = SCRIPT_DIR

# 下载日志文件（记录已处理的文件和图片）
LOG_FILE = SCRIPT_DIR / "download_log.json"

# 支持的图片格式
IMAGE_EXTENSIONS = {'.jpg', '.jpeg', '.png', '.gif', '.webp', '.svg', '.bmp'}

# 请求超时时间（秒）
REQUEST_TIMEOUT = 15

# 请求头（模拟浏览器）
HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
    'Accept': 'image/avif,image/webp,image/apng,image/svg+xml,image/*,*/*;q=0.8',
    'Accept-Language': 'zh-CN,zh;q=0.9,en;q=0.8',
}

# ==================== 工具函数 ====================

def load_log() -> Dict:
    """加载下载日志"""
    if LOG_FILE.exists():
        try:
            with open(LOG_FILE, 'r', encoding='utf-8') as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError):
            pass
    return {
        "last_update": None,
        "processed_files": {},  # {文件名: {最后修改时间, 图片数量}}
        "downloaded_urls": {}   # {url: 本地文件名}
    }

def save_log(log: Dict) -> None:
    """保存下载日志"""
    log["last_update"] = datetime.now().isoformat()
    with open(LOG_FILE, 'w', encoding='utf-8') as f:
        json.dump(log, f, ensure_ascii=False, indent=2)

def get_file_hash(filepath: Path) -> str:
    """获取文件的 MD5 哈希值（用于检测文件是否修改）"""
    hasher = hashlib.md5()
    with open(filepath, 'rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            hasher.update(chunk)
    return hasher.hexdigest()

def extract_image_urls(content: str) -> List[str]:
    """
    从 Markdown 内容中提取图片 URL
    
    支持的格式：
    - ![alt](url)
    - ![alt](url "title")
    - <img src="url" ...>
    """
    urls = []
    
    # Markdown 图片语法: ![alt](url) 或 ![alt](url "title")
    md_pattern = r'!\[([^\]]*)\]\(([^)\s]+)(?:\s+"[^"]*")?\)'
    for match in re.finditer(md_pattern, content):
        url = match.group(2)
        if url.startswith('http'):
            urls.append(url)
    
    # HTML img 标签: <img src="url" ...>
    html_pattern = r'<img[^>]+src=["\']([^"\']+)["\']'
    for match in re.finditer(html_pattern, content, re.IGNORECASE):
        url = match.group(1)
        if url.startswith('http'):
            urls.append(url)
    
    return urls

def url_to_filename(url: str, index: int, note_name: str) -> str:
    """
    将 URL 转换为本地文件名
    
    格式：{笔记编号}_{序号}_{原始文件名或哈希}
    例如：01_001_853e3a67-9126-4927-907c-9f48eade8561.png
    """
    # 解析 URL
    parsed = urlparse(url)
    path = unquote(parsed.path)
    
    # 获取原始文件名
    original_name = Path(path).name
    
    # 获取扩展名
    ext = Path(path).suffix.lower()
    if ext not in IMAGE_EXTENSIONS:
        ext = '.png'  # 默认扩展名
    
    # 提取笔记编号（从文件名开头提取数字）
    note_num_match = re.match(r'^(\d+)', note_name)
    note_num = note_num_match.group(1).zfill(2) if note_num_match else '00'
    
    # 生成文件名
    # 如果原始文件名有意义（不是纯哈希），保留一部分
    if len(original_name) > 40:
        # 可能是 UUID 或哈希，截取一部分
        name_part = Path(original_name).stem[:20]
    else:
        name_part = Path(original_name).stem
    
    # 清理文件名中的特殊字符
    name_part = re.sub(r'[<>:"/\\|?*]', '_', name_part)
    
    return f"{note_num}_{str(index).zfill(3)}_{name_part}{ext}"

def download_image(url: str, save_path: Path) -> bool:
    """
    下载图片到指定路径
    
    返回：是否成功
    """
    try:
        response = requests.get(
            url, 
            headers=HEADERS, 
            timeout=REQUEST_TIMEOUT,
            stream=True
        )
        response.raise_for_status()
        
        # 检查内容类型
        content_type = response.headers.get('Content-Type', '')
        if not content_type.startswith('image/'):
            print(f"  ⚠ 警告：非图片内容类型 ({content_type})")
        
        # 保存图片
        with open(save_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        return True
        
    except requests.exceptions.RequestException as e:
        print(f"  ✗ 下载失败：{e}")
        return False

# ==================== 主逻辑 ====================

def process_markdown_file(md_file: Path, log: Dict) -> Tuple[int, int]:
    """
    处理单个 Markdown 文件
    
    返回：(新下载数量, 跳过数量)
    """
    note_name = md_file.stem
    print(f"\n📄 处理文件：{note_name}")
    
    # 读取文件内容
    with open(md_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # 提取图片 URL
    urls = extract_image_urls(content)
    
    if not urls:
        print("  (无图片)")
        return 0, 0
    
    print(f"  发现 {len(urls)} 张图片")
    
    downloaded = 0
    skipped = 0
    
    for i, url in enumerate(urls, 1):
        # 检查是否已下载
        if url in log["downloaded_urls"]:
            local_file = log["downloaded_urls"][url]
            if (IMAGES_DIR / local_file).exists():
                skipped += 1
                continue
        
        # 生成本地文件名
        filename = url_to_filename(url, i, note_name)
        save_path = IMAGES_DIR / filename
        
        # 避免文件名冲突
        counter = 1
        while save_path.exists() and url not in log["downloaded_urls"]:
            stem = Path(filename).stem
            ext = Path(filename).suffix
            filename = f"{stem}_{counter}{ext}"
            save_path = IMAGES_DIR / filename
            counter += 1
        
        # 下载图片
        print(f"  [{i}/{len(urls)}] 下载：{filename[:50]}...")
        if download_image(url, save_path):
            log["downloaded_urls"][url] = filename
            downloaded += 1
            print(f"  ✓ 成功")
        else:
            print(f"  ✗ 失败")
    
    return downloaded, skipped

def main():
    """主函数"""
    print("=" * 60)
    print("📥 图片下载脚本 - 自组织与模式形成课程笔记")
    print("=" * 60)
    
    # 确保图片目录存在
    IMAGES_DIR.mkdir(parents=True, exist_ok=True)
    
    # 加载日志
    log = load_log()
    
    # 获取所有 Markdown 文件
    md_files = sorted(DOCS_DIR.glob("*.md"))
    
    if not md_files:
        print(f"\n⚠ 未找到 Markdown 文件：{DOCS_DIR}")
        return
    
    print(f"\n找到 {len(md_files)} 个 Markdown 文件")
    
    # 统计
    total_downloaded = 0
    total_skipped = 0
    processed_count = 0
    
    for md_file in md_files:
        # 获取文件修改时间
        mtime = md_file.stat().st_mtime
        file_key = md_file.name
        
        # 检查文件是否已处理且未修改
        if file_key in log["processed_files"]:
            last_mtime = log["processed_files"][file_key].get("mtime", 0)
            if mtime <= last_mtime:
                print(f"\n⏭ 跳过未修改文件：{md_file.stem}")
                continue
        
        # 处理文件
        downloaded, skipped = process_markdown_file(md_file, log)
        
        # 更新日志
        log["processed_files"][file_key] = {
            "mtime": mtime,
            "image_count": downloaded + skipped,
            "last_processed": datetime.now().isoformat()
        }
        
        total_downloaded += downloaded
        total_skipped += skipped
        processed_count += 1
    
    # 保存日志
    save_log(log)
    
    # 输出统计
    print("\n" + "=" * 60)
    print("📊 统计信息")
    print("=" * 60)
    print(f"  处理文件数：{processed_count}")
    print(f"  新下载图片：{total_downloaded}")
    print(f"  跳过已有图片：{total_skipped}")
    print(f"  图片保存目录：{IMAGES_DIR}")
    print(f"  日志文件：{LOG_FILE}")
    print("=" * 60)

if __name__ == "__main__":
    main()
