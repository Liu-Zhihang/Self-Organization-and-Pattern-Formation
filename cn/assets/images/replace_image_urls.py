#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
图片路径替换脚本 - 将云端图片链接替换为本地路径

功能：
1. 扫描 docs/cn/ 目录下所有 Markdown 文件
2. 将云端图片 URL 替换为本地相对路径
3. 生成替换报告

使用方法：
    python replace_image_urls.py

作者：Zhihang Liu
日期：2026-02
"""

import os
import re
import json
from pathlib import Path
from typing import Dict, List, Tuple

# ==================== 配置 ====================

# 脚本所在目录
SCRIPT_DIR = Path(__file__).parent.absolute()

# Markdown 文件目录
DOCS_DIR = SCRIPT_DIR.parent.parent

# 图片保存目录
IMAGES_DIR = SCRIPT_DIR

# 下载日志文件
LOG_FILE = SCRIPT_DIR / "download_log.json"

# 本地图片路径前缀（相对于 Markdown 文件的路径）
LOCAL_IMAGE_PREFIX = "assets/images/"

# ==================== 主逻辑 ====================

def load_url_mapping() -> Dict[str, str]:
    """
    从下载日志和实际文件构建 URL -> 本地文件名的映射
    """
    mapping = {}
    
    # 从日志加载
    if LOG_FILE.exists():
        with open(LOG_FILE, 'r', encoding='utf-8') as f:
            log = json.load(f)
            for url, filename in log.get("downloaded_urls", {}).items():
                # 去掉 _1 后缀（如果有）
                clean_filename = re.sub(r'_1(\.[^.]+)$', r'\1', filename)
                # 检查文件是否存在
                if (IMAGES_DIR / clean_filename).exists():
                    mapping[url] = clean_filename
                elif (IMAGES_DIR / filename).exists():
                    mapping[url] = filename
    
    # 从实际文件补充映射
    for img_file in IMAGES_DIR.glob("*.*"):
        if img_file.suffix.lower() in {'.png', '.jpg', '.jpeg', '.gif', '.webp', '.svg'}:
            # 从文件名提取原始 UUID
            match = re.search(r'_\d{3}_([a-f0-9-]+)(?:_\d+)?\.', img_file.name)
            if match:
                uuid = match.group(1)
                # 尝试匹配常见图床 URL 模式
                for user_id in ['141272', '129153']:
                    for ext in ['.png', '.jpg', '.jpeg', '.gif']:
                        url = f"https://files.mdnice.com/user/{user_id}/{uuid}{ext}"
                        if url not in mapping:
                            mapping[url] = img_file.name
    
    return mapping

def replace_urls_in_file(md_file: Path, url_mapping: Dict[str, str]) -> Tuple[int, List[str]]:
    """
    替换单个 Markdown 文件中的图片 URL
    
    返回：(替换数量, 未找到的 URL 列表)
    """
    with open(md_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    replaced_count = 0
    not_found = []
    new_content = content
    
    # 查找所有图片 URL
    pattern = r'!\[([^\]]*)\]\((https://[^)\s]+)(?:\s+"[^"]*")?\)'
    
    for match in re.finditer(pattern, content):
        alt_text = match.group(1)
        url = match.group(2)
        
        if url in url_mapping:
            local_filename = url_mapping[url]
            local_path = LOCAL_IMAGE_PREFIX + local_filename
            # 替换 URL
            old_str = match.group(0)
            new_str = f'![{alt_text}]({local_path})'
            new_content = new_content.replace(old_str, new_str, 1)
            replaced_count += 1
        else:
            not_found.append(url)
    
    # 写回文件
    if replaced_count > 0:
        with open(md_file, 'w', encoding='utf-8') as f:
            f.write(new_content)
    
    return replaced_count, not_found

def main():
    print("=" * 60)
    print("🔄 图片路径替换脚本")
    print("=" * 60)
    
    # 加载 URL 映射
    url_mapping = load_url_mapping()
    print(f"\n已加载 {len(url_mapping)} 个 URL 映射")
    
    # 获取所有 Markdown 文件
    md_files = sorted(DOCS_DIR.glob("*.md"))
    
    if not md_files:
        print(f"\n⚠ 未找到 Markdown 文件：{DOCS_DIR}")
        return
    
    print(f"找到 {len(md_files)} 个 Markdown 文件")
    
    # 统计
    total_replaced = 0
    all_not_found = []
    
    for md_file in md_files:
        print(f"\n📄 处理：{md_file.name}")
        replaced, not_found = replace_urls_in_file(md_file, url_mapping)
        print(f"   替换了 {replaced} 个链接")
        if not_found:
            print(f"   ⚠ 未找到本地文件的 URL: {len(not_found)} 个")
            all_not_found.extend([(md_file.name, url) for url in not_found])
        total_replaced += replaced
    
    # 输出统计
    print("\n" + "=" * 60)
    print("📊 统计信息")
    print("=" * 60)
    print(f"  总替换数：{total_replaced}")
    print(f"  未找到的 URL 数：{len(all_not_found)}")
    
    if all_not_found:
        print("\n⚠ 以下 URL 没有找到对应的本地文件（需要手动下载）：")
        for filename, url in all_not_found:
            print(f"  [{filename}] {url}")

if __name__ == "__main__":
    main()
