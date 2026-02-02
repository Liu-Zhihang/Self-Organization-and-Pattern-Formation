// Video initialization for proper playback
document.addEventListener('DOMContentLoaded', function() {
  const videos = document.querySelectorAll('video');
  videos.forEach(video => {
    video.setAttribute('playsinline', '');
    video.setAttribute('preload', 'metadata');
  });
});
