// Language navigation helper
document.addEventListener('DOMContentLoaded', function() {
  // Add language switcher functionality if needed
  const languageLinks = document.querySelectorAll('[hreflang]');
  languageLinks.forEach(link => {
    link.addEventListener('click', function(e) {
      // Store language preference
      localStorage.setItem('preferredLanguage', this.hreflang);
    });
  });
});
