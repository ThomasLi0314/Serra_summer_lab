// LCStool docs — sidebar navigation logic
(function() {
  // Mark the currently-open page and its parent as active
  document.addEventListener('DOMContentLoaded', function() {
    var here = (location.pathname.split('/').pop() || 'index.html').toLowerCase();
    if (here === '') here = 'index.html';

    document.querySelectorAll('aside.sidebar nav a').forEach(function(a) {
      var href = (a.getAttribute('href') || '').toLowerCase();
      // strip leading ./ or hash-only anchors
      var hrefFile = href.split('#')[0];
      if (hrefFile === here) {
        a.classList.add('active');
        // open parent li
        var li = a.closest('li');
        while (li) {
          li.classList.add('active', 'open');
          li = li.parentElement.closest('li');
        }
      }
    });

    // Toggle submenus on click of parent-link arrow
    document.querySelectorAll('aside.sidebar li.has-sub > a').forEach(function(a) {
      a.addEventListener('click', function(e) {
        var li = a.parentElement;
        // Only toggle if user clicked the arrow area (right side) OR if href is '#'
        var href = a.getAttribute('href') || '';
        if (href === '#' || href === '') {
          e.preventDefault();
          li.classList.toggle('open');
        }
        // if href points somewhere real, let navigation happen
      });
    });
  });
})();
