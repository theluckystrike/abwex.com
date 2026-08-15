/* Shared mobile nav toggle.

   Every page ships the .mobile-toggle button, but before this file only index.html
   loaded any JavaScript that wired it up. That stayed invisible while the stylesheet
   was leaving the nav permanently visible on mobile. Once the nav collapses properly
   the missing handler would have left 56 pages with no mobile navigation at all, so
   the handler has to live somewhere every page can load it.

   Written to NASA Power of 10 habits. There is no recursion, the one ancestor walk
   has a fixed upper bound, every query result is checked before use, and the module
   allocates nothing after wiring.

   Idempotent on purpose. index.html also runs initNav in app.js and the FDR page has
   its own copy. Binding twice would toggle twice and cancel out, so the first one to
   run claims the button through data-nav-wired and the others return. */
(function () {
  'use strict';

  var MAX_ANCESTOR_HOPS = 8;   /* fixed bound on the walk from click target to <a> */
  var OPEN_CLASS = 'open';
  var WIRED_FLAG = '1';

  /* Walk up from a click target looking for the anchor that owns it. Bounded rather
     than using closest() so the loop has a visible fixed ceiling, and so a link that
     later wraps an icon or a span still closes the menu. */
  function anchorFor(node, root) {
    if (!node || !root) { return null; }
    var cur = node;
    var hops = 0;
    while (cur && cur !== root && hops < MAX_ANCESTOR_HOPS) {
      if (cur.tagName === 'A') { return cur; }
      cur = cur.parentElement;
      hops += 1;
    }
    return null;
  }

  function setExpanded(toggle, isOpen) {
    if (!toggle) { return false; }
    if (typeof isOpen !== 'boolean') { return false; }
    toggle.setAttribute('aria-expanded', isOpen ? 'true' : 'false');
    return true;
  }

  function wire() {
    var toggle = document.querySelector('.mobile-toggle');
    var nav = document.querySelector('header nav') || document.querySelector('nav');
    /* Both queries can return null on a page without the shell. Nothing to wire. */
    if (!toggle || !nav) { return false; }
    /* Another script already owns this button. */
    if (toggle.dataset && toggle.dataset.navWired === WIRED_FLAG) { return false; }
    if (toggle.dataset) { toggle.dataset.navWired = WIRED_FLAG; }

    setExpanded(toggle, false);

    toggle.addEventListener('click', function () {
      var isOpen = nav.classList.toggle(OPEN_CLASS);
      setExpanded(toggle, isOpen === true);
    });

    /* Close after following a link, so the menu is not left open behind the next page
       on a browser that restores DOM state on back navigation. */
    nav.addEventListener('click', function (e) {
      if (!anchorFor(e.target, nav)) { return; }
      nav.classList.remove(OPEN_CLASS);
      setExpanded(toggle, false);
    });

    return true;
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', wire);
  } else {
    wire();
  }
}());
