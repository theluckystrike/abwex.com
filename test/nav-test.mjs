/* Regression test for the shared mobile nav.

   Guards the three defects found on 2026-08-14.
     1. An unscoped `nav { display: flex }` was added to style.css AFTER the
        max-width 768px block, so it re-set display on mobile and the menu sat
        open on top of the h1 on every page.
     2. 56 of 57 pages shipped .mobile-toggle with no JavaScript wiring it, so
        once the CSS was fixed those pages would have had no mobile nav at all.
     3. The close-on-link handler matched e.target.tagName, which breaks the
        moment a link wraps an icon.

   Pure logic is pulled out of the IIFE and evaluated standalone, matching
   fdr-calculator-test.mjs. The rest are static checks over the built files. */
import fs from 'fs';
import path from 'path';

const ROOT = '/Users/mike/Desktop/sites/abwex-repo';
const navSrc = fs.readFileSync(path.join(ROOT, 'assets/js/nav.js'), 'utf8');
const cssSrc = fs.readFileSync(path.join(ROOT, 'assets/style.css'), 'utf8');

let fails = 0;
function ok(name, cond, detail) {
  if (!cond) { fails += 1; console.log('  FAIL ' + name + (detail ? '  | ' + detail : '')); }
}

/* ---- 1. pure logic, extracted from the IIFE ------------------------------ */
const names = ['anchorFor', 'setExpanded'];
let code = '';
for (const n of names) {
  const m = navSrc.match(new RegExp('\\n  function ' + n + '\\s*\\([^)]*\\)\\s*\\{'));
  if (!m) { fails += 1; console.log('  FAIL missing function ' + n); continue; }
  let i = m.index + m[0].length - 1, depth = 0;
  const start = m.index;
  while (i < navSrc.length) {
    if (navSrc[i] === '{') depth += 1;
    else if (navSrc[i] === '}') { depth -= 1; if (!depth) break; }
    i += 1;
  }
  code += navSrc.slice(start, i + 1) + '\n';
}
const F = new Function('var MAX_ANCESTOR_HOPS=8;\n' + code + '\nreturn {' + names.join(',') + '};')();

/* a tiny fake element graph, enough for the ancestor walk */
function el(tag, parent) { return { tagName: tag, parentElement: parent || null }; }
const root = el('NAV');
const anchor = el('A', root);
const span = el('SPAN', anchor);
const svg = el('SVG', span);

ok('anchorFor finds a direct anchor', F.anchorFor(anchor, root) === anchor);
ok('anchorFor finds anchor through a nested span', F.anchorFor(span, root) === anchor);
ok('anchorFor finds anchor two levels deep', F.anchorFor(svg, root) === anchor);
ok('anchorFor returns null for a non-anchor branch', F.anchorFor(el('DIV', root), root) === null);
ok('anchorFor returns null on null node', F.anchorFor(null, root) === null);
ok('anchorFor returns null on null root', F.anchorFor(anchor, null) === null);
ok('anchorFor stops at the root', F.anchorFor(root, root) === null);

/* the walk must terminate on a chain longer than the bound, not hang */
let deep = el('DIV');
for (let i = 0; i < 50; i += 1) { deep = el('DIV', deep); }
ok('anchorFor is bounded on a deep chain', F.anchorFor(deep, root) === null);

/* setExpanded checks its inputs and reports success */
const fake = { attrs: {}, setAttribute(k, v) { this.attrs[k] = v; } };
ok('setExpanded rejects a null toggle', F.setExpanded(null, true) === false);
ok('setExpanded rejects a non-boolean', F.setExpanded(fake, 'yes') === false);
ok('setExpanded writes true', F.setExpanded(fake, true) === true && fake.attrs['aria-expanded'] === 'true');
ok('setExpanded writes false', F.setExpanded(fake, false) === true && fake.attrs['aria-expanded'] === 'false');

/* ---- 2. the CSS cascade defect -------------------------------------------- */
const mobileBlock = cssSrc.indexOf('@media (max-width: 768px)');
ok('style.css still has the mobile nav block', mobileBlock > -1);
/* any bare `nav {` rule after that block, outside a min-width media query, would
   re-break the collapse the same way the V5 rule did */
const after = cssSrc.slice(mobileBlock);
const bareNavAfter = /\n\s*nav\s*\{[^}]*display\s*:\s*flex/.test(after)
  && !/@media \(min-width: 769px\)[\s\S]{0,200}nav\s*\{/.test(after);
ok('no unscoped nav display:flex after the mobile block', !bareNavAfter,
  'an unscoped nav rule after the media query re-opens the mobile menu site-wide');
ok('mobile block still collapses the nav', /nav\s*\{[^}]*display:\s*none/.test(cssSrc));
ok('open class still shows the nav', /nav\.open\s*\{[^}]*display:\s*flex/.test(cssSrc));

/* ---- 3. every page with the button loads a handler ------------------------ */
function walk(dir, out) {
  const entries = fs.readdirSync(dir, { withFileTypes: true });
  for (const e of entries) {
    if (e.name === '.git' || e.name === 'node_modules' || e.name === 'test') continue;
    const p = path.join(dir, e.name);
    if (e.isDirectory()) walk(p, out);
    else if (e.name.endsWith('.html')) out.push(p);
  }
  return out;
}
const pages = walk(ROOT, []);
const withToggle = pages.filter((p) => fs.readFileSync(p, 'utf8').includes('mobile-toggle'));
const unwired = withToggle.filter((p) => !fs.readFileSync(p, 'utf8').includes('assets/js/nav.js'));
ok('every page with .mobile-toggle loads nav.js', unwired.length === 0,
  unwired.length + ' unwired: ' + unwired.slice(0, 5).map((p) => path.relative(ROOT, p)).join(', '));
ok('found a meaningful number of pages', withToggle.length >= 50, 'pages with toggle=' + withToggle.length);

/* pages that ship their own inline nav rule must also collapse it on mobile */
const inlineNav = withToggle.filter((p) => /nav\s*\{\s*display\s*:\s*flex/.test(fs.readFileSync(p, 'utf8')));
const inlineUncollapsed = inlineNav.filter((p) => {
  const s = fs.readFileSync(p, 'utf8');
  return !/@media\s*\(max-width:\s*768px\)\s*\{[^]*?nav\s*\{[^}]*display\s*:\s*none/.test(s);
});
ok('inline nav pages collapse on mobile', inlineUncollapsed.length === 0,
  inlineUncollapsed.map((p) => path.relative(ROOT, p)).join(', '));

if (fails === 0) {
  console.log('NAV REGRESSION: all checks passed (' + withToggle.length + ' pages, ' + inlineNav.length + ' with inline nav CSS)');
  process.exit(0);
}
console.log('NAV REGRESSION: ' + fails + ' failure(s)');
process.exit(1);
