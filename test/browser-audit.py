#!/usr/bin/env python3
"""Headless browser audit for abwex.com, driven over the Chrome DevTools Protocol.

Why CDP and not AppleScript. Driving the user's own Chrome turned out to be
unusable for QA: something on this machine closes newly created windows, and
once the QA window is gone an AppleScript `front window` reference silently
falls through to whatever window is left, so a sweep can resize, reload and
navigate the user's real tabs. This launches its own headless Chrome with its
own profile, talks to it directly, and never touches the visible browser.

Checks per page, all measured in a real layout engine rather than inferred:
  - console errors and page exceptions
  - real horizontal scroll at a phone width, by actually trying to pan
  - site shell present, correct fonts, correct background
  - em dashes in rendered text, which the static file gate cannot see
  - a screenshot when asked

Usage
  python3 test/browser-audit.py --base http://localhost:8899            all pages
  python3 test/browser-audit.py --base https://abwex.com --width 375
  python3 test/browser-audit.py --base http://localhost:8899 --shot tools/index.html
"""
import argparse, json, os, re, shutil, signal, socket, subprocess, sys, tempfile, time, base64
import urllib.request
import websocket  # websocket-client

CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def free_port():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    p = s.getsockname()[1]
    s.close()
    return p


class Chrome:
    def __init__(self, width, height):
        self.port = free_port()
        self.profile = tempfile.mkdtemp(prefix="abwex-cdp-")
        self.proc = subprocess.Popen(
            # --remote-allow-origins is required or Chrome answers the CDP websocket
            # handshake with 403, because websocket-client sends an Origin header.
            [CHROME, "--headless=new", "--disable-gpu", "--no-sandbox", "--no-first-run",
             "--disable-extensions", "--disable-background-networking", "--mute-audio",
             "--remote-allow-origins=*",
             f"--remote-debugging-port={self.port}", f"--user-data-dir={self.profile}",
             f"--window-size={width},{height}", "about:blank"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.ws = None
        self._id = 0
        for _ in range(100):
            try:
                tabs = json.load(urllib.request.urlopen(f"http://127.0.0.1:{self.port}/json", timeout=1))
                page = [t for t in tabs if t.get("type") == "page"]
                if page:
                    self.ws = websocket.create_connection(page[0]["webSocketDebuggerUrl"], timeout=30)
                    break
            except Exception:
                time.sleep(0.2)
        if not self.ws:
            raise RuntimeError("could not attach to headless Chrome")
        self.send("Page.enable")
        self.send("Runtime.enable")
        self.send("Log.enable")

    def send(self, method, **params):
        self._id += 1
        mid = self._id
        self.ws.send(json.dumps({"id": mid, "method": method, "params": params}))
        while True:
            msg = json.loads(self.ws.recv())
            if msg.get("id") == mid:
                return msg.get("result", {})

    # Errors thrown by third party tags are not site defects and we cannot fix
    # them. AdSense in particular throws a minified "Uncaught Y" on any page with
    # no fillable slot, which is every page in a headless run. Report them
    # separately so a real regression is never buried in that noise.
    THIRD_PARTY = ("googlesyndication", "adsbygoogle", "doubleclick", "google-analytics",
                   "googletagmanager", "gstatic", "fonts.googleapis", "favicon.ico")

    def drain(self):
        """Collect any buffered console/exception events without blocking long."""
        out = []
        self.ws.settimeout(0.25)
        try:
            while True:
                msg = json.loads(self.ws.recv())
                m = msg.get("method")
                if m == "Runtime.exceptionThrown":
                    d = msg["params"]["exceptionDetails"]
                    frames = (d.get("stackTrace") or {}).get("callFrames", [])
                    origin = " ".join([str(d.get("url") or "")] + [str(f.get("url") or "") for f in frames])
                    line = ("exception: " + (d.get("text") or "") + " " +
                            str((d.get("exception") or {}).get("description", ""))[:120])
                    out.append(("third-party " if any(t in origin for t in self.THIRD_PARTY) else "") + line)
                elif m == "Log.entryAdded":
                    e = msg["params"]["entry"]
                    if e.get("level") == "error":
                        origin = str(e.get("url") or "") + " " + str(e.get("text") or "")
                        line = "console: " + str(e.get("text"))[:120]
                        out.append(("third-party " if any(t in origin for t in self.THIRD_PARTY) else "") + line)
        except Exception:
            pass
        finally:
            self.ws.settimeout(30)
        return out

    def goto(self, url, settle=1.6):
        self.send("Page.navigate", url=url)
        time.sleep(settle)
        return self.drain()

    def js(self, expr):
        r = self.send("Runtime.evaluate", expression=expr, returnByValue=True, awaitPromise=True)
        if "exceptionDetails" in r:
            return {"__error": str(r["exceptionDetails"].get("text"))}
        return r.get("result", {}).get("value")

    def shot(self, path):
        r = self.send("Page.captureScreenshot", format="png", captureBeyondViewport=True)
        if r.get("data"):
            open(path, "wb").write(base64.b64decode(r["data"]))
            return True
        return False

    def close(self):
        try:
            self.ws.close()
        except Exception:
            pass
        self.proc.send_signal(signal.SIGTERM)
        try:
            self.proc.wait(timeout=5)
        except Exception:
            self.proc.kill()
        shutil.rmtree(self.profile, ignore_errors=True)


AUDIT_JS = r"""
(function(){
  var d=document.documentElement;
  window.scrollTo(0,0); var x0=window.scrollX;
  window.scrollTo(400,0); var x1=window.scrollX; window.scrollTo(0,0);
  var cs=getComputedStyle(document.body);
  var h1=document.querySelector('h1');
  var txt=document.body.innerText||'';
  return {
    title: document.title,
    vw: d.clientWidth,
    docScrollW: d.scrollWidth,
    pans: x1>0,
    panBy: Math.round(x1),
    bg: cs.backgroundColor,
    bodyFont: /IBM Plex/.test(cs.fontFamily),
    h1Font: h1 ? /Space Mono/.test(getComputedStyle(h1).fontFamily) : null,
    h1Count: document.querySelectorAll('h1').length,
    nav: !!document.querySelector('header .header-inner'),
    footer: !!document.querySelector('.site-footer'),
    emDash: (txt.match(/—/g)||[]).length,
    tables: document.querySelectorAll('table').length,
    tablesFill: Array.prototype.every.call(document.querySelectorAll('table'), function(t){
      var p=t.parentElement; if(!p) return true;
      return t.getBoundingClientRect().width >= p.clientWidth-2 || t.scrollWidth>p.clientWidth;
    })
  };
})()
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="http://localhost:8899")
    ap.add_argument("--width", type=int, default=375)
    ap.add_argument("--height", type=int, default=1400)
    ap.add_argument("--shot", default=None, help="capture a screenshot of this path")
    ap.add_argument("--only", default=None)
    a = ap.parse_args()

    pages = subprocess.run(["git", "ls-files", "*.html"], cwd=REPO,
                           capture_output=True, text=True).stdout.split()
    if a.only:
        pages = [p for p in pages if a.only in p]

    c = Chrome(a.width, a.height)
    fails, thirdparty, checked = [], [], 0
    try:
        if a.shot:
            c.goto(f"{a.base}/{a.shot}", settle=3)
            out = os.path.join(tempfile.gettempdir(), "abwex-shot.png")
            print("screenshot:", out if c.shot(out) else "FAILED")
            return 0
        for p in pages:
            errs = c.goto(f"{a.base}/{p}")
            r = c.js(AUDIT_JS)
            checked += 1
            if not isinstance(r, dict) or r.get("__error"):
                fails.append((p, f"audit js failed: {r}")); continue
            bad = []
            if r["pans"]:
                bad.append(f"pans sideways {r['panBy']}px (vw={r['vw']} scrollW={r['docScrollW']})")
            if r["emDash"]:
                bad.append(f"{r['emDash']} em dash(es) in rendered text")
            if r["h1Count"] != 1:
                bad.append(f"{r['h1Count']} h1 elements")
            if not r["bodyFont"]:
                bad.append("body not IBM Plex Sans")
            if r["h1Font"] is False:
                bad.append("h1 not Space Mono")
            if r["bg"] != "rgb(10, 10, 10)":
                bad.append("background " + r["bg"])
            if not r["nav"]:
                bad.append("no site nav")
            if not r["footer"]:
                bad.append("no site footer")
            for e in errs:
                if e.startswith("third-party "):
                    thirdparty.append((p, e))
                else:
                    bad.append(e)
            if bad:
                fails.append((p, "; ".join(bad)))
        print(f"pages audited: {checked} at {a.width}px, base {a.base}")
        if fails:
            if thirdparty:
                print(f"({len(thirdparty)} third-party error(s) ignored, AdSense and friends)")
            print(f"FAIL on {len(fails)} page(s):")
            for p, why in fails:
                print(f"  {p}\n      {why}")
            return 1
        print("all pages clean: no sideways scroll, no rendered em dashes, shell and type correct, "
              "one h1 each, no first-party console errors")
        if thirdparty:
            print(f"({len(thirdparty)} third-party error(s) ignored, AdSense and friends)")
        return 0
    finally:
        c.close()


if __name__ == "__main__":
    sys.exit(main())
