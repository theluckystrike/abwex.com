#!/usr/bin/env bash
# Deploy abwex.com and prove it actually landed.
#
# Why this exists. On 2026-08-07 a push went out and the HTML went live before
# GitHub Pages had propagated a brand new asset. Cloudflare fetched the asset,
# got a 404, and cached that 404 under "cache-control: max-age=14400". The page
# was live and broken for about three minutes and nothing reported it. The only
# reason it was caught was a manual curl.
#
# So this script does three things the bare "git push" did not.
#   1. Stamps every same-origin JS and CSS reference with a content hash, so an
#      edited asset always lands on a fresh cache key and a stale copy can never
#      be served for the four hour TTL.
#   2. Waits for the deploy and then fetches EVERY same-origin asset the changed
#      pages reference, failing loudly if any of them is not 200.
#   3. Tries a Cloudflare purge, and if the token cannot purge, says exactly what
#      to purge by hand instead of pretending it succeeded.
#
# Usage
#   ./deploy.sh            verify current live state only, no push
#   ./deploy.sh --push     commit any pending change, push, wait, verify, try purge
#   ./deploy.sh --stamp    rewrite every same-origin js/css ref to carry a content
#                          hash. Touches a lot of files, so it is deliberately not
#                          part of --push. Run it on purpose, review the diff.
set -uo pipefail

SITE="https://abwex.com"
ZONE="68dca5d93a5bf1a1bc21345e23e48a3f"
UA="Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 Chrome/126 Safari/537.36"
MAX_WAIT=420          # seconds to wait for a deploy to appear
POLL=15               # seconds between polls
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT" || exit 1

red()  { printf '\033[31m%s\033[0m\n' "$*"; }
grn()  { printf '\033[32m%s\033[0m\n' "$*"; }
ylw()  { printf '\033[33m%s\033[0m\n' "$*"; }

code_of() { curl -s -o /dev/null -w '%{http_code}' -A "$UA" "$1"; }

# --- 1. stamp same-origin js/css references with a content hash ---------------
stamp() {
  local changed=0
  while IFS= read -r page; do
    local before after
    before="$(cat "$page")"
    after="$before"
    while IFS= read -r ref; do
      [ -z "$ref" ] && continue
      local file="${ref%%\?*}"
      local disk=".${file}"
      [ -f "$disk" ] || continue
      local h
      h="$(shasum -a 256 "$disk" | cut -c1-8)"
      after="$(printf '%s' "$after" | sed -E "s|(${file})(\?v=[0-9a-f]+)?|\1?v=${h}|g")"
    done < <(printf '%s' "$before" | grep -oE '(src|href)="(/assets/[^"]+\.(js|css))' | sed -E 's/^(src|href)="//' | sort -u)
    if [ "$after" != "$before" ]; then
      printf '%s' "$after" > "$page"
      echo "  stamped $page"
      changed=1
    fi
  done < <(git ls-files '*.html')
  return $changed
}

# --- 2. collect every same-origin asset referenced by the tracked pages -------
collect_assets() {
  git ls-files '*.html' | while IFS= read -r page; do
    grep -oE '(src|href)="(/assets/[^"]+)"' "$page" 2>/dev/null | sed -E 's/^(src|href)="//; s/"$//'
  done | sort -u
}

verify_live() {
  local fails=0 n=0
  echo
  echo "Verifying live assets on $SITE"
  while IFS= read -r a; do
    [ -z "$a" ] && continue
    n=$((n+1))
    local c
    c="$(code_of "${SITE}${a}")"
    if [ "$c" != "200" ]; then
      red "  $c  ${a}"
      fails=$((fails+1))
    else
      printf '  200 %s\n' "$a"
    fi
  done < <(collect_assets)
  echo "  checked $n asset(s), $fails failure(s)"
  return $fails
}

purge() {
  local token="${CLOUDFLARE_API_TOKEN:-}"
  if [ -z "$token" ]; then
    ylw "No CLOUDFLARE_API_TOKEN set, skipping purge."
    return 1
  fi
  local body resp ok
  body="$(collect_assets | sed "s|^|\"${SITE}|; s|$|\",|" | tr -d '\n' | sed 's/,$//')"
  resp="$(curl -s -X POST "https://api.cloudflare.com/client/v4/zones/${ZONE}/purge_cache" \
    -H "Authorization: Bearer ${token}" -H "Content-Type: application/json" \
    --data "{\"files\":[${body}]}")"
  ok="$(printf '%s' "$resp" | python3 -c 'import json,sys;print(json.load(sys.stdin).get("success"))' 2>/dev/null)"
  if [ "$ok" = "True" ]; then
    grn "Cloudflare purge OK."
    return 0
  fi
  ylw "Cloudflare purge FAILED. The token lacks the Zone.Cache Purge permission."
  ylw "Either add that permission to the token, or purge these by hand at"
  ylw "  https://dash.cloudflare.com -> abwex.com -> Caching -> Configuration -> Purge Custom URLs"
  collect_assets | sed "s|^|  ${SITE}|"
  return 1
}

# --- main ---------------------------------------------------------------------
if [ "${1:-}" = "--stamp" ]; then
  echo "Stamping same-origin asset references with content hashes"
  stamp || true
  echo "Review with: git diff"
  exit 0
fi

if [ "${1:-}" = "--push" ]; then
  if ! git diff --quiet || [ -n "$(git status --porcelain)" ]; then
    git add -A
    git commit -q -m "${2:-Deploy pending changes}"
    echo "Committed $(git rev-parse --short HEAD)"
  fi
  git push origin main || { red "push failed"; exit 1; }
  grn "Pushed $(git rev-parse --short HEAD)"

  echo
  echo "Waiting for the deploy to appear (up to ${MAX_WAIT}s)"
  waited=0
  target="$(collect_assets | head -1)"
  while [ "$waited" -lt "$MAX_WAIT" ]; do
    if verify_live >/dev/null 2>&1; then break; fi
    sleep "$POLL"; waited=$((waited+POLL))
    printf '  %ss\n' "$waited"
  done
  purge || true
  sleep 5
fi

if verify_live; then
  grn "DEPLOY VERIFIED. Every referenced asset returns 200."
  exit 0
fi
red "DEPLOY NOT CLEAN. One or more referenced assets is not 200."
red "If an asset 404s but exists in the repo, Cloudflare has cached the 404."
red "Purge it (see above) or wait out cache-control: max-age=14400."
exit 1
