from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path

import requests


def download(url: str, output: Path) -> None:
    session = requests.Session()
    session.headers.update({"User-Agent": "LegionellaHostDependencyAtlas/1.0"})
    response = session.get(url, timeout=60)
    ctype = response.headers.get("content-type", "").lower()
    if "text/html" in ctype or response.content.lstrip().startswith(b"<html"):
        text = response.text
        match = re.search(r'POW_CHALLENGE\s*=\s*"([^"]+)"', text)
        diff = re.search(r'POW_DIFFICULTY\s*=\s*"(\d+)"', text)
        name = re.search(r'POW_COOKIE_NAME\s*=\s*"([^"]+)"', text)
        if not (match and diff and name):
            raise RuntimeError("PMC download returned HTML without a recognized proof-of-work challenge")
        challenge, difficulty, cookie_name = match.group(1), int(diff.group(1)), name.group(1)
        prefix = "0" * difficulty
        nonce = 0
        while not hashlib.sha256((challenge + str(nonce)).encode()).hexdigest().startswith(prefix):
            nonce += 1
        session.cookies.set(cookie_name, f"{challenge},{nonce}", path="/")
        response = session.get(url, timeout=120)
        response.raise_for_status()
    if response.content.lstrip().startswith(b"<html"):
        raise RuntimeError("PMC proof-of-work cookie was not accepted")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(response.content)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("url")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    download(args.url, args.output)
    print(f"downloaded {args.output} ({args.output.stat().st_size} bytes)")
