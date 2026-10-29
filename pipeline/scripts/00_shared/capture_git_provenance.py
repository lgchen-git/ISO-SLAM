#!/usr/bin/env python3
"""Capture deterministic Git provenance for one pipeline stage."""

from __future__ import annotations

import argparse
import hashlib
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlsplit, urlunsplit


RELEASE_TAG = re.compile(r"^pipeline-v\d+\.\d+\.\d+$")


def git(repo: Path, *args: str, text: bool = True) -> str | bytes:
    result = subprocess.run(
        ["git", "-C", str(repo), *args],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=text,
    )
    if result.returncode != 0:
        message = result.stderr.strip() if text else result.stderr.decode(errors="replace").strip()
        raise RuntimeError(message or f"git {' '.join(args)} failed")
    return result.stdout


def sanitize_remote(value: str) -> str:
    value = value.strip()
    if not value:
        return value
    if "://" not in value:
        return value.split("@", 1)[1] if "@" in value else value
    parsed = urlsplit(value)
    hostname = parsed.hostname or ""
    if parsed.port:
        hostname = f"{hostname}:{parsed.port}"
    return urlunsplit((parsed.scheme, hostname, parsed.path, parsed.query, parsed.fragment))


def sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--dirty-policy", choices=("allow", "warn", "error"), default="warn"
    )
    args = parser.parse_args()

    repo = args.repo_root.resolve()
    try:
        inside = str(git(repo, "rev-parse", "--is-inside-work-tree")).strip()
        if inside != "true":
            raise RuntimeError("not inside a Git worktree")
        commit = str(git(repo, "rev-parse", "--verify", "HEAD")).strip()
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    try:
        branch = str(git(repo, "symbolic-ref", "--quiet", "--short", "HEAD")).strip()
    except RuntimeError:
        branch = "DETACHED"

    status_text = str(git(repo, "status", "--porcelain=v1", "--untracked-files=all"))
    dirty = bool(status_text)
    if dirty and args.dirty_policy == "error":
        print("ERROR: Git worktree is dirty", file=sys.stderr)
        return 3
    if dirty and args.dirty_policy == "warn":
        print("WARNING: Git worktree is dirty; provenance records this state", file=sys.stderr)

    diff = git(repo, "diff", "--binary", "HEAD", "--", text=False)
    describe = str(git(repo, "describe", "--always", "--dirty", "--tags")).strip()
    try:
        remote = sanitize_remote(str(git(repo, "remote", "get-url", "origin")))
    except RuntimeError:
        remote = ""
    tags = str(git(repo, "tag", "--points-at", "HEAD")).splitlines()
    release_tag = next((tag for tag in sorted(tags) if RELEASE_TAG.fullmatch(tag)), "")

    rows = (
        ("schema_version", "1"),
        ("captured_at_utc", datetime.now(timezone.utc).isoformat()),
        ("repo_name", repo.name),
        ("commit", commit),
        ("branch", branch or "DETACHED"),
        ("describe", describe),
        ("dirty", str(dirty).lower()),
        ("tracked_diff_sha256", sha256(diff)),
        ("status_sha256", sha256(status_text.encode("utf-8"))),
        ("changed_path_count", str(len(status_text.splitlines()))),
        ("remote_url", remote),
        ("pipeline_release_tag", release_tag),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "key\tvalue\n" + "".join(f"{key}\t{value}\n" for key, value in rows),
        encoding="utf-8",
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
