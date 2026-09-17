"""Shared, dependency-free contracts for the isolated CAAS validation workflow."""
from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERSION = "caas-validation-20260917.1"
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def identifier(value):
    require(bool(SAFE_ID.fullmatch(value)), f"Unsafe identifier: {value!r}")
    return value


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            digest.update(block)
    return digest.hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def tree_hash(path):
    root = Path(path)
    return digest({str(p.relative_to(root)): sha(p) for p in sorted(root.rglob("*"))
                   if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"})


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames is not None, f"Missing TSV header: {path}")
        return list(reader)


def write_tsv(path, rows, columns):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in columns} for row in rows)


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def read_json(path):
    return json.loads(Path(path).read_text())


def resolve(base, value):
    p = Path(value)
    return p.resolve() if p.is_absolute() else (Path(base) / p).resolve()


def checked(path, checksum):
    require(Path(path).is_file(), f"Missing file: {path}")
    require(sha(path) == checksum, f"SHA-256 mismatch: {path}")
    return Path(path)


def pool_read(path):
    pools = {"FG": [], "BG": []}
    seen = set()
    for line in Path(path).read_text().splitlines():
        fields = line.split("\t")
        require(len(fields) == 2 and fields[1] in ("0", "1"), f"Invalid pool adapter: {path}: {line}")
        s, state = fields
        require(s not in seen, f"Duplicate/conflicting species {s} in {path}")
        seen.add(s)
        pools["FG" if state == "1" else "BG"].append(s)
    return pools


def cycle_read(path, pools):
    rows, keys, ids = [], set(), set()
    for line in Path(path).read_text().splitlines():
        fields = line.split("\t")
        require(len(fields) == 3, f"Invalid cycle adapter: {path}")
        cid, fg, bg = fields
        f, b = fg.split(","), bg.split(",")
        require(len(f) == len(set(f)) == len(b) == len(set(b)) == 4, f"Not unique 4 vs 4: {cid}")
        require(set(f) <= set(pools["FG"]) and set(b) <= set(pools["BG"]), f"Cycle outside fixed pools: {cid}")
        key = (tuple(sorted(f)), tuple(sorted(b)))
        require(key not in keys and cid not in ids, f"Duplicate cycle: {cid}")
        keys.add(key); ids.add(cid)
        rows.append((cid, f, b))
    require(bool(rows), f"No cycles: {path}")
    require(set().union(*(set(f) for _, f, _ in rows)) == set(pools["FG"]), "FG pool species never sampled")
    require(set().union(*(set(b) for _, _, b in rows)) == set(pools["BG"]), "BG pool species never sampled")
    return rows


def verify_design(design):
    design = Path(design).resolve()
    lock = read_json(design / "design.lock.json")
    for rel, checksum in lock["files"].items():
        checked(design / rel, checksum)
    return lock
