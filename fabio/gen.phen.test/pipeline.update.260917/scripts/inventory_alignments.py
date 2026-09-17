#!/usr/bin/env python3
"""Create a hashed inventory; reject filename-to-gene collisions before discovery."""
import argparse
import glob
import re
from pathlib import Path
from common import ROOT, require, sha, write_tsv


def load_alignment(path):
    from Bio import AlignIO
    with Path(path).open() as handle:
        header = handle.readline().split()
    require(len(header) == 2 and all(x.isdigit() for x in header), f"Malformed PHYLIP header: {path}")
    n_species, n_positions = map(int, header)
    require(n_species > 0 and n_positions > 0, f"Empty PHYLIP alignment: {path}")
    msa = AlignIO.read(str(path), 'phylip-relaxed')
    require(len(msa) == n_species and msa.get_alignment_length() == n_positions,
            f"PHYLIP header/content dimensions disagree: {path}")
    require(len({r.id for r in msa}) == len(msa), f"Duplicate species in {path}")
    return msa


def inventory(pattern, output):
    paths = sorted(Path(p).resolve() for p in glob.glob(pattern) if Path(p).is_file())
    require(paths, f"No alignment files found: {pattern}")
    genes, basenames, rows = set(), set(), []
    for path in paths:
        # alimport.slice in the frozen CAAStools uses the first dot-separated token.
        gene = path.name.split(".", 1)[0]
        require(bool(re.fullmatch(r"[A-Za-z0-9_-]+", gene)), f"Invalid gene symbol: {path}")
        require(gene not in genes, f"Duplicate/ambiguous gene mapping: {gene}: {path}")
        require(path.name not in basenames, f"Duplicate alignment basename: {path.name}")
        genes.add(gene); basenames.add(path.name)
        alignment = load_alignment(path)
        species = [r.id for r in alignment]
        require(len(species) == len(set(species)), f"Duplicate species in {path}")
        require(alignment.get_alignment_length() > 0, f"Empty alignment: {path}")
        rows.append(dict(gene_id=gene, alignment_path=str(path), sha256=sha(path),
                         format="phylip-relaxed", n_species=len(species),
                         n_positions=alignment.get_alignment_length(), species=",".join(sorted(species))))
    write_tsv(output, rows, ["gene_id", "alignment_path", "sha256", "format", "n_species", "n_positions", "species"])
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--alignments", required=True, help="Quoted glob, never an implicit production default")
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    print(f"Inventoried {len(inventory(a.alignments, a.output))} alignments")
