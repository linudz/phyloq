# Input paths inherited from the previous workflows

## Alignments: reuse the old CAAS collection

Authoritative configuration: `../caas/conf/cluster.config` in the `phyloq`
repository. It declares `${projectDir}/inputs/alignments/*.phy`, where
`projectDir` is the **old CAAS directory**, not this new sibling bundle.

Therefore `launch_all_caas.sh` now defaults to:

```text
phyloq/fabio/gen.phen.test/caas/inputs/alignments/*.phy
```

From the new pipeline this is `../caas/inputs/alignments/*.phy`. No absolute
cluster home/scratch path is hardcoded; moving the whole checkout remains safe.
The directory and matching files must exist on the cluster, not on the local
Mac. They remain outside Git. Explicit `--alignments-dir`,
`--alignments-pattern` or `--inventory` takes precedence.

For the original research-project layout, the old pipeline is named `pipeline`,
so the equivalent default is `../pipeline/inputs/alignments/*.phy`. If neither
legacy configuration exists, the standalone fallback is
`inputs/alignments/*.phy`. An existing legacy configuration with missing inputs
is an error, not permission to silently switch to a different collection.

## Trees: use the bundled diagnostic phylogeny, not RER gene trees

Both the previous and new CAAS pooled commands take an alignment, complete
species pools and cycle configurations. Neither takes a gene-tree argument.

The validation bundle already contains its checksum-locked diagnostic tree:

```text
inputs/frozen/brain-260917-v1/sources/tree.nwk
SHA-256: 4a9546635cceb3b8520f111c78eaffa7369eadb0818f73251e69539462854189
```

It was copied from the relative-brain-mass primate analysis:
`trait.in.primates/results/relative_brain_mass.Kuderna_S4.analysis_tree.nwk`.
The old absolute Mac path in `inputs/source.lock.json` is provenance, not an
execution requirement on the cluster. `scripts/pool_diagnostics.py` reads the
bundled tree, and launch validation checks its frozen checksum. The launch
summary records its resolved path/checksum; no tree is refitted or substituted.

The older `phyloq/bootstrap/nextflow.config` also names these **RERconverge** inputs:

```text
bootstrap/rerconverge/inputs/trees/rerconverge_gene_trees.tsv
bootstrap/rerconverge/inputs/science.abn7829_data_s4.nex.tree
```

Those paths are relative to that older bootstrap project. They are not CAAS
inputs and are deliberately not wired into this brain-trait validation.
