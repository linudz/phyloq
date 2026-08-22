# PSS benchmark with pooled CAAStools

This workflow benchmarks eight strategies for defining genome–phenome
hypotheses from relative brain mass. Here, *benchmark* refers to the direct
comparison of these predefined hypothesis-building strategies under the
same pooled CAAStools settings.

## Biological hypothesis

The hypothesis under test is that the Pairwise Shift Score (PSS) helps define
better foreground/background comparisons for genome–phenome analyses. If PSS
captures informative phenotypic shifts, PSS-informed species pairs should
yield stronger, more coherent, or more reproducible CAAStools results than
groups based only on taxonomic extrema or absolute trait tails.

This stage generates discovery results. The downstream benchmark must compare
the approaches with prespecified result-level metrics; raw numbers of
positive sites alone should not be interpreted as proof that one strategy is
biologically superior.

## Common experimental design

- Eight independent approaches: three original benchmarks and five additional
  PSS-informed pooling hypotheses.
- 100 unique cycles per approach.
- Four foreground (FG) versus four background (BG) species per cycle.
- A fixed seed (`260821`) and SHA-256 ranking make selection reproducible
  across Python versions and operating systems.
- Higher phenotype is encoded as FG; lower phenotype is encoded as BG.
- Species have fixed side membership within each approach, as required for
  pooled event reconstruction.
- The eight configurations are crossed with all alignments by Nextflow and
  can therefore run in parallel.

Each pooled configuration is headerless and uses the format:

```text
cycle_id<TAB>FG_species_separated_by_commas<TAB>BG_species_separated_by_commas
```

## Approach 1: family minima and maxima

Source: `inputs/config.creation/01_family_trait_extrema.tsv`.

Each cycle samples four families as paired units. For every selected family,
the species with the maximum relative brain mass enters FG and the species
with the minimum value enters BG. This preserves the family-level contrast in
every cycle.

Families represented by only one species cannot define a directional
contrast: the same species would be both the minimum and maximum and would
therefore occur in FG and BG. `Aotidae` (`Aotus_trivirgatus`) and
`Daubentoniidae` (`Daubentonia_madagascariensis`) are excluded for this reason.

## Approach 2: absolute trait tails

Source: `inputs/config.creation/02_trait_distribution_tails.tsv`.

FG species are sampled from the high absolute-trait tail already marked `FG`;
BG species are sampled from the low tail marked `BG`. Four species are chosen
independently for each side. Within FG and within BG, no two selected species
may belong to the same genus. The stated restriction is within each group; it
does not impose a cross-side genus-matching rule.

## Approach 3: PSS-informed genus pairs

Source: `inputs/config.creation/04_best_top1pct_pair_per_genus.tsv`.

The source contains the highest-PSS pair from each represented genus, provided
that the pair belongs to the global top 1% of the PSS test. Each cycle samples
four complete pairs. Within each pair, the species with the higher absolute
phenotype enters FG and the species with the lower phenotype enters BG. Pair
membership is never broken during sampling.

This is the PSS-informed arm of the benchmark.

## Approach 4: all congeneric pairs in the global PSS top 1%

Source: `inputs/config.creation/05_all_global_top1pct_pss_pairs.tsv`.

Every cycle samples four complete congeneric pairs from the global upper 1%
of PSS. The higher-phenotype endpoint of each pair enters FG and the lower
endpoint enters BG. Four pairs are retained only when they provide four
distinct FG and four distinct BG species.

*Macaca leonina* and *M. maura* occur on both sides across the source pairs.
All pairs involving either species are excluded because pooled event
reconstruction requires fixed side membership. This leaves 39 eligible source
pairs, 47,402 valid cycles, and a selected-cycle discovery pool of 24 FG and
22 BG species.

## Approach 5: top-1% pairs within Macaca, Papio, and Trachypithecus

Source: `inputs/config.creation/05_all_global_top1pct_pss_pairs.tsv`.

Both endpoints must belong to *Macaca*, *Papio*, or *Trachypithecus*; the
filter does not require them to belong to the same genus. In the observed
global top 1%, however, all 32 qualifying pairs happen to be congeneric. After
the same fixed-side exclusion, 26 source pairs yield 3,425 valid cycles and a
selected-cycle discovery pool of 14 FG and 9 BG species.

## Approach 6: top-1% pairs within Macaca and Papio

Source: `inputs/config.creation/05_all_global_top1pct_pss_pairs.tsv`.

This focal Papionini subset includes only the two prespecified genera
*Macaca* and *Papio*; it is not intended to represent every genus in the
tribe. After the fixed-side exclusion, 19 source pairs yield 371 valid cycles
and a selected-cycle discovery pool of 8 FG and 7 BG species.

## Approach 7: hierarchical PSS matching, 13 versus 13

Source: `inputs/config.creation/08_pss_ranked_endpoint_disjoint_13_pairs.tsv`.

The global top-1% PSS pairs are traversed from highest to lowest score. A pair
is retained only if neither endpoint has already been selected. The endpoint
with higher relative brain mass enters FG and its linked lower endpoint enters
BG. Selection stops at 13 pairs, producing 13 unique FG and 13 unique BG
species without cross-group ambiguities. With the current data, rank 20 is
sufficient to recover all 13 endpoint-disjoint pairs.

The 13 pairs define 715 possible linked 4-vs-4 cycles. The common seeded
selection retains 100 cycles. FG and BG are never sampled independently: each
cycle is always assembled from four intact PSS pairs.

This design tests paired directional shifts. It does not enforce separation
between the complete FG and BG trait distributions: an FG endpoint from one
pair may have a lower trait value than a BG endpoint from another pair.

Recreate the ranked source tables before regenerating all configurations with:

```bash
python3 scripts/create_pss_ranked_13x13_groups.py
python3 scripts/create_pss_benchmark_configs.py
```

## Approach 8: `pss.cercopithecidae.random.pools`

Sources:

- `inputs/config.creation/09_cercopithecidae_pss_bottom_linked_pairs.tsv`;
- `inputs/config.creation/09_cercopithecidae_pss_random_pools.tsv`.

The PSS analysis is restricted to Cercopithecidae. A top-1% pair contributes
to pool definition when at least one endpoint also occurs in a bottom-1% pair.
The higher endpoint defines FG membership and the lower endpoint defines BG
membership. Deduplication yields nine FG and ten BG species without cross-side
ambiguities.

The pairs justify the two pools but do not constrain individual comparisons.
Each cycle independently samples four FG and four BG species, with no genus or
linked-pair restriction. This produces 26,460 possible comparisons; 100 are
selected reproducibly with the common seed.

## Planned benchmark: `pss.with.separation`

This future benchmark will combine top-1% PSS support with complete separation
of the final trait distributions. The current maximum endpoint-disjoint
solution is 12 FG versus 12 BG around a relative-brain-mass threshold of
approximately 1.4565. It is distinct from approach 7 and has no pooled config
yet.

## Recreate the configurations

From the `CAAS` directory, run:

```bash
python3 scripts/create_pss_benchmark_configs.py
```

Useful options are:

```text
--cycles 100       number of cycles per approach
--group-size 4     number of FG and BG species per cycle
--seed 260821      deterministic selection seed
--examples 3       examples printed for each approach
```

The script validates group sizes, duplicate species, fixed FG/BG membership,
and the within-group genus restriction for the absolute-tail strategy. It
prints example cycles to the terminal and writes SHA-256 checksums to the
manifest.

## Generated inputs

The eight pooled configurations are:

```text
inputs/benchmark-configs/01_family_extrema.pooled.caas.cfg
inputs/benchmark-configs/02_absolute_trait_tails.pooled.caas.cfg
inputs/benchmark-configs/03_pss_best_pair_per_genus.pooled.caas.cfg
inputs/benchmark-configs/04_pss_all_congeneric_top1pct.pooled.caas.cfg
inputs/benchmark-configs/05_pss_macaca_papio_trachypithecus_top1pct.pooled.caas.cfg
inputs/benchmark-configs/06_pss_macaca_papio_top1pct.pooled.caas.cfg
inputs/benchmark-configs/07_pss_ranked_endpoint_disjoint_13x13.pooled.caas.cfg
inputs/benchmark-configs/08_pss_cercopithecidae_random_pools.pooled.caas.cfg
```

CAAStools also requires a complete fixed-side pool for reconstructing coherent
events and defining event denominators. These supporting files are stored in:

```text
inputs/benchmark-pools/
```

`inputs/benchmark.configs.tsv` links each approach to its source table, pooled
configuration and complete pool. It also records the number of possible and
selected cycles and the configuration checksum. Paths in the manifest are
relative to this project, so the directory can be transferred to the cluster.

`inputs/benchmark.configs.pss-ranked-13x13.tsv` contains only approach 7 and
can be used to launch this strategy without resubmitting the previous arms.

`inputs/benchmark.configs.pss-cercopithecidae-random-pools.tsv` contains only
approach 8 for a focused launch of the random-pool benchmark.

## Example selections

With the default seed, the first selected cycle of each strategy is:

```text
01 family extrema
FG: Propithecus_edwardsi, Nycticebus_coucang, Galago_senegalensis, Lepilemur_ruficaudatus
BG: Avahi_laniger, Loris_tardigradus, Otolemur_garnettii, Lepilemur_dorsalis

02 absolute trait tails
FG: Ateles_paniscus, Macaca_nigra, Pongo_abelii, Mandrillus_sphinx
BG: Cheirogaleus_medius, Otolemur_crassicaudatus, Loris_tardigradus, Perodicticus_potto

03 PSS best pair per genus
FG: Trachypithecus_cristatus, Semnopithecus_entellus, Colobus_polykomos, Cercopithecus_cephus
BG: Trachypithecus_auratus, Semnopithecus_priam, Colobus_guereza, Cercopithecus_petaurista
```

## Pooled CAAStools analysis

For every alignment and approach, CAAStools evaluates all 100 cycles. A
position is retained when at least one cycle is positive. Positive cycles are
then combined into amino-acid-compatible events: species are counted once per
event, while incompatible amino-acid signatures remain separate events.

The positional hypergeometric prefilter defaults to `0.05`. Each 4-vs-4 cycle
requires at least three observed, non-gap species on each side. These settings
are shared across the eight benchmark arms in `conf/cluster.config`.

## Nextflow execution

Alignment files are installed directly on Correfoc under
`inputs/alignments/` and are excluded from Git. The configured input glob is
`inputs/alignments/*.phy`. Set the remaining SLURM settings in
`conf/cluster.config`, then launch:

```bash
sbatch submit_pipeline_slurm.sh
```

Nextflow reads `inputs/benchmark.configs.tsv`, forms the Cartesian product of
the eight approaches and all alignments, and submits the resulting jobs in
parallel subject to the configured queue limit. Use a fresh run for this new
graph; subsequently, the same run can be resumed with:

```bash
sbatch submit_pipeline_slurm.sh -resume
```

For the first test of the hierarchical 13-vs-13 strategy, launch only its
one-row manifest as a fresh run:

```bash
sbatch submit_pipeline_slurm.sh \
  --benchmark_manifest "$PWD/inputs/benchmark.configs.pss-ranked-13x13.tsv"
```

To launch only the cercopithecid random-pool strategy:

```bash
sbatch submit_pipeline_slurm.sh \
  --benchmark_manifest "$PWD/inputs/benchmark.configs.pss-cercopithecidae-random-pools.tsv"
```

The original three hypothesis and pool files are byte-identical to those used
in the first run. On resume, their gene-level tasks retain the same cache keys;
Nextflow therefore reuses them and schedules the three newly added approaches.

Results are separated by approach:

```text
results/RUN_ID/APPROACH/caas-pooled/GENE.pooled.caas.tsv
results/RUN_ID/APPROACH/caas-pooled-events/GENE.pooled.caas.events.tsv
```

The event table is the primary biological output. The legacy position table
retains positive-cycle counts and identifiers for traceability; overlapping
cycles are not statistically independent replicates.
