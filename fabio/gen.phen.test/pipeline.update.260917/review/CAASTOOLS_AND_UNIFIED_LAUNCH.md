# Modified CAAStools and unified launch verification

Date: 17 September 2026. This update follows Fabio's explicit request to launch
CAAS directly on the cluster without a mandatory pilot or staged approvals.

## The bundled tool is the previously used modified version

The 20 non-cache files in `bin/caastools` were compared with the previous
workflow's `../pipeline/bin/caastools`. The tree hashes match exactly and also
match the frozen design lock:

```text
e5ca882e4c5f7e85029b842a420aa31ee5dee63bf3bbf4060b60227564afc97f
```

The `ct` entry-point SHA-256 is:

```text
a665dc1f3128c2c14bca5861ebb9050159ca86ecf0c557bedf535f143c7ed971
```

The selected command in `scripts/run_caas.py` is `ct pooled-discovery`, using
the complete pool (`-t`) and existing frozen cycles (`-s`), and writing both
the legacy pooled table and the species-level event table. It does not download
another release, use a PATH-installed `ct`, or replace linked-pair cycles with
independently sampled endpoints. CAAStools itself was not modified for this
launcher update. The frozen tool hash is checked again during launch and tasks.

## Current direct interface

`launch_all_caas.sh` submits one SLURM driver, which creates a run-specific
inventory and explicit selection and starts one Nextflow `benchmark` run.
Defaults: N3/N4/N5, R0_001–099, R1_001–099, P1/N6 = 203 hypotheses / 20,130 cycles.
References P0/N0/N1/N2 and auxiliary P2 require explicit inclusion flags.

The direct stage requires no pilot report, approval JSON or reviewed inferential
metric. This is implemented as a separately recorded discovery-only mode, not
by fabricating approvals or setting `statistics.approved` to true. Existing
frozen inputs, cycle rules, filtering, tool checks, output separation and resume
identity checks remain enforced. The original staged interface remains optional.

Enrichment and gene-list/inferential null calculations are not requested during
this CAAS-only run. Output gene lists and descriptive summaries are retained for
later analysis. The functional status records `not_requested_discovery_only`.

## Verification

Seven focused tests cover the complete 203/208/102-hypothesis selections,
99-replicate inclusion, stable launch inputs and refusal of changed resume
inputs, direct discovery of replicate 99 without pilot approval, rejection of
altered frozen cycles, skipping enrichment, and the shell driver/plan paths.
SLURM submission was tested using a fake `sbatch`; no real cluster job was sent.

An actual local Nextflow run, `unified-discovery-check`, exercised nine fixture
hypotheses against two synthetic alignments (18 CAAS tasks), including the
new direct-launch validation, assembly/filtering and figure generation. Inputs
from `tests/fixtures` are explicitly recorded as synthetic. The output check
is saved to `review/tests/unified_discovery_smoke.json`.

This is bounded software verification, not an obligatory user pilot or a test
of full cluster throughput. Biological discovery has not been launched here.
The existing per-hypothesis/per-alignment scheduling granularity is unchanged;
203 hypotheses against 16,133 alignments means 3,274,999 CAAS tasks and requires
corresponding runtime and shared storage. See `RUN_ON_CLUSTER.md` for launching,
resource overrides and resume.
