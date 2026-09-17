# Validation summary figures

Generated as vector PDF and 300-dpi PNG, using the adjacent standalone base-R script.

1. **01_query_background_fraction:** query genes divided by the cohort's common coverage-testable background. Each point is a complete hypothesis, never a cycle. Boxes describe between-hypothesis distributions. P0 is highlighted. The denominator is explicit in `matched/strategy_comparison.tsv`.
2. **02_positions_before_after_filter:** nominally significant primary positions before dense-cluster pruning (light); surviving positions (dark overlay of identical width). These are positions, not genes. P1 and N6 have 15 linked cycles and should only be compared internally for discovery volume.
3. **03_conditional_benchmark:** only produced when a prespecified inferential summary exists. Points are plus-one upper-tail benchmark fractions. Intervals are exact binomial Monte Carlo exceedance intervals, not biological effect confidence intervals. Zero interval bounds are clipped to 1e-5 for the log axis. Multiple-testing adjustment, when prespecified, is in the table, not these unadjusted points.

Smoke figures are marked as synthetic software tests. They cannot establish biological performance.

Regenerate independently of Nextflow:

```bash
Rscript plot_validation.R ../matched/strategy_comparison.tsv . ../functional
```
