# Cercopithecidae absolute-trait tails

Approach 09 applies the same upper/lower-tail principle as benchmark 02 after
restricting the trait table to Cercopithecidae. The source contains 55 species.
Each tail contains `ceil(55 x 0.10) = 6` species:

- FG: the six highest relative-brain-mass values;
- BG: the six lowest relative-brain-mass values.

CAAStools cycles sample four species independently from each fixed pool. No
within-group genus restriction is imposed because the lower tail contains only
three genera and could not otherwise support a 4-vs-4 comparison. The 6-by-6
pools define 225 possible comparisons; 100 are selected reproducibly.

Regenerate the tail table and pooled inputs with:

```bash
python3 scripts/create_cercopithecidae_absolute_trait_tails.py
python3 scripts/create_pss_benchmark_configs.py
```
