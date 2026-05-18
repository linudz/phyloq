#!/bin/bash -ue
python "/Users/fabio/pCloud Drive/Bio/Projects/active.research/traits.evolution.omar/phyloq/bootstrap/caastools/ct" discovery \
    -a "ali2.fasta" \
    -t "b_3.txt" \
    -o "ali2.b_3.caas" \
    --fmt fasta

if [ ! -f "ali2.b_3.caas" ]; then
    touch "ali2.b_3.caas"
fi
