#!/bin/bash -ue
python "/Users/fabio/pCloud Drive/Bio/Projects/active.research/traits.evolution.omar/phyloq/bootstrap/caastools/ct" discovery \
    -a "ali2.fasta" \
    -t "b_9.txt" \
    -o "ali2.b_9.caas" \
    --fmt fasta

if [ ! -f "ali2.b_9.caas" ]; then
    touch "ali2.b_9.caas"
fi
