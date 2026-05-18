#!/bin/bash -ue
python "/Users/fabio/pCloud Drive/Bio/Projects/active.research/traits.evolution.omar/phyloq/bootstrap/caastools/ct" discovery \
    -a "ali1.fasta" \
    -t "b_6.txt" \
    -o "ali1.b_6.caas" \
    --fmt fasta

if [ ! -f "ali1.b_6.caas" ]; then
    touch "ali1.b_6.caas"
fi
