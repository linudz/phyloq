#!/bin/bash -ue
python "/Users/fabio/pCloud Drive/Bio/Projects/active.research/traits.evolution.omar/phyloq/bootstrap/caastools/ct" discovery \
    -a "ali3.fasta" \
    -t "b_3.txt" \
    -o "ali3.b_3.caas" \
    --fmt fasta

if [ ! -f "ali3.b_3.caas" ]; then
    touch "ali3.b_3.caas"
fi
