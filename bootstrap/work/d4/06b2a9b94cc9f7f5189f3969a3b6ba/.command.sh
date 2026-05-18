#!/bin/bash -ue
python "/Users/fabio/pCloud Drive/Bio/Projects/active.research/traits.evolution.omar/phyloq/bootstrap/caastools/ct" discovery \
    -a "ali3.fasta" \
    -t "b_5.txt" \
    -o "ali3.b_5.caas" \
    --fmt fasta

if [ ! -f "ali3.b_5.caas" ]; then
    touch "ali3.b_5.caas"
fi
