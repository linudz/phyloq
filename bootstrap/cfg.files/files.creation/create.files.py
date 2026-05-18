#!/usr/bin/env python3

input_file = "BM_resampling.txt"

with open(input_file, "r") as f:
    lines = f.read().splitlines()

def line2file(line):
    columns = line.split("\t")
    zeros = [x + "\t" + "0" for x in columns[1].split(",")]
    unos = [x + "\t" + "1" for x in columns[2].split(",")]

    outstring = "\n".join(zeros) + "\n" + "\n".join(unos)
    return outstring

for line in lines[:10]:
    with open(line.split("\t")[0] + ".txt", "w") as f:
        f.write(line2file(line))
