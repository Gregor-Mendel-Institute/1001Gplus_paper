#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

def read_stats(file_stats):
    """
    Read the stats file
    """
    data = []
    with open(file_stats) as file:
        for line in file.readlines():
            if not line.startswith("#"):
                ls = line.split()
                sizes = [int(x) for x in ls[2].split(",")]
                if max(sizes) >= 15:
                    data.append(ls[0])

    data = set(data)
    return data

def read_bed(file_bed, svs):
    """
    Read the bed file 
    
    """
    bed = []
    with open(file_bed) as file:
        for line in file.readlines():
            ls = line.split()
            if ls[3] in svs:
                bed.append(line)
    return bed




def write_file(file_output, data):
    """
    Write the output in a new file
    """
    with open(file_output, "w") as file:
        for x in data:
            file.write(x)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--stats", help="stats file", required=True)
    parser.add_argument("-b", "--bed", help="bed file", required=True)

    parser.add_argument("-o", "--output", help="output dataframe", required=True)
    args = parser.parse_args()

    data = read_stats(args.stats)
    bed = read_bed(args.bed, data)
    write_file(args.output, bed)
