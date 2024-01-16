#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Run this script to convert from anna output gff to 

"""

import argparse
def read_bed(file_bed):
    """

    :param file_bed:
    :return:
    """
    data = []
    with open(file_bed) as file:
        for line in file.readlines():
            ls = line.split()
            if not line.startswith("#"): 
                if not line.startswith("Pan"):
                    if line.startswith("0"):
                        chr = ls[0].split("_")[1]
                        data.append(["TAIR10_" + chr, float(ls[1]), float(ls[2])])
                    elif line.startswith("220011"):
                        chr = ls[0].split("_")[1]
                        data.append(["22001_"+ chr + "_mod2", float(ls[1]), float(ls[2])])
                    else:
                        data.append([ls[0], float(ls[1]), float(ls[2])])
    return data

def write_newBED(data, file_output):
    """

    :param data: modified BED file
    :param file_output: output file name
    :return:
    """
    with open(file_output, "w") as file:
         for x in data:
             print("\t".join([x[0], str(int(x[1])), str(int(x[2]))]), file = file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input file", required=True)
    parser.add_argument("-o", "--out", help = "outfile", required=True)
    args = parser.parse_args()

    print("Read")
    data = read_bed(args.input)
    print("Write")
    write_newBED(data, args.out)