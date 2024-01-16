#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

def read_vcf(filename):
    """
    Reading VCF file. 

    Adjusted to this: https://www.biostars.org/p/84686/
    vcf 1-based, bed 0-based
    """
    data = []
    with open(filename) as file:
        for line in file.readlines():
            if not line.startswith("#"):
                ls = line.split()
                start = int(ls[1])

                allele = [ls[3]] + ls[4].split(",")

                bo = check_allele(allele)

                # length of the ref allele
                l1 = len(ls[3])

                # Get the allele (node order) with start and stop
                trav = ls[7].split(";")[3][3:].split(",")[0]

                # Data is [chr_name, start, length, ref_nodes, alt_nodes ]
                if bo:
                    data.append([ls[0], start, start + l1 - 1, ls[2], trav])
                else:
                    data.append([ls[0], start-1, start + l1 -1 , ls[2], trav])

    return data


def check_allele(l):
    """
    Check if the first character of the ref and alle are the same (this leads to different size)
    """
    oo = l[0][0]
    for x in l:
        if x[0] != oo:
            return False
    return True


def write_vcf(data, filename):
    """
    Write a BED file
    """
    with open(filename, "w") as file:
        for x in data:
            print("\t".join([str(y) for y in x]), file = file)

if __name__ == "__main__":
    """
    Main function
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="gfa file", required=True)
    parser.add_argument("-o", "--output", help="output dataframe", required=True)
    args = parser.parse_args()

    data = read_vcf(args.input)
    write_vcf(data, args.output)