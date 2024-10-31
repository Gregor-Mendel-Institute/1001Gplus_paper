#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys

def read_wrapper(dir, chr):
    """
    Iterate in this folder over each file which starts with *chr*
    """
    data = dict()
    names = []
    for path, currentDirectory, files in os.walk(dir):
        for file in files:
            if file.startswith(chr):
                read_vcf(path + "/" + file, data, names)
    return data, names


def read_vcf(filename, data, k):
    """
    Reading the vcf file 
    """
    print("Filename", filename, file = sys.stderr)
    print("Size of the data", len(data), file = sys.stderr)
    with open(filename) as file:

        # Parsing the name of the file
        name = filename.split(".")[-2]
        k.append(name)
        for line in file.readlines():
            if not line.startswith("#"):
                ls = line.split()
                
                # These are the anchors
                start_stop = ls[2]


                if start_stop not in data:
                    pa = ls[9:]
                    allele = [ls[3]] + ls[4].split(",")
                    bo = check_allele(allele)
                    allels = ls[7].split(";")[3][3:].split(",")
                    ssizes = [len(ls[3])] + [len(x) for x in ls[4].split(",")]

                    if bo:
                        # Adjust sizes (reduce by 1)
                        s2 = [x-1 for x in ssizes]
                        data[start_stop] = (name, allels, pa, s2)
                    else:
                        data[start_stop] = (name, allels, pa, ssizes)
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



def clean_data(data, names):
    """
    We add the allele information for the missing stuff
    """
    data2 = dict()
    names = sorted(names)
    for k,v in data.items():
        # What is the index for the reference this
        index = names.index(v[0])

        # This is the old traversal
        p = v[2]
        # Add 
        p.insert(index, 0)
        # Remove the name, add the sizes, and the new traversal
        data2[k] = (v[1], v[3], p)
    return data2, names

def write_stats(data2, names, file_out):
    """
    Write statistics to a file
    """
    with open(file_out, "w") as file:
        print("#" + ",".join(names), file = file)
        for k,v in data2.items():
            print(k + "\t" + ",".join(v[0]) + "\t" + ",".join([str(y) for y in v[1]]) + "\t" + "\t".join([str(y) for y in v[2]]), file = file)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="gfa file", required=True)
    parser.add_argument("-c", "--chromosome", help="chromosome", required=True)
    parser.add_argument("-o", "--output", help="output dataframe")

    args = parser.parse_args()

    data, names = read_wrapper(args.directory, args.chromosome)
    data2, names2 = clean_data(data, names)
    write_stats(data2, names2, args.output);
