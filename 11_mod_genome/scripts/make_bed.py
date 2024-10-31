#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/23/22

@author: moinSebi

"""
import argparse


# Read Paf and check for smallest and biggest alignment
def read_file(filename, m1, m2):
    data = 10000000000000
    p = 0
    with open(filename) as file:
        for line in file.readlines():
            ds = line.split()
            mi = min(int(ds[2]), int(ds[3]))
            ma = max(int(ds[2]), int(ds[3]))
            if mi < data:
                data = mi
            if ma > p:
                p = ma

    if m1:
        print(ds[0], 0, data, sep = "\t")
    if m2:
        print(ds[0], p, 10000000000000, sep = "\t")




if __name__ == "__main__":
    """
    Main function
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--paf", help="paf file", required=True)
    parser.add_argument("--min", help = "min", action = "store_true", default=False)
    parser.add_argument("--max", help = "max", action = "store_true", default=False)
    args = parser.parse_args()

    o = read_file(args.paf, args.min, args.max)

