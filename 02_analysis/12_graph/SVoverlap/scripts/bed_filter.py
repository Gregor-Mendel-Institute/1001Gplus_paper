#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse


def read_bed(file_bed):
    '''
    Read the bed file
    '''
    data = []
    with open(file_bed) as file:
        for line in file.readlines():
            data.append(line.split())
    return data

def filter(data):
    '''
    Filter bed file, only return SVs which are bigger than 14 bp
    '''
    data2 = []
    for x in data:
        if abs(int(float(x[1])) - int(float(x[2]))) > 14:
            data2.append([x[0], min(int(float(x[1])), int(float(x[2]))), max(int(float(x[1])), int(float(x[2])))])
    return data2

def write_data(data, file_output):
    '''
    Write data to a new bed file
    '''
    with open(file_output, "w") as file:
        for x in data:
            print("\t".join([str(y) for y in x]), file = file)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input file", required=True)
    parser.add_argument("-o", "--out", help = "outfile", required=True)
    args = parser.parse_args()

    print("Reading")
    data = read_bed(args.input)
    print("Filter")
    data2 = filter(data)

    print("Writing")
    write_data(data2, args.out)
