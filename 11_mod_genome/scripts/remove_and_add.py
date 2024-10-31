#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/23/22

@author: moinSebi

"""
import argparse



def read_fasta(filename):
    """ 
    Read the faste file
    """
    data = []
    seq = ""
    header = ""
    with open(filename) as file:
        for line in file.readlines():
            if line.startswith(">"):
                if len(header) != 0:
                    data.append([header, seq])
                seq = ""
                header = line.replace("\n", "")
            else:
                seq += line.replace("\n", "")

    data.append([header, seq])
    return data


def read_bed(bedfile):
    """
    Read the bed file
    """
    data = dict()
    with open(bedfile) as file:
        for line in file.readlines():
            line = line.split()
            data[line[0]] =  line[1:]
    return data

def remove_add(fasta, bed):
    """
    Remove parts of a genome and add new one
    """
    print("djasd")
    seq = ""
    seq2 = ""
    for x in fasta:
        if x[0][1:] == "22001_Chr3":
            seq = x[1][:int(bed["22001_Chr3"][1])]
            x[1] = x[1][int(bed["22001_Chr3"][1]):]
        if x[0][1:] == "22001_Chr5":
            seq2 = x[1][int(bed["22001_Chr5"][0]):]
            x[1] = x[1][:int(bed["22001_Chr5"][0])]

    for x in fasta:
        # This just checks
        if x[0] == ">22001_Chr5":
            x[1] = x[1] + rev_comp_st(seq)
        if x[0] == ">22001_Chr3":
            x[1] = rev_comp_st(seq2) + x[1]



def rev_comp_st(seq):
    '''
    This function returns a reverse complement
    of a DNA or RNA strand
    '''
        # complement strand
    seq = seq.replace("A", "t").replace(
        "C", "g").replace("T", "a").replace("G", "c")
    seq = seq.upper()

    # reverse strand
    seq = seq[::-1]
    return seq

def write_fasta(fasta, filename):
    """
    Write the fasta file
    """
    with open(filename, "w") as file:
        for x in fasta:
            print(x[0].split("_")[0] + "f_" + x[0].split("_")[1], file = file)
            print(x[1], file = file)



if __name__ == "__main__":
    """
    Main function
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", help="fasta file", required=True)
    parser.add_argument("-b", "--bed", help = "bed file with coordinates", required=True)
    parser.add_argument("-o", "--output", help = "output file name", required = True)
    args = parser.parse_args()

    fasta = read_fasta(args.fasta)
    bed = read_bed(args.bed)
    remove_add(fasta, bed)
    write_fasta(fasta, args.output)

