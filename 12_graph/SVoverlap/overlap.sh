#!/bin/bash

file1=$1
file2=$2
genome=$3

/ebio/abt6_projects9/abt6_software/bin/bedtools/bin/multiIntersectBed -empty -header -i $file1 $file2 -names graph, annagram -g $genome > overlap.txt