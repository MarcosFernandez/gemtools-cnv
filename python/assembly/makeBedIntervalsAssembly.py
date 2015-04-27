#!/usr/bin/env python
# makeBedIntervals.py Creates a bed of chromosome intervals

import sys, argparse
import os.path

def is_valid_file(file):
    """Checks if file exists"""
    if not os.path.exists(file):
       return False
    else:
       return True

def is_valid_path(path):
    """Check if paths exists"""
    if not os.path.exists(path):
        return False
    else:
        return True

#PARSING ARGUMENTS
parser = argparse.ArgumentParser(description='Creates a bed of chromosome intervals.')
parser.add_argument('-chrName', metavar='chrName', help='chromosome name',
                    required=True)
parser.add_argument('-chrLen', metavar='chrLen',type=int, help='chromosome length',
                    required=True)
parser.add_argument('-outFile', metavar='outFile', help='Path to store bed file',
                    required=True)

args = parser.parse_args()

#SIZE OF kmers AND step
kmerLen = 36
step = 5;

#OPEN BED FILE RESULTING
fileOutput = args.outFile
with open(fileOutput, 'w') as bedFile:
    begin = 1
    end = begin + kmerLen - 1
    while end <= args.chrLen:
        bedFile.write(args.chrName + "\t" + str(begin-1) + "\t" + str(end) + "\n");
        begin += step
        end = begin + kmerLen - 1 
