#!/usr/bin/env python3

import sys
import getopt
import re
import os
import statistics
import random
import csv
import pysam
import numpy

class RefNotFasta(Exception):
    pass

class SeqNotInReference(Exception):
    pass

class IncorrectCigarChar(Exception):
    pass

def main(args):
    RefFile, NewFile, Coverage, Target = None, None, None, None
    seed = 0
    try:
        opts, args = getopt.getopt(args, "r:n:c:t:s:h", ["reference=", "new=", "coverage=", "target=", "seed=", "help"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    if opts == []:
        usage()
        sys.exit(2)
    
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-r", "--reference"):
            RefFile = arg
        elif opt in ("-n", "--new"):
            NewFile = arg
        elif opt in ("-c", "--coverage"):
            Coverage = arg
        elif opt in ("-t", "--target"):
            Target = int(arg)
        elif opt in ("-s", "--seed"):
            seed = int(arg)

    Ranges, Depths = readDepth(Coverage, Target)
    Keep, Header, Depths = GetReads(Ranges, Depths, RefFile, Target, seed)
    WriteFile(Header, NewFile, Keep)
    exit(0)

def usage():
    print("DownsampleToCoverage.py -r Reference.bam -n Out.bam -c Coverage.txt -t 1000 -s 12345678")
    exit(2)

def readDepth(Coverage, Target):
    ranges = []
    depths = []
    name = ""
    with open(Coverage, "r") as tsvin:
        tsvin = csv.reader(tsvin, delimiter = "\t")
        start, stop, prev = 0, 0, 0
        for row in tsvin:
            if row[2] != "coverage":
                coverage = int(row[2])
                position = int(row[1])
                name = row[0]
                depths.append(coverage)
                if coverage > Target:
                    if start == 0:
                        start = position
                        prev = position
                    else:
                        if position <= (prev + 3):
                            prev = position
                            continue
                        elif position > (prev + 3):
                            stop = prev
                            ranges.append([name, start, stop])
                            start = position
                            prev = position
                            stop = 0
        ranges.append([name, start, prev])
    [print(range) for range in ranges]
    return(ranges, numpy.array(depths))


def GetReads(Ranges, Depths, RefFile, Target, Seed):
    if Seed > 0:
        random.seed(Seed)
    keep = []
    with pysam.AlignmentFile(RefFile, "rb") as infile:
        print(infile.count())
        prev_end = 1
        if Ranges[0][2] > 0:
            for overSeq in Ranges:
                if overSeq[1] > prev_end:
                    end = overSeq[1] - 1
                    keep += list(infile.fetch(overSeq[0], prev_end, end))
                    prev_end = overSeq[2] + 1
                temp, Depths = Downsample(overSeq, Depths, Target, infile)
                keep += temp
                print(len(keep))
        if prev_end < infile.get_reference_length(Ranges[0][0]):
            end = infile.get_reference_length(Ranges[0][0])
            keep += list(infile.fetch(Ranges[0][0], prev_end, end))
        header = infile.header
        print(len(keep))
        return(list(set(keep)), infile.header, Depths)


def Downsample(Range, Depths, Target, InFile):
    reads = list(InFile.fetch(Range[0], Range[1], Range[2]))
    keep = []
    random.shuffle(reads)
    for read in reads:
        start = read.reference_start
        end = read.reference_end
        if read.cigarstring != None:
            if (start > Range[1]) and (end < Range[2]):
                number = random.uniform(0, 1)
                localDepth = statistics.mean(Depths[start:end])
                if number < (float(Target)/float(localDepth)):
                    keep.append(read)
                #else:
                    #Depths[start:end] -= 1
    return(keep, Depths)


def WriteFile(Header, NewFile, Keep):
    print("Writing %d reads." % (len(Keep),))
    with pysam.AlignmentFile(NewFile, "wb", header = Header) as outfile:
        [outfile.write(read) for read in Keep]
    print("reads written")


if __name__ == "__main__":
    main(sys.argv[1:])
    