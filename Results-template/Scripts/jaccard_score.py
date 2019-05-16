#!/usr/bin/env python3

"""
Name: jaccard_score.py
Created by: Tovah Markowitz
Date: 1/23/19

Purpose: To do all pairwise comparisons of bed/peak files given. Uses bedtools
to calculate a jaccard score for every comparison. All data is saved in a 
single tab-delimited file.
"""

##########################################
# Modules
import optparse
from pybedtools import BedTool

##########################################
# Functions

def split_infiles(infiles):
    # breaks the infile string with space-delimited file names and creates a list
    infileList = infiles.strip("\'").strip('\"').split(" ")
    if len(infileList) == 1:
        infileList = infileList[0].split(";")
    return(infileList)

def loop_jaccard(infileList, genomefile):
    """ Running bedtools. Uses two loops to do all possible pairwise comparisons of files 
    in a list. Reads in a bedtools approved file types, sorts the files, and calculates a 
    jaccard score.
    """
    nfiles = len(infileList)
    out = []
    for z in range(nfiles):
        fileA = infileList[z]
        print("fileA is: " + fileA) 
        a = BedTool(fileA)
        a = a.sort(g=genomefile)
        for y in range(z+1,nfiles):
            fileB = infileList[y]
            b = BedTool(fileB)
            b = b.sort(g=genomefile)
            j = a.jaccard(b,g=genomefile)
            j["fileA"] = fileA.split("/")[-1]
            j["fileB"] = fileB.split("/")[-1]
            if len(out) == 0:
                keylist = list(j.keys())
                keylist.sort()
                out.append( "\t".join(keylist) )
            tmp = [ str(j[key]) for key in keylist ]
            out.append( "\t".join(tmp) )
    return(out)

def write_out(out, outFile):
    f = open(outFile, 'w')
    f.write( "\n".join(out) )
    f.close()

##########################################
# Main

def main():
    desc="""
    This function takes a space-delimited list of files (bed, bedgraph, gff, gtf, etc.)
    and calculates all possible pairwise jaccard scores. From bedtools: 'Jaccard is the 
    length of the intersection over the union. Values range from 0 (no intersection) to 
    1 (self intersection)'. The columns of the output file are: fileA, fileB, 
    intersection, jaccard, n_intersections, and union-intersection.
    """

    parser = optparse.OptionParser(description=desc)

    parser.add_option('-i', dest='infiles', default='', help='A space- or semicolon-delimited list of \
input files for analysis.')
    parser.add_option('-o', dest='outfile', default='', help='The name of the output file \
where all the jaccard score information will be saved.')
    parser.add_option('-g', dest='genomefile', default='', help='The name of the .genome file.')

    (options,args) = parser.parse_args()
    infiles = options.infiles
    outfile = options.outfile
    genomefile = options.genomefile

    infileList = split_infiles(infiles)
    out = loop_jaccard(infileList, genomefile)
    write_out(out, outfile)

if __name__ == '__main__':
    main()


