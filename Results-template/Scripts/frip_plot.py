#!/usr/bin/env python3

"""
Name: frip_plot.py
Created by: Tovah Markowitz
Date: 5/14/19

Purpose: To create visually appealing plots of FRiP scores
Currently only works with python/3.5
"""

##########################################
# Modules
import optparse
from pybedtools import BedTool
import pysam
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

##########################################
# Functions

def split_infiles(infiles):
    """ 
    breaks the infile string with space-delimited file names and 
    creates a list
    """
    infileList = infiles.strip("\'").strip('\"').split(" ")
    if len(infileList) == 1:
        infileList = infileList[0].split(";")
    return(infileList)

def count_reads_in_bed(bam, bedfile, genomefile):
    """
    some of this comes directly from the pybedtools site; read in 
    bed (or bed-like) file, sort it, and then count the number of 
    reads within the regions
    """
    bedinfo = BedTool(bedfile)
    bedinfo.sort(g=genomefile)
    return (
        BedTool(bam).intersect( bedinfo, bed=True, stream=True, )
        ).count()

def count_reads_in_bam(bam):
    """ count the number of reads in a given bam file """
    return( pysam.AlignmentFile(bam).mapped )

def calculate_frip(nreads, noverlaps):
    """ calculate FRiP score from nreads and noverlaps """
    return( float(noverlaps) / nreads )

def measure_bedfile_coverage(bedfile, genomefile):
    """ calculate the number of bases covered by a given bed file """
    bedinfo = BedTool(bedfile)
    return( bedinfo.sort(g=genomefile).total_coverage() )

def clip_bamfile_name(bamfile):
    """ 
    clip bam file name for table/plotting purposes; assumes file 
    naming system matches that of Pipeliner
    """
    sample = bamfile.split("/")[-1].split(".")[0]
    condition =  ".".join(bamfile.split("/")[-1].split(".")[1:-1])
    return( sample, condition )

def clip_bedfile_name(bedfile):
    """
    clip bed file name for table/plotting purposes; assumes file 
    naming system matches that of Pipeliner
    """
    toolused = bedfile.split("/")[-3]
    sample = bedfile.split("/")[-2]
    return( toolused, sample )

def process_files(bamfiles, bedfiles, genome):
    """ 
    this is the main function to take in list of input files and 
    put out an array containing key file name information, read 
    counts, and FRiP scores
    """
    bamfileL = split_infiles(bamfiles)
    bedfileL = split_infiles(bedfiles)
    out = [[ "bedtool", "bedsample", "bamsample", "bamcondition", 
    "n_reads", "n_overlap_reads", "FRiP", "n_basesM" ]]
    for bam in bamfileL:
        nreads = count_reads_in_bam(bam)
        (bamsample, condition) = clip_bamfile_name(bam)
        for bed in bedfileL:
            (bedtool, bedsample) = clip_bedfile_name(bed)
            noverlaps = count_reads_in_bed(bam, bed, genome)
            frip = calculate_frip(nreads, noverlaps)
            nbases = measure_bedfile_coverage(bed, genome) / 1000000
            out.append( [bedtool, bedsample, bamsample, condition, 
                nreads, noverlaps, frip, nbases] )
    out2 = pd.DataFrame(out[1:], columns=out[0])
    return(out2)

def create_outfile_names(outroot):
    """ uses outroot to create the output file names """
    outtable = "FRiP_table.txt"
    outbar = "FRiP_barplot.png"
    outscatter = "FRiP_scatterplot.png"
    if outroot != "":
        outtable = outroot + "." + outtable
        outbar = outroot + "." + outbar
        outscatter = outroot + "." + outscatter
    return(outtable, outbar, outscatter)

def write_table(out2, outtable):
    out2.to_csv(outtable,sep='\t',index=False)

def create_barplot(out2,outbar):
    """ 
    creates a barplot of FRiP scores where the x-axis is the bam
    sample, the y-axis is the FRiP score, the color is the sample 
    used to create the peak file, and the panels are the different 
    peak calling tools
    """
    sns.set(style="whitegrid", palette="Set1", font_scale=1.5)
    bp = sns.catplot(x="bamsample", y="FRiP", hue="bedsample",
             col="bedtool", data=out2, kind="bar", col_wrap=2)
    bp.set_axis_labels("Bam File", 'Fraction Reads in Peaks (FRiP)')
    bp.set_titles("{col_name}")
    bp._legend.set_title("Peak File")
    bp.set_xticklabels(rotation=10)
    #plt.show(bp)
    plt.savefig(outbar, bbox_inches='tight')
    plt.close("all")

def create_scatter(out2, outscatter):
    """
    Create a scatterplot of FRiP scores relative to number of bases
    under the peaks. Each panel is for a single bam file. Dots are
    colored based upon sample used for peak calling and given symbols
    based upon peak caller used.
    """
    bams= out2.loc[:,'bamsample'].unique()
    nplots= len(bams)
    sns.set(style="whitegrid", palette="Set1")
    f, axes = plt.subplots(1, nplots, sharey=True)
    for bi in range(nplots-1):
        tmp = out2.loc[ out2['bamsample'] == bams[bi] ]
        sns.scatterplot(data=tmp, x="n_basesM", y="FRiP",
                        hue="bedsample", style="bedtool", ax=axes[bi],
                        markers=['o','s','v','X'], legend=False)
        axes[bi].set(xlabel="# bases in peaks (M)", 
                     ylabel='Fraction Reads in Peaks (FRiP)')
        axes[bi].set_title( bams[bi] )
        axes[bi].get_xaxis().set_minor_locator( mpl.ticker.AutoMinorLocator() )
        axes[bi].grid(b=True, which='major', color="gray", linewidth=1)
        axes[bi].grid(b=True, which='minor', linestyle="--", 
                      color="gray", linewidth=0.5)
    tmp = out2.loc[ out2['bamsample'] == bams[nplots-1] ]
    sns.scatterplot(data=tmp, x="n_basesM", y="FRiP", hue="bedsample",
                    style="bedtool", ax=axes[nplots-1],
                    markers=['o','s','v','X'])
    axes[nplots-1].set(xlabel="# bases in peaks (M)", 
                       ylabel='Fraction Reads in Peaks (FRiP)')
    axes[nplots-1].set_title( bams[nplots-1] )
    axes[nplots-1].get_xaxis().set_minor_locator( mpl.ticker.AutoMinorLocator() )
    axes[nplots-1].grid(b=True, which='major', color="gray", linewidth=1)
    axes[nplots-1].grid(b=True, which='minor', linestyle="--",
                        color="gray", linewidth=0.5)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
    #plt.show()
    plt.savefig(outscatter, bbox_inches='tight')
    plt.close("all")

###############################################
# Main

def main():
    desc="""
This function takes a space-delimited or semi-colon delimited list
of bed-like files (extensions must be recognizable by bedtools)
and a list of bam files. It will then calculate the FRiP score for
all possible combinations of files and save the information in a
txt file and a barplot. It will also calculate the number of bases
covered by each bed-like file and create a scatterplot. Note: this
function assumes that the file naming system of the input files
matches that of Pipeliner.
    """

    parser = optparse.OptionParser(description=desc)

    parser.add_option('-p', dest='peakfiles', default='', 
           help='A space- or semicolon-delimited list of peakfiles \
(or bed-like files).')
    parser.add_option('-b', dest='bamfiles', default='', 
           help='A space- or semicolon-delimited list of bamfiles.')
    parser.add_option('-g', dest='genomefile', default='', 
           help='The name of the .genome file so bedtools knows the \
size of every chromosome.')
    parser.add_option('-o', dest='outroot', default='', 
           help='The root name of the multiple output files. Default: ""')

    (options,args) = parser.parse_args()
    bedfiles = options.peakfiles
    bamfiles = options.bamfiles
    genomefile = options.genomefile
    outroot = options.outroot

    out2 = process_files(bamfiles, bedfiles, genomefile)
    (outtable, outbar, outscatter) = create_outfile_names(outroot)
    write_table(out2, outtable)
    create_barplot(out2,outbar)
    create_scatter(out2, outscatter)

if __name__ == '__main__':
    main()

###############################################
# example cases

#bedfiles = "macs_broad/mWT_HCF1_mm_i81/mWT_HCF1_mm_i81_peaks.broadPeak macs_broad/mWT_HCF1_mm_i89/mWT_HCF1_mm_i89_peaks.broadPeak"
#bamfiles = "bam/Input_mm_i95.sorted.Q5DD.bam bam/mWT_HCF1_mm_i81.sorted.Q5DD.bam bam/mWT_HCF1_mm_i89.sorted.Q5DD.bam"
#genomefile = "/data/CCBR_Pipeliner/db/PipeDB/Indices/mm10_basic/indexes/mm10.fa.sizes"
#out2 = pd.read_csv("FRIP_test.txt",sep="\t")
