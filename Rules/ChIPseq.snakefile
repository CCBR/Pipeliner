from snakemake.utils import R
from os.path import join
import re,os

from os import listdir

configfile: "run.json"
    
workpath = config['project']['workpath']    
filetype = config['project']['filetype']
readtype = config['project']['readtype']

se=""
pe=""

if config['project']['nends'] == 2 :
    pe="yes"
elif config['project']['nends'] == 1 :
    se="yes"

def peaks_per_chrom(file, chrom):
    """takes a peak file as input and counts how many peaks there are 
    on the chromosome of interest"""
    f = open(file, 'r')
    datain = f.readlines()
    f.close()
    data = [row.strip().split('\t')[0] for row in datain]
    return(data.count(chrom))

def outputIDR(groupswreps, groupdata, chip2input, tools):
    """
    Produces the correct output files for IDR. All supposed replicates
    should be directly compared when possible using IDR. IDR malfunctions
    with bed files and GEM so it will not run with either of those.
    Because there is no q-value calculated for SICER when there is no 
    input file, those samples are also ignored.
    """
    IDRgroup, IDRsample1, IDRsample2, IDRpeaktool = [], [], [], []
    for group in groupswreps:
        nsamples = len(groupdata[group])
        for i in range(nsamples):
            ctrlTF = chip2input[groupdata[group][i]] != ""
            for j in range(i+1,nsamples):
                if ctrlTF == (chip2input[groupdata[group][j]] != ""):
                    if ctrlTF == False:
                        tooltmp = [ tool for tool in tools if tool != "sicer" ]
                    else:
                        tooltmp = tools			           
                    IDRgroup.extend([group] * len(tooltmp))
                    IDRsample1.extend([groupdata[group][i]] * len(tooltmp))
                    IDRsample2.extend([groupdata[group][j]] * len(tooltmp))
                    IDRpeaktool.extend(tooltmp)
    return( IDRgroup, IDRsample1, IDRsample2, IDRpeaktool )

def zip_peak_files(chips, PeakTools, PeakExtensions):
    """ making input file names for FRiP """
    zipSample, zipTool, zipExt = [], [], []
    for chip in chips:
        for PeakTool in PeakTools:
            zipSample.append(chip)
            zipTool.append(PeakTool)
            zipExt.append(PeakExtensions[PeakTool])
    return( zipSample, zipTool, zipExt )

def calc_effective_genome_fraction(effectivesize, genomefile):
    """
    calculate the effective genome fraction by calculating the
    actual genome size from a .genome-like file and then dividing
    the effective genome size by that number
    """
    lines=list(map(lambda x:x.strip().split("\t"),open(genomefile).readlines()))
    genomelen=0
    for chrom,l in lines:
        if not "_" in chrom and chrom!="chrX" and chrom!="chrM" and chrom!="chrY":
            genomelen+=int(l)
    return( str( float(effectivesize)/ genomelen ) )

effectivegenomesize=config['references'][pfamily]['EFFECTIVEGENOMESIZE']
reflen = config['references'][pfamily]['REFLEN']

def zip_contrasts(contrast, PeakTools):
    """making output file names for differential binding analyses"""
    zipGroup1, zipGroup2, zipTool, contrasts = [], [], [], []
    for g1, g2 in contrast:
        for PeakTool in PeakTools:
            zipGroup1.append(g1)
            zipGroup2.append(g2)
            zipTool.append(PeakTool)
            contrasts.append( g1 + "_vs_" + g2 + "-" + PeakTool )
    return( zipGroup1, zipGroup2, zipTool, contrasts )

##########
# DEFINING SAMPLES

chips = config['project']['peaks']['chips']
chip2input = config['project']['peaks']['inputs']
uniq_inputs = list(sorted(set([v for v in chip2input.values() if v])))

sampleswinput = []

for input in chip2input:
        if chip2input[input] != 'NA' and chip2input[input] != '':
                sampleswinput.append(input)

groupdata = config['project']['groups']

groupdatawinput = {}
groupswreps = []

for group, chipsamples in groupdata.items() :
    tmp = [ ]
    if len(chipsamples) > 1:
        groupswreps.append(group)
    for chip in chipsamples :
        if chip in samples:
            tmp.append(chip)
            input = chip2input[chip]
            if input != 'NA' and input != '':
                tmp.append(input)
    if len(tmp) != 0:
        groupdatawinput[group]=set(tmp)

groups = list(groupdatawinput.keys())

reps=""
if len(groupswreps) > 0:
    reps="yes"

contrast = config['project']['contrast']

#########
# PREPARING TO DEAL WITH A VARIED SET OF PEAKCALL TOOLS

PeakTools_narrow = [ "macs_narrow", "gem" ]
PeakTools_broad = [ "macs_broad", "sicer" ]

PeakTools = PeakTools_narrow + PeakTools_broad
PeakToolsNG = [ tool for tool in PeakTools if tool != "gem" ]

PeakExtensions = { 'macs_narrow': '_peaks.narrowPeak', 'macs_broad': '_peaks.broadPeak',
                   'sicer': '_broadpeaks.bed', 'gem': '.GEM_events.narrowPeak' ,
                   'MANorm': '_all_MA.bed'}

FileTypesDiffBind = { 'macs_narrow': 'narrowPeak', 'macs_broad': 'narrowPeak',
                    'sicer': 'bed', 'gem': 'narrowPeak' }

PeakExtensionsIDR = { 'macs_narrow': '_peaks.narrowPeak', 'macs_broad': '_peaks.broadPeak',
                      'sicer': '_sicer.broadPeak' }

FileTypesIDR = { 'macs_narrow': 'narrowPeak', 'macs_broad': 'broadPeak',
                 'sicer': 'broadPeak' }

RankColIDR = { 'macs_narrow': 'q.value', 'macs_broad': 'q.value',
               'sicer': 'q.value' }

IDRgroup, IDRsample1, IDRsample2, IDRpeaktool =	outputIDR(groupswreps, groupdata, chip2input, PeakToolsNG)

zipSample, zipTool, zipExt = zip_peak_files(chips, PeakTools, PeakExtensions)
zipGroup1, zipGroup2, zipTool, contrasts = zip_contrasts(contrast, PeakTools)

#########
# CREATING DIRECTORIES

bam_dir='bam'
qc_dir='PeakQC'

idr_dir = 'IDR'
memechip_dir = "MEME"
homer_dir = "HOMER_motifs"
homer2_dir = "HOMER_annotations"
uropa_dir = "UROPA_annotations"
#diffbind_dir = "DiffBind"
manorm_dir = "MANorm"

otherDirs = [qc_dir, homer_dir, homer2_dir, uropa_dir]
if reps == "yes":
    otherDirs.append(idr_dir)
#    otherDirs.append(diffbind_dir)

for d in PeakTools + otherDirs:
        if not os.path.exists(join(workpath,d)):
                os.mkdir(join(workpath,d))

##########
# RULE ALL

if reps == "yes":
    rule ChIPseq:
        input:
            expand(join(workpath,"macs_narrow","{name}","{name}_peaks.narrowPeak"),name=chips),
            expand(join(workpath,"macs_broad","{name}","{name}_peaks.broadPeak"),name=chips),
            expand(join(workpath,"sicer","{name}","{name}_broadpeaks.bed"),name=chips),
            expand(join(workpath,"gem","{name}","{name}.GEM_events.narrowPeak"),name=chips),
            join(workpath,"FRiP_barplot.png"),
            expand(join(workpath,qc_dir,'{PeakTool}_jaccard.txt'),PeakTool=PeakTools),
            expand(join(workpath,homer_dir,'{PeakTool}',"{name}_{PeakTool}"),PeakTool=PeakToolsNG,name=chips),
            expand(join(workpath,homer2_dir,'{PeakTool}',"{name}_{PeakTool}_annotations.txt"),PeakTool=PeakTools_narrow,name=chips),
            expand(join(workpath, uropa_dir,'{PeakTool}','{name}_{PeakTool}_uropa_allhits.txt'),PeakTool=PeakTools,name=chips),
            expand(join(workpath,idr_dir,'{PeakTool}','{group}','{sample1}_vs_{sample2}.idrValue.txt'),zip,PeakTool=IDRpeaktool,group=IDRgroup,sample1=IDRsample1,sample2=IDRsample2),
#            expand(join(workpath,diffbind_dir,"{PeakTool}_diffbind_prep.csv"),PeakTool=PeakTools),
else:
    rule ChIPseq:
        input:
            expand(join(workpath,"macs_narrow","{name}","{name}_peaks.narrowPeak"),name=chips),
            expand(join(workpath,"macs_broad","{name}","{name}_peaks.broadPeak"),name=chips),
            expand(join(workpath,"sicer","{name}","{name}_broadpeaks.bed"),name=chips),
            expand(join(workpath,"gem","{name}","{name}.GEM_events.narrowPeak"),name=chips),
            expand(join(workpath,qc_dir,'{PeakTool}_jaccard.txt'),PeakTool=PeakTools),
            expand(join(workpath,homer_dir,'{PeakTool}',"{name}_{PeakTool}"),PeakTool=PeakToolsNG,name=chips),
            expand(join(workpath,homer2_dir,'{PeakTool}',"{name}_{PeakTool}_annotations.txt"),PeakTool=PeakTools_narrow,name=chips),
            expand(join(workpath, uropa_dir,'{PeakTool}','{name}_{PeakTool}_uropa_allhits.txt'),PeakTool=PeakTools,name=chips),
            expand(join(workpath, uropa_dir,'{PeakTool}','{name}_{PeakTool}_uropa_allhits.txt'),PeakTool="MANorm",name=contrasts),
            expand(join(workpath,manorm_dir,"{group1}_vs_{group2}-{tool}","{group1}_vs_{group2}-{tool}_all_MAvalues.xls"),zip,group1=zipGroup1,group2=zipGroup2,tool=zipTool),


##########
# INDIVIDUAL RULES

if se == "yes":
    rule MACS2_narrow:
        input:
            chip = join(workpath,bam_dir,"{name}.sorted.Q5DD.tagAlign.gz"),
            ppqt = join(workpath,bam_dir,"{name}.sorted.Q5DD.ppqt")
        output:
            join(workpath,"macs_narrow","{name}","{name}_peaks.narrowPeak"),
        params:
            rname='pl:MACS2_narrow',
            gsize=config['references'][pfamily]['EFFECTIVEGENOMESIZE'],
            macsver=config['bin'][pfamily]['tool_versions']['MACSVER'],
            macsn_dir="macs_narrow",
            ctrl = lambda w : join(workpath,bam_dir,chip2input[w.name] + ".sorted.Q5DD.tagAlign.gz"),
        run:
            commoncmd = "module load {params.macsver}; "
            file=list(map(lambda z:z.strip().split(),open(input.ppqt,'r').readlines()))
            extenders = []
            for ppqt_value in file[0][2].split(","):
                if int(ppqt_value) > 1:
                    extenders.append(ppqt_value)
            try:
                extsize = extenders[0]
            except IndexError:
                extsize = "{} {}".format(file[0][2].split(",")[0], "# Negative Value which will cause pipeline to fail (wrong ref genome selected or low starting DNA)")               
            if params.ctrl != join(workpath,bam_dir,".sorted.Q5DD.tagAlign.gz"):
                cmd = "macs2 callpeak -t " + input.chip + " -c " + params.ctrl + " -g " + params.gsize + " -n " + wildcards.name + " --outdir " + join(workpath,params.macsn_dir,wildcards.name) + " -q 0.01 --keep-dup='all' --nomodel --extsize " + extsize
            else:
                cmd = "macs2 callpeak -t " + input.chip + " -g " + params.gsize + " -n " + wildcards.name + " --outdir " + join(workpath,params.macsn_dir,wildcards.name) + " -q 0.01 --keep-dup='all' --nomodel --extsize " + extsize
            shell(commoncmd+cmd)

if pe == "yes":
    rule MACS2_narrow:
        input:
            chip = join(workpath,bam_dir,"{name}.sorted.Q5DD.bam"),
        output:
            join(workpath,"macs_narrow","{name}","{name}_peaks.narrowPeak"),
        params:
            rname='pl:MACS2_narrow',
            gsize=config['references'][pfamily]['EFFECTIVEGENOMESIZE'],
            macsver=config['bin'][pfamily]['tool_versions']['MACSVER'],
            macsn_dir="macs_narrow",
            ctrl = lambda w : join(workpath,bam_dir,chip2input[w.name] + ".sorted.Q5DD.bam"),
        shell: """
module load {params.macsver};
if [ {params.ctrl} != "{workpath}/{bam_dir}/.sorted.Q5DD.bam" ]; then
    macs2 callpeak -t {input.chip} -c {params.ctrl} -g {params.gsize} -n {wildcards.name} \
          --outdir {workpath}/{params.macsn_dir}/{wildcards.name} -q 0.01 --keep-dup="all" -f "BAMPE";
else
    macs2 callpeak -t {input.chip} -g {params.gsize} -n {wildcards.name} \
          --outdir {workpath}/{params.macsn_dir}/{wildcards.name} -q 0.01 --keep-dup="all" -f "BAMPE";
fi
"""

if se == "yes":
    rule MACS2_broad:
        input:
            chip = join(workpath,bam_dir,"{name}.sorted.Q5DD.tagAlign.gz"),
            ppqt = join(workpath,bam_dir,"{name}.sorted.Q5DD.ppqt")
        output:
            join(workpath,"macs_broad","{name}","{name}_peaks.broadPeak"),
        params:
            rname='pl:MACS2_broad',
            gsize=config['references'][pfamily]['EFFECTIVEGENOMESIZE'],
            macsver=config['bin'][pfamily]['tool_versions']['MACSVER'],
            macsb_dir="macs_broad",
            ctrl = lambda w : join(workpath,bam_dir,chip2input[w.name] + ".sorted.Q5DD.tagAlign.gz"),
        run:
            commoncmd = "module load {params.macsver}; "
            file=list(map(lambda z:z.strip().split(),open(input.ppqt,'r').readlines()))
            extenders = []
            for ppqt_value in file[0][2].split(","):
                if int(ppqt_value) > 1:
                    extenders.append(ppqt_value)
            try:
                extsize = extenders[0]
            except IndexError:
                extsize = "{} {}".format(file[0][2].split(",")[0], "# Negative Value which will cause pipeline to fail (wrong ref genome selected or low starting DNA)")               
            if params.ctrl != join(workpath,bam_dir,".sorted.Q5DD.tagAlign.gz"):
                cmd = "macs2 callpeak -t " + input.chip + " -c " + params.ctrl + " -g " + params.gsize + " -n " + wildcards.name + " --outdir " + join(workpath,params.macsb_dir,wildcards.name) + " --broad --broad-cutoff 0.01 --keep-dup='all' --nomodel --extsize " + extsize
            else:
                cmd = "macs2 callpeak -t " + input.chip + " -g " + params.gsize + " -n " + wildcards.name + " --outdir " + join(workpath,params.macsb_dir,wildcards.name) + " --broad --broad-cutoff 0.01 --keep-dup='all' --nomodel --extsize " + extsize
            shell(commoncmd+cmd)

if pe == "yes":
    rule MACS2_broad:
        input:
            chip = join(workpath,bam_dir,"{name}.sorted.Q5DD.bam"),
        output:
            join(workpath,"macs_broad","{name}","{name}_peaks.broadPeak"),
        params:
            rname='pl:MACS2_broad',
            gsize=config['references'][pfamily]['EFFECTIVEGENOMESIZE'],
            macsver=config['bin'][pfamily]['tool_versions']['MACSVER'],
            macsb_dir="macs_broad",
            ctrl = lambda w : join(workpath,bam_dir,chip2input[w.name] + ".sorted.Q5DD.bam"),
        shell: """
module load {params.macsver};
if [ {params.ctrl} != "{workpath}/{bam_dir}/.sorted.Q5DD.bam" ]; then
    macs2 callpeak -t {input.chip} -c {params.ctrl} -g {params.gsize} -n {wildcards.name} \
          --outdir {workpath}/{params.macsb_dir}/{wildcards.name} --broad --broad-cutoff 0.01 \
          --keep-dup="all" -f "BAMPE";
else
    macs2 callpeak -t {input.chip} -g {params.gsize} -n {wildcards.name} \
          --outdir {workpath}/{params.macsb_dir}/{wildcards.name} --broad --broad-cutoff 0.01 \
          --keep-dup="all" -f "BAMPE";
fi
"""

if se == "yes":
    rule SICER:
        input:
            chip = join(workpath,bam_dir,"{name}.sorted.Q5DD.tagAlign.gz"),
            ppqt = join(workpath,bam_dir,"{name}.sorted.Q5DD.ppqt"),
        output:
            txt = join(workpath,"sicer","{name}","{name}_broadpeaks.txt"),
# output columns: chrom, start, end, ChIP tag count, control tag count, p-value, fold-enrichment, q-value
        params:
            rname='pl:SICER',
            sicerver=config['bin'][pfamily]['tool_versions']['SICERVER'],
            bedtoolsver=config['bin'][pfamily]['tool_versions']['BEDTOOLSVER'],
            genomever = config['project']['annotation'],
            frac = calc_effective_genome_fraction(effectivegenomesize, reflen),
            ctrl = lambda w : join(workpath,bam_dir,chip2input[w.name] + ".sorted.Q5DD.tagAlign.gz"),
        run:
            commoncmd1 = "if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID; fi "
            commoncmd2 = "cd /lscratch/$SLURM_JOBID; "
            commoncmd3 = "module load {params.sicerver}; module load {params.bedtoolsver}; "
            cmd1 = "cp {input.chip} chip.bed.gz; gzip -d chip.bed.gz; "
            file=list(map(lambda z:z.strip().split(),open(input.ppqt,'r').readlines()))
            extenders = []
            for ppqt_value in file[0][2].split(","):
                if int(ppqt_value) > 1:
                    extenders.append(ppqt_value)
            try:
                extsize = extenders[0]
            except IndexError:
                extsize = "{} {}".format(file[0][2].split(",")[0], "# Negative Value which will cause pipeline to fail (wrong ref genome selected or low starting DNA)")               
            if params.ctrl != join(workpath,bam_dir,".sorted.Q5DD.tagAlign.gz"):
                cmd2 = "cp {params.ctrl} input.bed.gz; gzip -d input.bed.gz; "
                cmd3 =  "sh $SICERDIR/SICER.sh . chip.bed input.bed . {params.genomever} 100 300 " + extsize + " {params.frac} 600 1E-2 ; mv chip-W300-G600-islands-summary-FDR1E-2 {output.txt}"
                shell(commoncmd1)
                shell(commoncmd2 + commoncmd3 + cmd1 + cmd2 + cmd3)
            else:
                cmd2 = "sh $SICERDIR/SICER-rb.sh . chip.bed . {params.genomever} 100 300 " + extsize + " 0.75 600 100 ; mv chip-W300-G600-E100.scoreisland {output.txt}"
                shell(commoncmd1)
                shell(commoncmd2 + commoncmd3 + cmd1 + cmd2)

if pe =="yes":
    rule SICER:
        input: 
            chip = join(workpath,bam_dir,"{name}.sorted.Q5DD.bam"),
            ppqt = join(workpath,bam_dir,"{name}.sorted.Q5DD.ppqt"),
        output:
            txt = join(workpath,"sicer","{name}","{name}_broadpeaks.txt"),
# output columns: chrom, start, end, ChIP tag count, control tag count, p-value, fold-enrichment, q-value
        params:
            rname='pl:SICER',
            sicerver=config['bin'][pfamily]['tool_versions']['SICERVER'],
            bedtoolsver=config['bin'][pfamily]['tool_versions']['BEDTOOLSVER'],
            genomever = config['project']['annotation'],
            ctrl = lambda w : join(workpath,bam_dir,chip2input[w.name] + ".sorted.Q5DD.bam"),
        run:
            commoncmd1 = "if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID; fi "
            commoncmd2 = "cd /lscratch/$SLURM_JOBID; "
            commoncmd3 = "module load {params.sicerver}; module load {params.bedtoolsver}; "
            cmd1 = "bamToBed -i {input.chip} > chip.bed; "
            file=list(map(lambda z:z.strip().split(),open(input.ppqt,'r').readlines()))
            extenders = []
            for ppqt_value in file[0][2].split(","):
                if int(ppqt_value) > 1:
                    extenders.append(ppqt_value)
            try:
                extsize = extenders[0]
            except IndexError:
                extsize = "{} {}".format(file[0][2].split(",")[0], "# Negative Value which will cause pipeline to fail (wrong ref genome selected or low starting DNA)")               
            if params.ctrl != join(workpath,bam_dir,".sorted.Q5DD.bam"):
                cmd2 = "bamToBed -i {params.ctrl} > input.bed; "
                cmd3 =  "sh $SICERDIR/SICER.sh . chip.bed input.bed . {params.genomever} 100 300 " + extsize + " 0.75 600 1E-2 ; mv chip-W300-G600-islands-summary-FDR1E-2 {output.txt}"
                shell(commoncmd1)
                shell(commoncmd2 + commoncmd3 + cmd1 + cmd2 + cmd3)
            else:
                cmd2 = "sh $SICERDIR/SICER-rb.sh . chip.bed . {params.genomever} 100 300 " + extsize + " 0.75 600 100 ; mv chip-W300-G600-E100.scoreisland {output.txt}"
                shell(commoncmd1)
                shell(commoncmd2 + commoncmd3 + cmd1 + cmd2)

rule convertSICER:
    input:
         txt = join(workpath,"sicer","{name}","{name}_broadpeaks.txt"),
# input columns if input-normalized: chrom, start, end, ChIP tag count, control tag count, p-value, fold-enrichment, q-value
# input columns if no input: chrom, start, end, score
    output:
         bed = join(workpath,"sicer","{name}","{name}_broadpeaks.bed"),
# output broadPeak columns: chrom, start, end, name, ChIP tag count, strand, fold-enrichment, -log10 p-value, -log10 q-value
    params:
         rname='pl:convertSICER',
    run:
        import math
        f = open(input.txt,'r')
        intxt = f.readlines()
        f.close()
        outBroadPeak = [None] * len(intxt)
        outBed = [None] * len(intxt)
        for i in range(len(intxt)):
            tmp = intxt[i].strip().split('\t')
            if len(tmp) == 8:
# assuming that a p-value/q-value of 0 is super significant, -log10(1e-500)
                if tmp[5] == "0.0":
                    pval="500"
                else:
                    pval = str(-(math.log10(float(tmp[5]))))
                if tmp[7] == "0.0":
                    qval="500"
                    qvalScore="5000"
                else:
                    qval = str(-(math.log10(float(tmp[7]))))
                    qvalScore = str(int(-10*math.log10(float(tmp[7]))))
                outBroadPeak[i] = "\t".join(tmp[0:3] + ["Peak"+str(i+1),tmp[3],".", tmp[6], pval, qval])
                outBed[i] = "\t".join(tmp[0:3] + ["Peak"+str(i+1),qvalScore])
            else:
                score = str(int(float(tmp[3])))
                outBed[i] = "\t".join(tmp[0:3] + ["Peak"+str(i+1),score])
        g = open(output.bed,'w')
        g.write( "\n".join(outBed) )
        g.close()
        if outBroadPeak[0] != None:
            h = open( join(workpath, "sicer", wildcards.name, wildcards.name + "_sicer.broadPeak"), 'w')
            h.write( "\n".join(outBroadPeak) )
            h.close()

if se =="yes":
    rule GEM:
        input:
            chip = join(workpath,bam_dir,"{name}.sorted.Q5DD.tagAlign.gz"),
        output:
            join(workpath,"gem","{name}","{name}.GEM_events.narrowPeak"),         
        params:
            rname='pl:GEM',
            gemver=config['bin'][pfamily]['tool_versions']['GEMVER'],
            readDist=config['bin'][pfamily]['GEMREADDIST'],
            genome = config['references']['ChIPseq']['REFLEN'],
            fastas = config['references']['ChIPseq']['GENOMECHR'],
            gem_dir = "gem",
            ctrl = lambda w : join(workpath,bam_dir,chip2input[w.name] + ".sorted.Q5DD.tagAlign.gz"),
        threads: 32
        shell: """
if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID; fi
cd /lscratch/$SLURM_JOBID;
module load {params.gemver};
cp {input.chip} chip.tagAlign.gz
gzip -d chip.tagAlign.gz
if [ "{params.ctrl}" != "{workpath}/{bam_dir}/.sorted.Q5DD.tagAlign.gz" ]; then
    cp {params.ctrl} input.tagAlign.gz
    gzip -d input.tagAlign.gz
    java -Xmx30g -jar $GEMJAR --t {threads} --d {params.readDist} --g {params.genome} \
         --genome {params.fastas}  --expt chip.tagAlign --ctrl input.tagAlign \
         --out {workpath}/{params.gem_dir}/{wildcards.name} --k_min 6 --k_max 13 --outNP --nrf
else
    java -Xmx30g -jar $GEMJAR --t {threads} --d {params.readDist} --g {params.genome} \
         --genome {params.fastas}  --expt chip.tagAlign \
         --out {workpath}/{params.gem_dir}/{wildcards.name} --k_min 6 --k_max 13 --outNP --nrf
fi
"""

if pe == "yes":
    rule GEM:
        input:
            chip = join(workpath,bam_dir,"{name}.sorted.Q5DD.bam"),
        output:
            join(workpath,"gem","{name}","{name}.GEM_events.narrowPeak"),         
        params:
            rname='pl:GEM',
            gemver=config['bin'][pfamily]['tool_versions']['GEMVER'],
            readDist=config['bin'][pfamily]['GEMREADDIST'],
            genome = config['references']['ChIPseq']['REFLEN'],
            fastas = config['references']['ChIPseq']['GENOMECHR'],
            gem_dir = "gem",
            ctrl = lambda w : join(workpath,bam_dir,chip2input[w.name] + ".sorted.Q5DD.bam"),
        threads: 32
        shell: """
if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID; fi
cd /lscratch/$SLURM_JOBID;
module load {params.gemver};
if [ {params.ctrl} != "{workpath}/{bam_dir}/.sorted.Q5DD.bam" ]; then
    java -Xmx30g -jar $GEMJAR --t {threads} --d {params.readDist} --g {params.genome} \
         --genome {params.fastas}  --expt {input.chip} --ctrl {params.ctrl} --f SAM \
         --out {workpath}/{params.gem_dir}/{wildcards.name} --k_min 6 --k_max 13 --outNP --nrf
else
    java -Xmx30g -jar $GEMJAR --t {threads} --d {params.readDist} --g {params.genome} \
         --genome {params.fastas}  --expt {input.chip} --f SAM \
         --out {workpath}/{params.gem_dir}/{wildcards.name} --k_min 6 --k_max 13 --outNP --nrf
fi
"""

if False:
  rule ChIPQC:
    input:
        lambda w: [ join(workpath, w.PeakTool, chip, chip + PeakExtensions[w.PeakTool]) for chip in chips ]
    output:
        join(workpath,chipQC_dir,'{PeakTool}','ChIPQCreport.html'),
    params:
        rname="pl:ChIPQC",
        genomever = config['project']['annotation'],
        Rver = config['bin'][pfamily]['tool_versions']['RVER'],
        chrom = "chr1"
    run:
        samplesheet = ["\t".join(["SampleID","Condition", "Replicate", "bamReads", "ControlID", "bamControl", "Peaks", "PeakCaller"])]
        for chip in chips:
            file = join(workpath, wildcards.PeakTool, chip, chip + PeakExtensions[wildcards.PeakTool])
            peaksNum = peaks_per_chrom(file, params.chrom)
            if peaksNum > 1:
                condition = [ key for key,value in groupdata.items() if chip in value ][0]
                replicate = str([ i + 1 for i in range(len(groupdata[condition])) if groupdata[condition][i]== chip ][0])
                bamReads = join(workpath, bam_dir, chip + ".sorted.Q5DD.bam")
                controlID = chip2input[chip]
                if controlID != "":
                    bamControl = join(workpath, bam_dir, controlID + ".sorted.Q5DD.bam")
                else:
                    bamControl = ""
                peaks = join(workpath, wildcards.PeakTool, chip, chip + PeakExtensions[wildcards.PeakTool])
                peakcaller = FileTypesChIPQC[wildcards.PeakTool]
                samplesheet.append("\t".join([chip, condition, replicate, bamReads, controlID, bamControl, peaks, peakcaller]))

        csvfile = join(workpath,chipQC_dir,wildcards.PeakTool + "_ChIPQC_prep.csv") 
        f = open(csvfile, 'w')
        f.write ("\n".join(samplesheet))
        f.close()

        R_string = """
        library(ChIPQC)
        samples <- read.table('{tab}', header=1, sep="\\t")
        samples
        result <- ChIPQC(samples, annotation='{genome}', chromosomes='{chr}')
        result
        ChIPQCreport( result, reportName="ChIPQCreport", reportFolder="ChIPQC/{caller}" )
        """.format( tab=csvfile, caller=wildcards.PeakTool, genome=params.genomever , chr=params.chrom )
        
        rscript_fn = join(workpath, chipQC_dir, "chipqc_run_" + wildcards.PeakTool + ".R")
        of=open(rscript_fn,'w')
        of.write(R_string)
        of.close()
        
        shell("module load {params.Rver}; Rscript {rscript_fn}")

# double list of sample1 and sample2 with different names are designed as a work around to deal with snakemake flow rules

rule IDR:
    input:
        sample1 = lambda w: [join(workpath, w.PeakTool, w.sample1, w.sample1 + PeakExtensions[w.PeakTool])],
        sample2 = lambda w: [join(workpath, w.PeakTool, w.sample2, w.sample2 + PeakExtensions[w.PeakTool])],
    output:
        join(workpath,idr_dir,'{PeakTool}','{group}','{sample1}_vs_{sample2}.idrValue.txt')
    params:
        rname = "pl:IDR",
        idrver = config['bin'][pfamily]['tool_versions']['IDRVER'],
        intype = lambda w: FileTypesIDR[w.PeakTool],
        fldr = join(workpath,idr_dir, '{PeakTool}'),
        rank = lambda w: RankColIDR[w.PeakTool],
        sample1 = lambda w: [join(workpath, w.PeakTool, w.sample1, w.sample1 + PeakExtensionsIDR[w.PeakTool])],
        sample2 = lambda w: [join(workpath, w.PeakTool, w.sample2, w.sample2 + PeakExtensionsIDR[w.PeakTool])],
    shell: """
module load {params.idrver}
test -d  {params.fldr} || mkdir {params.fldr}
idr -s {params.sample1} {params.sample2} -o {output} --input-file-type {params.intype} --plot --rank {params.rank}
"""

rule jaccard:
    input:
        lambda w: [ join(workpath, w.PeakTool, chip, chip + PeakExtensions[w.PeakTool]) for chip in chips ],
    output:
        join(workpath,qc_dir,'{PeakTool}_jaccard.txt'),
    params:
        rname="pl:jaccard",
        script=join(workpath,"Scripts","jaccard_score.py"),
        genome = config['references']['ChIPseq']['REFLEN']
    shell: """
python {params.script} -i "{input}" -o {output} -g {params.genome}
"""

rule HOMER_motif:
    input:
        lambda w: [ join(workpath, w.PeakTool, w.name, w.name + PeakExtensions[w.PeakTool]) ]
    output:
        join(workpath, homer_dir,'{PeakTool}',"{name}_{PeakTool}")
    params:
        rname="pl:HOMER_motif",
        homerver = config['bin'][pfamily]['tool_versions']['HOMERVER'],
        genomever = config['project']['annotation'],
    threads: 16
    run:
        commoncmd3="module load {params.homerver}; "
        if wildcards.PeakTool in PeakTools_broad:
            cmd="findMotifsGenome.pl {input} {params.genomever} {output} -size given -p {threads} -preparsedDir /lscratch/$SLURM_JOBID"
        else:
            cmd="findMotifsGenome.pl {input} {params.genomever} {output} -p {threads} -preparsedDir /lscratch/$SLURM_JOBID"
        shell(commoncmd3 + cmd)

rule UROPA:
    input:
        lambda w: [ join(workpath, w.PeakTool, w.name, w.name + PeakExtensions[w.PeakTool]) ]
    output:
        join(workpath, uropa_dir, '{PeakTool}', '{name}_{PeakTool}_uropa_allhits.txt')
    params:
        rname="pl:uropa",
        uropaver = config['bin'][pfamily]['tool_versions']['UROPAVER'],
        fldr = join(workpath, uropa_dir, '{PeakTool}'),
        json = join(workpath, uropa_dir, '{PeakTool}','{name}.json'),
        outroot = join(workpath, uropa_dir, '{PeakTool}','{name}_{PeakTool}_uropa'),
        gtf = config['references']['ChIPseq']['GTFFILE'],
        attributefilter = '"filter.attribute":"gene_type","attribute.value":"protein_coding",'
    threads: 4
    shell: """
module load {params.uropaver};
if [ ! -e {params.fldr} ]; then mkdir {params.fldr}; fi
echo '{{"queries":[ ' > {params.json}
echo '      {{ "feature":"gene","distance":5000,{params.attributefilter}"show.attributes":["gene_id", "gene_name","gene_type"] }},' >> {params.json}
echo '      {{ "feature":"gene",{params.attributefilter}"show.attributes":["gene_id", "gene_name","gene_type"] }}],' >> {params.json}
echo '"priority":"Yes",' >> {params.json}
echo '"gtf":"{params.gtf}",' >> {params.json}
echo '"bed": "{input}" }}' >> {params.json}
uropa -i {params.json} -p {params.outroot} -t {threads} -s
"""

rule FRiP:
     input:
        bam = expand(join(workpath,bam_dir,"{name}.sorted.Q5DD.bam"),name=samples),
        bed = expand(join(workpath,'{PeakTool}','{name}','{name}{PeakExt}'),zip,name=zipSample,PeakTool=zipTool,PeakExt=zipExt),
     output:
        join(workpath,"FRiP_barplot.png"),
     params:
        rname="pl:frip",
        pythonver="python/3.5",
        script=join(workpath,"Scripts","frip_plot.py"),
        genome = config['references']['ChIPseq']['REFLEN']
     shell: """
module load {params.pythonver}
python {params.script} -p "{input.bed}" -b "{input.bam}" -g {params.genome}
mv FRiP_table.txt PeakQC/
mv FRiP_scatterplot.png PeakQC/
"""

rule HOMER_annot:
    input:
        lambda w: [ join(workpath, w.PeakTool, w.name, w.name + PeakExtensions[w.PeakTool]) ]
    output:
        join(workpath, homer2_dir,'{PeakTool}',"{name}_{PeakTool}_annotations.txt")
    params:
        rname="pl:HOMER_annot",
        homerver = config['bin'][pfamily]['tool_versions']['HOMERVER'],
        genomever = config['project']['annotation'],
    threads: 16
    shell: """
module load {params.homerver}
annotatePeaks.pl {input} {params.genomever} > {output}
"""

rule diffbind:
    input:
        lambda w: [ join(workpath, w.PeakTool, chip, chip + PeakExtensions[w.PeakTool]) for chip in chips ]
    output:
        join(workpath,diffbind_dir,"{PeakTool}","{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_DiffBind.html")
    params:
        rname="pl:diffbind",
        Rver = config['bin'][pfamily]['tool_versions']['RVER'],
        rscript1 = join(workpath,"Scripts","runDiffBind.R"),
        rscript2 = join(workpath,"Scripts","DiffBind_pipeliner.Rmd"),
        projectId = config['project']['id'],
        projDesc=config['project']['description'].rstrip('\n'),
        outdir = join(workpath,diffbind_dir,"{PeakTool}","{group1}_vs_{group2}-{PeakTool}"),
        contrast = "{group1}_vs_{group2}"
    run:
        samplesheet = [",".join(["SampleID","Condition", "Replicate", "bamReads", "ControlID", "bamControl", "Peaks", "PeakCaller"])]
        for chip in chips:
            file = join(workpath, wildcards.PeakTool, chip, chip + PeakExtensions[wildcards.PeakTool])
            condition = [ key for key,value in groupdata.items() if chip in value ][0]
            replicate = str([ i + 1 for i in range(len(groupdata[condition])) if groupdata[condition][i]== chip ][0])
            bamReads = join(workpath, bam_dir, chip + ".sorted.Q5DD.bam")
            controlID = chip2input[chip]
            if controlID != "":
                bamControl = join(workpath, bam_dir, controlID + ".sorted.Q5DD.bam")
            else:
                bamControl = ""
            peaks = join(workpath, wildcards.PeakTool, chip, chip + PeakExtensions[wildcards.PeakTool])
            peakcaller = FileTypesDiffBind[wildcards.PeakTool]
            samplesheet.append(",".join([chip, condition, replicate, bamReads, controlID, bamControl, peaks, peakcaller]))

        csvfile = join(workpath,diffbind_dir,wildcards.PeakTool,wildcards.PeakTool + "_diffbind_prep.csv") 
        f = open(csvfile, 'w')
        f.write ("\n".join(samplesheet))
        f.close()
        cmd1 = "module load {params.Rver}; cp {params.rscript1} {params.outdir}; cp {params.rscript2} {params.outdir}; "
        cmd2 = "Rscript runDiffBind.R '.' '{output.html}' '{params.contrast}' '{wildcards.PeakTool}' '{params.projectID}' '{params.projDesc}'"
        shell( cmd1 + cmd2 )

if pe == "yes":
    rule manorm:
        input: 
            bam1 = lambda w: join(workpath,bam_dir, groupdata[w.group1][0] + ".sorted.Q5DD.bam"),
            ppqt1 = lambda w: join(workpath,bam_dir, groupdata[w.group1][0] + ".sorted.Q5DD.ppqt"),
            bam2 = lambda w: join(workpath,bam_dir, groupdata[w.group2][0] + ".sorted.Q5DD.bam"),
            ppqt2 = lambda w: join(workpath,bam_dir, groupdata[w.group2][0] + ".sorted.Q5DD.ppqt"),
            peak1 = lambda w: join(workpath, w.tool, groupdata[w.group1][0], groupdata[w.group1][0] + PeakExtensions[w.tool]),
            peak2 = lambda w: join(workpath, w.tool, groupdata[w.group2][0], groupdata[w.group2][0] + PeakExtensions[w.tool]),
        output:
            fldr = join(workpath,manorm_dir,"{group1}_vs_{group2}-{tool}"),
            xls = join(workpath,manorm_dir,"{group1}_vs_{group2}-{tool}","{group1}_vs_{group2}-{tool}_all_MAvalues.xls"),
            bed = tmp(join(workpath,manorm_dir,"{group1}_vs_{group2}-{tool}","{group1}_vs_{group2}-{tool}_all_MA.bed"))
        params:
            rname='pl:manorm',
            bedtoolsver=config['bin'][pfamily]['tool_versions']['BEDTOOLSVER'],
            manormver="manorm/1.1.4"
        run:
            commoncmd1 = "if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID; fi "
            commoncmd2 = "cd /lscratch/$SLURM_JOBID; "
            commoncmd3 = "module load {params.manormver}; module load {params.bedtoolsver}; "
            cmd1 = "bamToBed -i {input.bam1} > bam1.bed; "
            cmd2 = "bamToBed -i {input.bam2} > bam2.bed; "
            cmd3 = "cut -f 1,2,3 {input.peak1} > peak1.bed; "
            cmd4 = "cut -f 1,2,3 {input.peak2} > peak2.bed; "
            file=list(map(lambda z:z.strip().split(),open(input.ppqt1,'r').readlines()))
            extenders = []
            for ppqt_value in file[0][2].split(","):
                if int(ppqt_value) > 1:
                    extenders.append(ppqt_value)
            try:
                extsize1 = extenders[0]
            except IndexError:
                extsize1 = "{} {}".format(file[0][2].split(",")[0], "# Negative Value which will cause pipeline to fail (wrong ref genome selected or low starting DNA)")               
            file=list(map(lambda z:z.strip().split(),open(input.ppqt2,'r').readlines()))
            extenders = []
            for ppqt_value in file[0][2].split(","):
                if int(ppqt_value) > 1:
                    extenders.append(ppqt_value)
            try:
                extsize2 = extenders[0]
            except IndexError:
                extsize2 = "{} {}".format(file[0][2].split(",")[0], "# Negative Value which will cause pipeline to fail (wrong ref genome selected or low starting DNA)")               
            cmd5 = "manorm --p1 peak1.bed --p2 peak2.bed --r1 bam1.bed --r2 bam2.bed --s1 " + extsize1  + " --s2 " + extsize2 + " -o {output.fldr} --name1 '" + wildcards.group1 + "' --name2 '" + wildcards.group2 + "'; "
            cmd6 = "gzip {output.fldr}/output_tracks/*wig; "
            cmd7 = "mv {output.fldr}/" + wildcards.group1 + "_vs_" + wildcards.group2 + "_all_MAvalues.xls {output.xls}; "
            cmd8 = "tail -n +2 {output.xls} | nl -w2 | awk -v OFS='\t' '{{print $2,$3,$4,$9$1,$6}}' > {output.bed}"
            shell(commoncmd1)
            shell( commoncmd2 + commoncmd3 + cmd1 + cmd2 + cmd3 + cmd4 + cmd5 + cmd6 + cmd7 + cmd8 )