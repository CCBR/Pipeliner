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

for group, chips in groupdata.items() :
    tmp = [ ]
    if len(chips) > 1:
        groupswreps.append(group)
    for chip in chips :
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


#########
# PREPARING TO DEAL WITH A VARIED SET OF PEAKCALL TOOLS

PeakTools_narrow = [ "macs_narrow", "gem" ]
PeakTools_broad = [ "macs_broad", "sicer" ]

PeakTools = PeakTools_narrow + PeakTools_broad

PeakExtensions = { 'macs_narrow': '_peaks.narrowPeak', 'macs_broad': '_peaks.broadPeak',
                   'sicer': '_broadpeaks.bed', 'gem': '.GEM_events.narrowPeak' }

FileTypesChIPQC = { 'macs_narrow': 'narrowPeak', 'macs_broad': 'narrowPeak',
              'sicer': 'bed', 'gem': 'narrowPeak' }

PeakExtensionsIDR = { 'macs_narrow': '_peaks.narrowPeak', 'macs_broad': '_peaks.broadPeak',
                   'sicer': '_sicer.broadPeak', 'gem': '.GEM_events.narrowPeak' }

FileTypesIDR = { 'macs_narrow': 'narrowPeak', 'macs_broad': 'broadPeak',
              'sicer': 'broadPeak', 'gem': 'narrowPeak' }


#########
# CREATING DIRECTORIES

bam_dir='bam'
qc_dir='QC'

chipQC_dir = "ChIPQC"
idr_dir = 'IDR'
memechip_dir = "MEME"
homer_dir = "HOMER_motifs"
uropa_dir = "UROPA_annotations"

otherDirs = [chipQC_dir, homer_dir, uropa_dir]
if reps == "yes":
    otherDirs.append(idr_dir)

for d in PeakTools + otherDirs:
        if not os.path.exists(join(workpath,d)):
                os.mkdir(join(workpath,d))

##########
# RULE ALL

if reps == "yes":
    rule ChIPseq:
        params:
            batch='--time=168:00:00'
        input:
            expand(join(workpath,"macs_narrow","{name}","{name}_peaks.narrowPeak"),name=chips),
            expand(join(workpath,"macs_broad","{name}","{name}_peaks.broadPeak"),name=chips),
            expand(join(workpath,"sicer","{name}","{name}_sicer.broadPeak"),name=chips),
            expand(join(workpath,"gem","{name}","{name}.GEM_events.narrowPeak"),name=chips),
            expand(join(workpath,chipQC_dir,"{PeakTool}","ChIPQCreport.html"),PeakTool=PeakTools),
            expand(join(workpath,qc_dir,'{PeakTool}_jaccard.txt'),PeakTool=PeakTools),
            expand(join(workpath,homer_dir,'{PeakTool}',"{name}_{PeakTool}"),PeakTool=PeakTools,name=chips),
#            expand(join(workpath, uropa_dir,'{PeakTool}','{name}_{PeakTool}_uropa_allhits.txt'),PeakTool=PeakTools,name=chips),
	    expand(join(workpath,idr_dir,'{PeakTool}','{group}.idrValue.txt'),PeakTool=PeakTools, group=groupswreps),
else:
    rule ChIPseq:
        params:
            batch='--time=168:00:00'
        input:
            expand(join(workpath,"macs_narrow","{name}","{name}_peaks.narrowPeak"),name=chips),
            expand(join(workpath,"macs_broad","{name}","{name}_peaks.broadPeak"),name=chips),
            expand(join(workpath,"sicer","{name}","{name}_sicer.broadPeak"),name=chips),
            expand(join(workpath,"gem","{name}","{name}.GEM_events.narrowPeak"),name=chips),
            expand(join(workpath,chipQC_dir,"{PeakTool}","ChIPQCreport.html"),PeakTool=PeakTools),
            expand(join(workpath,qc_dir,'{PeakTool}_jaccard.txt'),PeakTool=PeakTools),

##########
# INDIVIDUAL RULES

rule MACS2_narrow:
    input:
         chip = join(workpath,bam_dir,"{name}.sorted.Q5MDD.tagAlign.gz") if \
            se == "yes" else \
            join(workpath,bam_dir,"{name}.sorted.Q5DD.bam"),
         ctrl = lambda w : join(workpath,bam_dir,chip2input[w.name] + ".sorted.Q5MDD.tagAlign.gz") if \
            se == "yes" else \
            lambda w : join(workpath,bam_dir,chip2input[w.name] + ".sorted.Q5DD.bam"),
         ppqt = join(workpath,bam_dir,"{name}.sorted.Q5MDD.ppqt") if \
            se == "yes" else join(workpath,bam_dir,"{name}.sorted.Q5DD.ppqt"),
    output:
         join(workpath,"macs_narrow","{name}","{name}_peaks.narrowPeak"),
    params:
         rname='pl:MACS2_narrow',
         gsize=config['project']['gsize'],
         macsver=config['bin'][pfamily]['tool_versions']['MACSVER'],
         macsn_dir="macs_narrow"
    shell: """
module load {params.macsver};
ppqt=`cut -f 3 {input.ppqt} | cut -f 1 -d ","`;
# if se == "yes"
if [ {se} == "yes" ]; then 
    if [ "{input.ctrl}" != "{workpath}/{bam_dir}/.sorted.Q5MDD.tagAlign.gz" ]; then
        macs2 callpeak -t {input.chip} -c {input.ctrl} -g {params.gsize} -n {wildcards.name} \
              --outdir {workpath}/{params.macsn_dir}/{wildcards.name} -q 0.01 --keep-dup="all" \
              --extsize $ppqt --nomodel;
    else
        macs2 callpeak -t {input.chip} -g {params.gsize} -n {wildcards.name} \
              --outdir {workpath}/{params.macsn_dir}/{wildcards.name} -q 0.01 --keep-dup="all" \
              --extsize $ppqt --nomodel;
    fi
else
    if [ {input.ctrl} != "{workpath}/{bam_dir}/.sorted.Q5DD.bam" ]; then
        macs2 callpeak -t {input.chip} -c {input.ctrl} -g {params.gsize} -n {wildcards.name} \
              --outdir {workpath}/{params.macsn_dir}/{wildcards.name} -q 0.01 --keep-dup="all" -f "BAMPE";
    else
        macs2 callpeak -t {input.chip} -g {params.gsize} -n {wildcards.name} \
              --outdir {workpath}/{params.macsn_dir}/{wildcards.name} -q 0.01 --keep-dup="all" -f "BAMPE";
    fi
fi
"""

rule MACS2_broad:
    input:
         chip = join(workpath,bam_dir,"{name}.sorted.Q5MDD.tagAlign.gz") if \
            se == "yes" else \
            join(workpath,bam_dir,"{name}.sorted.Q5DD.bam"),
         ctrl = lambda w : join(workpath,bam_dir,chip2input[w.name] + ".sorted.Q5MDD.tagAlign.gz") if \
            se == "yes" else \
            lambda w : join(workpath,bam_dir,chip2input[w.name] + ".sorted.Q5DD.bam"),
         ppqt = join(workpath,bam_dir,"{name}.sorted.Q5MDD.ppqt") if \
            se == "yes" else join(workpath,bam_dir,"{name}.sorted.Q5DD.ppqt"),
    output:
         join(workpath,"macs_broad","{name}","{name}_peaks.broadPeak"),
    params:
         rname='pl:MACS2_broad',
         gsize=config['project']['gsize'],
         macsver=config['bin'][pfamily]['tool_versions']['MACSVER'],
         macsb_dir="macs_broad"
    shell: """
module load {params.macsver};
# if se == "yes"
if [ {se} == "yes" ]; then 
    ppqt=`cut -f 3 {input.ppqt} | cut -f 1 -d ","`;
    if [ "{input.ctrl}" != "{workpath}/{bam_dir}/.sorted.Q5MDD.tagAlign.gz" ]; then
        macs2 callpeak -t {input.chip} -c {input.ctrl} -g {params.gsize} -n {wildcards.name} \
              --outdir {workpath}/{params.macsb_dir}/{wildcards.name} --broad --broad-cutoff 0.01 \
              --keep-dup="all" --extsize $ppqt --nomodel;
    else
        macs2 callpeak -t {input.chip} -g {params.gsize} -n {wildcards.name} \
              --outdir {workpath}/{params.macsb_dir}/{wildcards.name} --broad --broad-cutoff 0.01 \
              --keep-dup="all" --extsize $ppqt --nomodel;
    fi
else
    if [ {input.ctrl} != "{workpath}/{bam_dir}/.sorted.Q5DD.bam" ]; then
        macs2 callpeak -t {input.chip} -c {input.ctrl} -g {params.gsize} -n {wildcards.name} \
              --outdir {workpath}/{params.macsb_dir}/{wildcards.name} --broad --broad-cutoff 0.01 \
              --keep-dup="all" -f "BAMPE";
    else
        macs2 callpeak -t {input.chip} -g {params.gsize} -n {wildcards.name} \
              --outdir {workpath}/{params.macsb_dir}/{wildcards.name} --broad --broad-cutoff 0.01 \
              --keep-dup="all" -f "BAMPE";
    fi
fi
"""

rule SICER:
    input:
         chip = join(workpath,bam_dir,"{name}.sorted.Q5MDD.tagAlign.gz") if \
            se == "yes" else \
            join(workpath,bam_dir,"{name}.sorted.Q5DD.bam"),
         ctrl = lambda w : join(workpath,bam_dir,chip2input[w.name] + ".sorted.Q5MDD.tagAlign.gz") if \
            se == "yes" else \
            lambda w : join(workpath,bam_dir,chip2input[w.name] + ".sorted.Q5DD.bam"),
         ppqt = join(workpath,bam_dir,"{name}.sorted.Q5MDD.ppqt") if \
            se == "yes" else join(workpath,bam_dir,"{name}.sorted.Q5DD.ppqt"),
    output:
         bed = join(workpath,"sicer","{name}","{name}_broadpeaks.bed"),
# output columns: chrom, start, end, ChIP tag count, control tag count, p-value, fold-enrichment, q-value
    params:
         rname='pl:SICER',
         sicerver=config['bin'][pfamily]['tool_versions']['SICERVER'],
         bedtoolsver=config['bin'][pfamily]['tool_versions']['BEDTOOLSVER'],
         genomever = config['project']['annotation'],
    shell: """
if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID; fi
cd /lscratch/$SLURM_JOBID;
module load {params.sicerver};
module load {params.bedtoolsver};
ppqt=`cut -f 3 {input.ppqt} | cut -f 1 -d ","`;
# if se == "yes"
if [ {se} == "yes" ]; then
    cp {input.chip} chip.bed.gz
    gzip -d chip.bed.gz
    if [ "{input.ctrl}" != "{workpath}/{bam_dir}/.sorted.Q5MDD.tagAlign.gz" ]; then
        cp {input.ctrl} input.bed.gz
        gzip -d input.bed.gz
        sh $SICERDIR/SICER.sh . chip.bed input.bed . {params.genome} 100 300 $ppqt 0.75 600 1E-2
	mv chip-W300-G600-islands-summary-FDR1E-2 {output.bed}
    else
        sh $SICERDIR/SICER.sh . chip.tagAlign . {params.genome} 100 300 $ppqt 0.75 600 100
        mv chip-W300-G600-E100.scoreisland {output.bed}
    fi
else
    bamToBed -i {input.chip} > chip.bed
    if [ {input.ctrl} != "{workpath}/{bam_dir}/.sorted.Q5DD.bam" ]; then
        bamToBed -i {input.ctrl} > input.bed
        sh $SICERDIR/SICER.sh . chip.bed input.bed . {params.genomever} 100 300 $ppqt 0.75 600 1E-2
        mv chip-W300-G600-islands-summary-FDR1E-2 {output.bed}
    else
        sh $SICERDIR/SICER.sh . chip.bed . {params.genomever} 100 300 $ppqt 0.75 600 100
        mv chip-W300-G600-E100.scoreisland {output.bed}
    fi
fi
"""

rule convertSICER:
    input:
         bed = join(workpath,"sicer","{name}","{name}_broadpeaks.bed"),
# input columns: chrom, start, end, ChIP tag count, control tag count, p-value, fold-enrichment, q-value
    output:
         broadPeak = join(workpath,"sicer","{name}","{name}_sicer.broadPeak"),
# output columns: chrom, start, end, name, fold-enrichment, strand, ChIP tag count, -log10 p-value, -log10 q-value
    params:
         rname='pl:convertSICER',
    run:
        import math
        f = open(input.bed,'r')
        inbed = f.readlines()
        f.close()
        outPeak = [None] * len(inbed)
        for i in range(len(inbed)):
            tmp = inbed[i].strip().split('\t')
# assuming that a p-value/q-value of 0 is super significant, -log10(1e-500)
            if tmp[5] == "0.0":
                pval="500"
            else:
                pval = str(-(math.log10(float(tmp[5]))))
            if tmp[7] == "0.0":
                qval="500"
            else:
                qval = str(-(math.log10(float(tmp[7]))))
            outPeak[i] = "\t".join(tmp[0:3]+["Peak"+str(i+1),tmp[6],".", tmp[3], pval, qval])
        g = open(output.broadPeak,'w')
        g.write( "\n".join(outPeak) )
        g.close()

rule GEM:
    input:
         chip = join(workpath,bam_dir,"{name}.sorted.Q5MDD.tagAlign.gz") if \
            se == "yes" else \
            join(workpath,bam_dir,"{name}.sorted.Q5DD.bam"),
         ctrl = lambda w : join(workpath,bam_dir,chip2input[w.name] + ".sorted.Q5MDD.tagAlign.gz") if \
            se == "yes" else \
            lambda w : join(workpath,bam_dir,chip2input[w.name] + ".sorted.Q5DD.bam"),
         ppqt = join(workpath,bam_dir,"{name}.sorted.Q5MDD.ppqt") if \
            se == "yes" else join(workpath,bam_dir,"{name}.sorted.Q5DD.ppqt"),
    output:
         join(workpath,"gem","{name}","{name}.GEM_events.narrowPeak"),         
         join(workpath,"gem","{name}","{name}.GEM_events.bed"),         
    params:
         rname='pl:GEM',
         gemver=config['bin'][pfamily]['tool_versions']['GEMVER'],
         readDist=config['bin'][pfamily]['GEMREADDIST'],
         genome = config['references']['ChIPseq']['REFLEN'],
         fastas = config['references']['ChIPseq']['GENOMECHR'],
         gem_dir = "gem"
    threads: 32
    shell: """
if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID; fi
cd /lscratch/$SLURM_JOBID;
module load {params.gemver};
ppqt=`cut -f 3 {input.ppqt} | cut -f 1 -d ","`;
# if se == "yes"
if [ {se} == "yes" ]; then
    cp {input.chip} chip.tagAlign.gz
    gzip -d chip.tagAlign.gz
    if [ "{input.ctrl}" != "{workpath}/{bam_dir}/.sorted.Q5MDD.tagAlign.gz" ]; then
        cp {input.ctrl} input.tagAlign.gz
        gzip -d input.tagAlign.gz
        java -Xmx30g -jar $GEMJAR --t {threads} --d {params.readDist} --g {params.genome} \
             --genome {params.fastas}  --expt chip.tagAlign --ctrl input.tagAlign \
             --out {workpath}/{params.gem_dir}/{wildcards.name} --k_min 6 --k_max 13 --outNP --nrf
    else
        java -Xmx30g -jar $GEMJAR --t {threads} --d {params.readDist} --g {params.genome} \
             --genome {params.fastas}  --expt chip.tagAlign \
             --out {workpath}/{params.gem_dir}/{wildcards.name} --k_min 6 --k_max 13 --outNP --nrf        
    fi
else
    if [ {input.ctrl} != "{workpath}/{bam_dir}/.sorted.Q5DD.bam" ]; then
        java -Xmx30g -jar $GEMJAR --t {threads} --d {params.readDist} --g {params.genome} \
             --genome {params.fastas}  --expt {input.chip} --ctrl {input.ctrl} --f SAM \
             --out {workpath}/{params.gem_dir}/{wildcards.name} --k_min 6 --k_max 13 --outNP --nrf
    else
        java -Xmx30g -jar $GEMJAR --t {threads} --d {params.readDist} --g {params.genome} \
             --genome {params.fastas}  --expt {input.chip} --f SAM \
             --out {workpath}/{params.gem_dir}/{wildcards.name} --k_min 6 --k_max 13 --outNP --nrf
    fi
fi
"""

rule ChIPQC:
    input:
        lambda w: [ join(workpath, w.PeakTool, chip, chip + PeakExtensions[w.PeakTool]) for chip in chips ]
    output:
        join(workpath,chipQC_dir,'{PeakTool}','ChIPQCreport.html'),
    params:
        rname="pl:ChIPQC",
        genomever = config['project']['annotation'],
        Rver = config['bin'][pfamily]['tool_versions']['RVER'],
    run:
        samplesheet = ["\t".join(["SampleID","Condition", "Replicate", "bamReads", "ControlID", "bamControl", "Peaks", "PeakCaller"])]
        for chip in chips:
            condition = [ key for key,value in groupdata.items() if chip in value ][0]
            replicate = str([ i + 1 for i in range(len(groupdata[condition])) if groupdata[condition][i]== chip ][0])
            if se == "yes":
                bamReads = join(workpath, bam_dir, chip + ".sorted.Q5MDD.bam")
            else:
                bamReads = join(workpath, bam_dir, chip + ".sorted.Q5DD.bam")
            controlID = chip2input[chip]
            if controlID != "":
                if se == "yes":
                    bamControl = join(workpath, bam_dir, controlID + ".sorted.Q5MDD.bam")
                else:
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
        """.format( tab=csvfile, caller=wildcards.PeakTool, genome=params.genomever , chr="chr1")
        
        rscript_fn = join(workpath, chipQC_dir, "chipqc_run_" + wildcards.PeakTool + ".R")
        of=open(rscript_fn,'w')
        of.write(R_string)
        of.close()
        
        shell("module load {params.Rver}; Rscript {rscript_fn}")

rule IDR:
    input:
        lambda w: [join(workpath, w.PeakTool, sample, sample + PeakExtensionsIDR[w.PeakTool]) for sample in groupdata[w.group]]
    output:
        join(workpath,idr_dir,'{PeakTool}','{group}.idrValue.txt')
    params:
        rname="pl:IDR",
        idrver=config['bin'][pfamily]['tool_versions']['IDRVER'],
        intype= lambda w: FileTypesIDR[w.PeakTool],
        fldr= join(workpath,idr_dir, '{PeakTool}')
    shell: """
module load {params.idrver}
test -d  {params.fldr} || mkdir {params.fldr}
idr -s {input} -o {output} --input-file-type {params.intype} --plot
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
    threads: 4
    run:
        commoncmd1="if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID; fi "
        commoncmd2="cd /lscratch/$SLURM_JOBID; "
        commoncmd3="module load {params.homerver}; "
        if wildcards.PeakTool in PeakTools_broad:
            cmd="findMotifsGenome.pl {input} {params.genomever} {output} -size given -p {threads} -preparsedDir '.'"
        else:
            cmd="findMotifsGenome.pl {input} {params.genomever} {output} -p {threads} -preparsedDir '.'"
        shell( commoncmd1 + commoncmd2 + commoncmd3 + cmd )

rule UROPA:
    input:
        lambda w: [ join(workpath, w.PeakTool, w.name, w.name + PeakExtensions[w.PeakTool]) ]
    output:
        join(workpath, uropa_dir, '{PeakTool}', '{name}_{PeakTool}_uropa_allhits.txt')
    params:
        rname="pl:uropa",
        Rver = config['bin'][pfamily]['tool_versions']['RVER'],
        fldr = join(workpath, uropa_dir),
        outroot = join(workpath, uropa_dir, '{PeakTool}', '{name}_{PeakTool}_uropa'),
        gtf = config['references']['ChIPseq']['GTFFILE']
    threads: 4
    shell: """
module load {params.Rver};
echo '{{"queries":[
          {{ "feature":"gene","distance":5000,"show.attributes":["gene_id", "gene_name","gene_type"] }},
          {{ "feature":"gene","show.attributes":["gene_id", "gene_name","gene_type"] }}],
       "priority":"Yes",
       "gtf":{params.gtf},
       "bed":{input} }}' > {params.fldr}/{wildcards.name}.json
uropa -i {params.fldr}/{wildcards.name}.json -p  {params.outroot} -s -t {threads}
"""