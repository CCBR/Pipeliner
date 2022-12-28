from snakemake.utils import R
from os.path import join
import re,os

from os import listdir


configfile: "run.json"

workpath = config['project']['workpath']
filetype = config['project']['filetype']
readtype = config['project']['readtype']


def outputfiles2(groupslist, inputnorm):
    '''
    Produces correct output filenames based on group information.
    Names will be:
    Inputnorm.Q5DD.RPGC.metagene_heatmap.pdf
    {groupName}.Q5DD.RPGC.metagene_heatmap.pdf
    {groupName}.sorted.RPGC.metagene_heatmap.pdf
    Note: Inputnorm will only be included when there are input samples.
    '''
    dtoolgroups, dtoolext = [], []
    extensions = ["sorted.RPGC", "Q5DD.RPGC"]
    if len(inputnorm) == 2:
            dtoolgroups.extend(["InputNorm"])
            dtoolext.extend([extensions[1]])
    for group in groupslist:
            dtoolgroups.extend([group] * 2)
            dtoolext.extend([extensions[1], extensions[0]])
    if len(inputnorm) == 2:
            dtoolgroups.extend(["InputNorm.prot"])
            dtoolext.extend([extensions[1]])
    for group in groupslist:
            dtoolgroups.extend([group + ".prot"] * 2)
            dtoolext.extend([extensions[1], extensions[0]])
    return dtoolgroups, dtoolext

se=""
pe=""

if config['project']['nends'] == 2 :
    pe="yes"
elif config['project']['nends'] == 1 :
    se="yes"

if pe == "yes":
    extensions = [ "sorted.RPGC", "Q5DD.RPGC" ]
    extensions2 = list(map(lambda x:re.sub(".RPGC","",x),extensions))
    extensions3 = { extensions2[i] + "." : "bam" for i in range(len(extensions2)) }
    extensions4 = [ extensions2[i] + ".bam" for i in range(len(extensions2)) ]
else:
    extensions = [ "sorted.RPGC", "Q5DD.RPGC" ]
    extensions2 = list(map(lambda x:re.sub(".RPGC","",x),extensions))
    types = [ "bam", "tagAlign.gz" ]
    extensions3 = { extensions2[i] + "." : types[i] for i in range(len(extensions2)) }
    extensions4 = [ extensions2[i] + "." + types[i] for i in range(len(extensions2)) ]

d = {i: True for i in [1,2,3]}

#########
# DEFINING SAMPLES

chip2input = config['project']['peaks']['inputs']
uniq_inputs = list(sorted(set([v for v in chip2input.values() if v])))

sampleswinput = []
for input in chip2input:
	if chip2input[input] != 'NA' and chip2input[input] != '':
		sampleswinput.append(input)

if len(sampleswinput) == 0:
    inputnorm = [""]
else:
    inputnorm = ["",".inputnorm"]

groupdata = config['project']['groups']

groupdatawinput = {}
for group, chips in groupdata.items() :
    tmp = [ ]
    for chip in chips :
        if chip in samples:
            tmp.append(chip)
            input = chip2input[chip]
            if input != 'NA' and input != '':
                tmp.append(input)
    if len(tmp) != 0:
        groupdatawinput[group]=set(tmp)

groups = list(groupdatawinput.keys())
deepgroups, deepexts = outputfiles2(groups,inputnorm)

##########
# CREATING DIRECTORIES

trim_dir='trim'
kraken_dir='kraken'
bam_dir='bam'
bw_dir='bigwig'
deeptools_dir='deeptools'
extra_fingerprint_dir='deeptools/sorted_fingerprint'
qc_dir="QC"

for d in [trim_dir,kraken_dir,bam_dir,bw_dir,deeptools_dir,extra_fingerprint_dir,qc_dir]:
	if not os.path.exists(join(workpath,d)):
		os.mkdir(join(workpath,d))

# Checking to see if output filenames are files are being generated correctly
#print("DEEPTOOLS OUTPUT FILESNAMES\n",[file.split('/')[-1] for file in expand(join(workpath,deeptools_dir,"{group}.metagene_heatmap.{ext}.pdf"), zip, group=deepgroups,ext=deepexts)])
#print(samples)

##########
# RULES

if se == 'yes' :
    rule InitialChIPseqQC:
        params: 
            batch='--time=168:00:00'
        input: 
            # Multiqc Report
            join(workpath,"Reports","multiqc_report.html"),
            # FastqScreen
            expand(join(workpath,"FQscreen","{name}.R1.trim_screen.png"),name=samples),
            expand(join(workpath,"FQscreen2","{name}.R1.trim_screen.png"),name=samples),
            # Kraken
            expand(join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.taxa.txt"),name=samples),
            expand(join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.krona.html"),name=samples),
            join(workpath,kraken_dir,"kraken_bacteria.taxa.summary.txt"),
            # Align using BWA and dedup with Picard or MACS2
            expand(join(workpath,bam_dir,"{name}.{ext}"),name=samples,ext=extensions4),
            # BWA --> BigWig
            expand(join(workpath,bw_dir,"{name}.{ext}.bw",),name=samples,ext=extensions),
            # Input Normalization
            expand(join(workpath,bw_dir,"{name}.Q5DD.RPGC.inputnorm.bw",),name=sampleswinput), 
            # PhantomPeakQualTools
            expand(join(workpath,bam_dir,"{name}.{ext}.ppqt"),name=samples,ext=extensions2),
            expand(join(workpath,bam_dir,'{ext}.ppqt.txt'),ext=extensions2),
            # deeptools
            expand(join(workpath,deeptools_dir,"spearman_heatmap.{ext}.pdf"),ext=extensions),
            expand(join(workpath,deeptools_dir,"spearman_scatterplot.{ext}.pdf"),ext=extensions),
            expand(join(workpath,deeptools_dir,"pca.{ext}.pdf"),ext=extensions),
            expand(join(workpath,deeptools_dir,"{group}.fingerprint.{ext}.pdf"),ext=extensions2,group=groups),
            expand(join(workpath,deeptools_dir,"{group}.metagene_heatmap.{ext}.pdf"), zip, group=deepgroups,ext=deepexts),
            expand(join(workpath,deeptools_dir,"{group}.TSS_heatmap.{ext}.pdf"), zip, group=deepgroups,ext=deepexts),
            expand(join(workpath,deeptools_dir,"{group}.metagene_profile.{ext}.pdf"), zip, group=deepgroups,ext=deepexts),
            expand(join(workpath,deeptools_dir,"{group}.TSS_profile.{ext}.pdf"), zip, group=deepgroups,ext=deepexts),
            # QC Table
            join(workpath,qc_dir,"QCTable.txt"),


    rule trim_se: # actually trim, filter polyX and remove black listed reads
        input:
            infq=join(workpath,"{name}.R1.fastq.gz"),
        output:
            outfq=temp(join(workpath,trim_dir,"{name}.R1.trim.fastq.gz")),
        params:
            rname="pl:trim",
            cutadaptver=config['bin'][pfamily]['tool_versions']['CUTADAPTVER'],
            workpath=config['project']['workpath'],
            adaptersfa=config['bin'][pfamily]['FASTAWITHADAPTERSETD'],
            blacklistbwaindex=config['references'][pfamily]['BLACKLISTBWAINDEX'],
            picardver=config['bin'][pfamily]['tool_versions']['PICARDVER'],
            bwaver=config['bin'][pfamily]['tool_versions']['BWAVER'],
            samtoolsver=config['bin'][pfamily]['tool_versions']['SAMTOOLSVER'],
            minlen=35,
            leadingquality=10,
            trailingquality=10,
            javaram="100g",
        threads: 16
        shell: """
module load {params.cutadaptver};
if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID ;fi
cd /lscratch/$SLURM_JOBID
sample=`echo {input.infq}|awk -F "/" '{{print $NF}}'|awk -F ".R1.fastq" '{{print $1}}'`
cutadapt --nextseq-trim=2 --trim-n -n 5 -O 5 -q {params.leadingquality},{params.trailingquality} -m {params.minlen} -b file:{params.adaptersfa} -j {threads} -o ${{sample}}.cutadapt.fastq {input.infq}
module load {params.bwaver};
module load {params.samtoolsver};
module load {params.picardver};
bwa mem -t {threads} {params.blacklistbwaindex} ${{sample}}.cutadapt.fastq | samtools view -@{threads} -f4 -b -o ${{sample}}.bam
java -Xmx{params.javaram} -jar $PICARDJARPATH/picard.jar SamToFastq VALIDATION_STRINGENCY=SILENT INPUT=${{sample}}.bam FASTQ=${{sample}}.cutadapt.noBL.fastq
pigz -p 16 ${{sample}}.cutadapt.noBL.fastq;
mv ${{sample}}.cutadapt.noBL.fastq.gz {output.outfq};
            """

    rule kraken_se:
        input:
            fq = join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        output:
            krakentaxa = join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.taxa.txt"),
            kronahtml = join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.krona.html"),
        params: 
            rname='pl:kraken',
            prefix="{name}",
            outdir=join(workpath,kraken_dir),
            bacdb=config['bin'][pfamily]['KRAKENBACDB'],
            krakenver=config['bin'][pfamily]['tool_versions']['KRAKENVER'],
            kronatoolsver=config['bin'][pfamily]['tool_versions']['KRONATOOLSVER'],
        threads: 16
        shell: """
module load {params.krakenver};
module load {params.kronatoolsver};
cd /lscratch/$SLURM_JOBID;
cp -rv {params.bacdb} /lscratch/$SLURM_JOBID/;
kraken --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` --fastq-input --gzip-compressed --threads {threads} --output /lscratch/$SLURM_JOBID/{params.prefix}.krakenout --preload {input.fq}
kraken-translate --mpa-format --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` /lscratch/$SLURM_JOBID/{params.prefix}.krakenout |cut -f2|sort|uniq -c|sort -k1,1nr > /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa
cut -f2,3 {params.prefix}.krakenout | ktImportTaxonomy - -o {params.prefix}.kronahtml
mv /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa {output.krakentaxa}
mv /lscratch/$SLURM_JOBID/{params.prefix}.kronahtml {output.kronahtml}
            """

    rule process_kraken:
        input:
            fq = expand(join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),name=samples),
            krakentaxa = expand(join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.taxa.txt"),name=samples),
        output:
            kraken_taxa_summary = join(workpath,kraken_dir,"kraken_bacteria.taxa.summary.txt"),
        params:
            rname = "pl:krakenProcess",
        run:
            cmd="echo -ne \"Sample\tPercent\tBacteria\n\" > "+output.kraken_taxa_summary
            for f,t in zip(input.fq,input.krakentaxa):
                cmd="sh Scripts/kraken_process_taxa.sh "+f+" "+t+" >> "+output.kraken_taxa_summary
                shell(cmd)

    rule BWA_SE:
        input:
            infq=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        params:
            d=join(workpath,bam_dir),
            rname='pl:bwa',
            reference=config['references'][pfamily]['BWA'],
            bwaver=config['bin'][pfamily]['tool_versions']['BWAVER'],
            samtoolsver=config['bin'][pfamily]['tool_versions']['SAMTOOLSVER'],
        output:
            outbam1=join(workpath,bam_dir,"{name}.sorted.bam"), 
            outbam2=temp(join(workpath,bam_dir,"{name}.Q5.bam")),
            flagstat1=join(workpath,bam_dir,"{name}.sorted.bam.flagstat"),
            flagstat2=join(workpath,bam_dir,"{name}.Q5.bam.flagstat"),
            idxstat1=join(workpath,bam_dir,"{name}.sorted.bam.idxstat"),
            idxstat2=join(workpath,bam_dir,"{name}.Q5.bam.idxstat"),
        threads: 32
        shell: """
module load {params.bwaver};
module load {params.samtoolsver};
bwa mem -t {threads} {params.reference} {input.infq} | \
samtools sort -@{threads} -o {output.outbam1}
samtools index {output.outbam1}
samtools flagstat {output.outbam1} > {output.flagstat1}
samtools idxstats {output.outbam1} > {output.idxstat1}
samtools view -b -q 6 {output.outbam1} -o {output.outbam2}
samtools index {output.outbam2}
samtools flagstat {output.outbam2} > {output.flagstat2}
samtools idxstats {output.outbam2} > {output.idxstat2}
            """  

    rule macs2_dedup:
        input:
            join(workpath,bam_dir,"{name}.Q5.bam")
        output:
            outtagalign=join(workpath,bam_dir,"{name}.Q5DD.tagAlign.gz"),
            outbam=join(workpath,bam_dir,"{name}.Q5DD.bam"),
            outflagstat=join(workpath,bam_dir,"{name}.Q5DD.bam.flagstat"),
            outidxstat=join(workpath,bam_dir,"{name}.Q5DD.bam.idxstat"),
        params:
            rname='pl:macs2_dedup',
            macsver=config['bin'][pfamily]['tool_versions']['MACSVER'],
            samtoolsver=config['bin'][pfamily]['tool_versions']['SAMTOOLSVER'],
            bedtoolsver=config['bin'][pfamily]['tool_versions']['BEDTOOLSVER'],
            gsize=config['references'][pfamily]['EFFECTIVEGENOMESIZE'],
            folder=join(workpath,bam_dir),
	    genomefile=config['references'][pfamily]['REFLEN']
        shell: """
if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID; fi
cd /lscratch/$SLURM_JOBID;
module load {params.macsver};
module load {params.samtoolsver};
module load {params.bedtoolsver};
macs2 filterdup -i {input} -g {params.gsize} --keep-dup="auto" -o TmpTagAlign;
awk -F"\\t" -v OFS="\\t" '{{if ($2>0 && $3>0) {{print}}}}' TmpTagAlign > TmpTagAlign2;
awk -F"\\t" -v OFS="\\t" '{{print $1,1,$2}}' {params.genomefile} | sort -k1,1 -k2,2n > GenomeFileBed;
bedtools intersect -wa -f 1.0 -a TmpTagAlign2 -b GenomeFileBed > TmpTagAlign3;
bedtools bedtobam -i TmpTagAlign3 -g {params.genomefile} | samtools sort -@4 -o {output.outbam};
gzip TmpTagAlign3;
mv TmpTagAlign3.gz {output.outtagalign};
samtools index {output.outbam};
samtools flagstat {output.outbam} > {output.outflagstat}
samtools idxstats {output.outbam} > {output.outidxstat}
"""

if pe == 'yes':
    rule InitialChIPseqQC:
        params: 
            batch='--time=168:00:00'
        input: 
            # Multiqc Report
            join(workpath,"Reports","multiqc_report.html"),
            # FastqScreen
            expand(join(workpath,"FQscreen","{name}.R{rn}.trim_screen.png"),name=samples,rn=[1,2]),
            expand(join(workpath,"FQscreen2","{name}.R{rn}.trim_screen.png"),name=samples,rn=[1,2]),
            # Kraken
            expand(join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.taxa.txt"),name=samples),
            expand(join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.krona.html"),name=samples),
            # Align using BWA and dedup with Picard
            expand(join(workpath,bam_dir,"{name}.{ext}"),name=samples,ext=extensions4),
            # BWA --> BigWig
            expand(join(workpath,bw_dir,"{name}.{ext}.bw",),name=samples,ext=extensions),
            # Input Normalization
            expand(join(workpath,bw_dir,"{name}.Q5DD.RPGC.inputnorm.bw",),name=sampleswinput),
            # PhantomPeakQualTools
            expand(join(workpath,bam_dir,"{name}.{ext}.ppqt"),name=samples,ext=extensions2),
            expand(join(workpath,bam_dir,"{name}.{ext}.pdf"),name=samples,ext=extensions2),
            # deeptools
            expand(join(workpath,deeptools_dir,"spearman_heatmap.{ext}.pdf"),ext=extensions),
            expand(join(workpath,deeptools_dir,"spearman_scatterplot.{ext}.pdf"),ext=extensions),
            expand(join(workpath,deeptools_dir,"pca.{ext}.pdf"),ext=extensions),
       	    expand(join(workpath,deeptools_dir,"{group}.fingerprint.{ext}.pdf"),group=groups,ext=extensions2),
            expand(join(workpath,deeptools_dir,"{group}.metagene_heatmap.{ext}.pdf"), zip, group=deepgroups,ext=deepexts),
            expand(join(workpath,deeptools_dir,"{group}.TSS_heatmap.{ext}.pdf"), zip, group=deepgroups,ext=deepexts),
            expand(join(workpath,deeptools_dir,"{group}.metagene_profile.{ext}.pdf"), zip, group=deepgroups,ext=deepexts),
            expand(join(workpath,deeptools_dir,"{group}.TSS_profile.{ext}.pdf"), zip, group=deepgroups,ext=deepexts),
            # QC Table
            join(workpath,qc_dir,"QCTable.txt"),
            expand(join(workpath,qc_dir,"{name}.{ext}.insert_size_metrics.txt"),name=samples,ext=extensions2),


    rule trim_pe: # trim, remove PolyX and remove BL reads
        input:
            file1=join(workpath,"{name}.R1.fastq.gz"),
            file2=join(workpath,"{name}.R2.fastq.gz"),
        output:
            outfq1=temp(join(workpath,trim_dir,"{name}.R1.trim.fastq.gz")),
            outfq2=temp(join(workpath,trim_dir,"{name}.R2.trim.fastq.gz")),
        params:
            rname="pl:trim",
            cutadaptver=config['bin'][pfamily]['tool_versions']['CUTADAPTVER'],
            workpath=config['project']['workpath'],
            fastawithadaptersetd=config['bin'][pfamily]['FASTAWITHADAPTERSETD'],
            blacklistbwaindex=config['references'][pfamily]['BLACKLISTBWAINDEX'],
            picardver=config['bin'][pfamily]['tool_versions']['PICARDVER'],
            bwaver=config['bin'][pfamily]['tool_versions']['BWAVER'],
            samtoolsver=config['bin'][pfamily]['tool_versions']['SAMTOOLSVER'],
            minlen=35,
            leadingquality=10,
            trailingquality=10,
            javaram="100g",
        threads: 16
        shell: """
module load {params.cutadaptver};
if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID ;fi
cd /lscratch/$SLURM_JOBID
sample=`echo {input.file1}|awk -F "/" '{{print $NF}}'|awk -F ".R1.fastq" '{{print $1}}'`
cutadapt --pair-filter=any --nextseq-trim=2 --trim-n -n 5 -O 5 -q {params.leadingquality},{params.trailingquality} -m {params.minlen}:{params.minlen} -b file:{params.fastawithadaptersetd} -B file:{params.fastawithadaptersetd} -j {threads} -o ${{sample}}.R1.cutadapt.fastq -p ${{sample}}.R2.cutadapt.fastq {input.file1} {input.file2}
module load {params.bwaver};
module load {params.samtoolsver};
module load {params.picardver};
bwa mem -t {threads} {params.blacklistbwaindex} ${{sample}}.R1.cutadapt.fastq ${{sample}}.R2.cutadapt.fastq | samtools view -@{threads} -f4 -b -o ${{sample}}.bam
java -Xmx{params.javaram} -jar $PICARDJARPATH/picard.jar SamToFastq \
VALIDATION_STRINGENCY=SILENT \
INPUT=${{sample}}.bam \
FASTQ=${{sample}}.R1.cutadapt.noBL.fastq \
SECOND_END_FASTQ=${{sample}}.R2.cutadapt.noBL.fastq \
UNPAIRED_FASTQ=${{sample}}.unpaired.noBL.fastq
pigz -p {threads} ${{sample}}.R1.cutadapt.noBL.fastq;
pigz -p {threads} ${{sample}}.R2.cutadapt.noBL.fastq;
mv /lscratch/$SLURM_JOBID/${{sample}}.R1.cutadapt.noBL.fastq.gz {output.outfq1};
mv /lscratch/$SLURM_JOBID/${{sample}}.R2.cutadapt.noBL.fastq.gz {output.outfq2};
"""

    rule kraken_pe:
        input:
            fq1 = join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
            fq2 = join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"),
        output:
            krakentaxa = join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.taxa.txt"),
            kronahtml = join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.krona.html"),
        params: 
            rname='pl:kraken',
            prefix="{name}",
            outdir=join(workpath,kraken_dir),
            bacdb=config['bin'][pfamily]['KRAKENBACDB'],
            krakenver=config['bin'][pfamily]['tool_versions']['KRAKENVER'],
            kronatoolsver=config['bin'][pfamily]['tool_versions']['KRONATOOLSVER'],
        threads: 32
        shell: """
module load {params.krakenver};
module load {params.kronatoolsver};
if [ ! -d {params.outdir} ];then mkdir {params.outdir};fi
cd /lscratch/$SLURM_JOBID;
cp -rv {params.bacdb} /lscratch/$SLURM_JOBID/;
kraken --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` --fastq-input --gzip-compressed --threads {threads} --output /lscratch/$SLURM_JOBID/{params.prefix}.krakenout --preload --paired {input.fq1} {input.fq2}
kraken-translate --mpa-format --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` /lscratch/$SLURM_JOBID/{params.prefix}.krakenout |cut -f2|sort|uniq -c|sort -k1,1nr > /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa
cut -f2,3 /lscratch/$SLURM_JOBID/{params.prefix}.krakenout | ktImportTaxonomy - -o /lscratch/$SLURM_JOBID/{params.prefix}.kronahtml
mv /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa {output.krakentaxa}
mv /lscratch/$SLURM_JOBID/{params.prefix}.kronahtml {output.kronahtml}
"""

    rule BWA_PE:
        input:
            infq1 = join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
            infq2 = join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"),
        params:
            d=join(workpath,bam_dir),
            rname='pl:bwa',
            reference=config['references'][pfamily]['BWA'],
            bwaver=config['bin'][pfamily]['tool_versions']['BWAVER'],
            samtoolsver=config['bin'][pfamily]['tool_versions']['SAMTOOLSVER'],
            script=join(workpath,"Scripts","bam_filter_by_mapq.py"),
            pythonver="python/3.5"
        output:
            outbam1=join(workpath,bam_dir,"{name}.sorted.bam"), 
            outbam2=temp(join(workpath,bam_dir,"{name}.Q5.bam")),
            flagstat1=join(workpath,bam_dir,"{name}.sorted.bam.flagstat"),
            idxstat1=join(workpath,bam_dir,"{name}.sorted.bam.idxstat"),
            flagstat2=join(workpath,bam_dir,"{name}.Q5.bam.flagstat"),
            idxstat2=join(workpath,bam_dir,"{name}.Q5.bam.idxstat"),
        threads: 32
        shell: """
module load {params.bwaver};
module load {params.samtoolsver};
module load {params.pythonver};
bwa mem -t {threads} {params.reference} {input.infq1} {input.infq2} | \
samtools sort -@{threads} -o {output.outbam1}
samtools index {output.outbam1}
samtools flagstat {output.outbam1} > {output.flagstat1}
samtools idxstats {output.outbam1} > {output.idxstat1}
#samtools view -b -q 6 {output.outbam1} -o {output.outbam2}
python {params.script} -i {output.outbam1} -o {output.outbam2} -q 6
samtools index {output.outbam2}
samtools flagstat {output.outbam2} > {output.flagstat2}
samtools idxstats {output.outbam2} > {output.idxstat2}
            """  

    rule picard_dedup:
        input: 
            bam2=join(workpath,bam_dir,"{name}.Q5.bam")
        output:
            out4=temp(join(workpath,bam_dir,"{name}.bwa_rg_added.Q5.bam")), 
            out5=join(workpath,bam_dir,"{name}.Q5DD.bam"),
            out5f=join(workpath,bam_dir,"{name}.Q5DD.bam.flagstat"),
            out5i=join(workpath,bam_dir,"{name}.Q5DD.bam.idxstat"),
            out6=join(workpath,bam_dir,"{name}.bwa.Q5.duplic"), 
        params:
            rname='pl:dedup',
            picardver=config['bin'][pfamily]['tool_versions']['PICARDVER'],
            samtoolsver=config['bin'][pfamily]['tool_versions']['SAMTOOLSVER'],
            javaram='16g',
        shell: """
module load {params.samtoolsver};
module load {params.picardver}; 
java -Xmx{params.javaram} \
  -jar $PICARDJARPATH/picard.jar AddOrReplaceReadGroups \
  INPUT={input.bam2} \
  OUTPUT={output.out4} \
  TMP_DIR=/lscratch/$SLURM_JOBID \
  RGID=id \
  RGLB=library \
  RGPL=illumina \
  RGPU=machine \
  RGSM=sample; 
java -Xmx{params.javaram} \
  -jar $PICARDJARPATH/picard.jar MarkDuplicates \
  INPUT={output.out4} \
  OUTPUT={output.out5} \
  TMP_DIR=/lscratch/$SLURM_JOBID \
  VALIDATION_STRINGENCY=SILENT \
  REMOVE_DUPLICATES=true \
  METRICS_FILE={output.out6}
samtools index {output.out5}
samtools flagstat {output.out5} > {output.out5f}
samtools idxstats {output.out5} > {output.out5i}
            """

    rule insert_size:
        input:
            bam = lambda w : join(workpath,bam_dir,w.name + "." + w.ext + "." + extensions3[w.ext + "."])
        output:
            txt= join(workpath,qc_dir,"{name}.{ext}.insert_size_metrics.txt"),
            pdf= temp(join(workpath,qc_dir,"{name}.{ext}.insert_size_histogram.pdf")),
        params:
            rname="pl:insert_size",
            picardver=config['bin'][pfamily]['tool_versions']['PICARDVER'],
            javaram='16g',
        shell: """
module load {params.picardver};
java -Xmx{params.javaram} -jar $PICARDJARPATH/picard.jar CollectInsertSizeMetrics \
INPUT={input.bam} OUTPUT={output.txt} H={output.pdf}
"""

rule ppqt:
    input:
        bam = lambda w : join(workpath,bam_dir,w.name + "." + w.ext + "." + extensions3[w.ext + "."])
    output:
        ppqt= join(workpath,bam_dir,"{name}.{ext}.ppqt"),
        pdf= join(workpath,bam_dir,"{name}.{ext}.pdf"),
    params:
        rname="pl:ppqt",
        samtoolsver=config['bin'][pfamily]['tool_versions']['SAMTOOLSVER'],
        rver=config['bin'][pfamily]['tool_versions']['RVER'],
    run:
        commoncmd="module load {params.samtoolsver};module load {params.rver};"
        if se=="yes":
            cmd="Rscript Scripts/phantompeakqualtools/run_spp.R \
                -c={input.bam} -savp={output.pdf} -out={output.ppqt} -tmpdir=/lscratch/$SLURM_JOBID -rf;"
        elif pe=="yes":
            cmd="samtools view -b -f 66 -o /lscratch/$SLURM_JOBID/bam1.f66.bam {input.bam}; \
                samtools index /lscratch/$SLURM_JOBID/bam1.f66.bam; \
                Rscript Scripts/phantompeakqualtools/run_spp.R \
                -c=/lscratch/$SLURM_JOBID/bam1.f66.bam -savp={output.pdf} -out={output.ppqt} \
                -tmpdir=/lscratch/$SLURM_JOBID -rf;"
        shell(commoncmd+cmd)

rule ppqt_process:
    input:
        lambda w: [ join(workpath,bam_dir,sample + "." + w.ext + ".ppqt") for sample in samples ],
    output:
        out=join(workpath,bam_dir,'{ext}.ppqt.txt'),
    params:
        rname="pl:ppqt_process",
    run:
        o=open(output.out,'w')
        for inppqt in input:
            sample_name = inppqt.split("." + wildcards.ext)[0].split("/")[-1]
            if sample_name in uniq_inputs:
                o.write(sample_name + "\t" + "200" + "\n")
            else:
                file = list(map(lambda z:z.strip().split(),open(inppqt,'r').readlines()))
                ppqt_values = file[0][2].split(",")
                extenders = []
                for ppqt_value in ppqt_values:
                    if int(ppqt_value) > 150:
                        extenders.append(ppqt_value)
                if len(extenders) > 0:
                    o.write(sample_name + "\t" + extenders[0] + "\n")
                else:
                    print(sample_name)
                    print("All estimated fragments lengths were less than 150bp which will may cause the pipeline to fail.")
                    print("Potential causes include: wrong ref genome selected or low starting DNA.")
                    print("Assuming default estimated fragment length of 200bp.\n")
                    o.write(sample_name + "\t" + "200" + "\n")
        o.close()

rule bam2bw:
    input:
        bam=join(workpath,bam_dir,"{name}.{ext}.bam"),
        ppqt=join(workpath,bam_dir,"{ext}.ppqt.txt"),
    output:
        outbw=join(workpath,bw_dir,"{name}.{ext}.RPGC.bw"),
    params:
        rname="pl:bam2bw",
        reflen=config['references'][pfamily]['REFLEN'],
        deeptoolsver=config['bin'][pfamily]['tool_versions']['DEEPTOOLSVER'],
        effectivegenomesize=config['references'][pfamily]['EFFECTIVEGENOMESIZE'],
    run:
        lines=list(map(lambda x:x.strip().split("\t"),open(params.reflen).readlines()))
        genomelen=0
        chrs=[]
        includedchrs=[]
        excludedchrs=[]
        for chrom,l in lines:
            chrs.append(chrom)
            if not "_" in chrom and chrom!="chrX" and chrom!="chrM" and chrom!="chrY":
                includedchrs.append(chrom)
                genomelen+=int(l)
        excludedchrs=list(set(chrs)-set(includedchrs))
        commoncmd="module load {params.deeptoolsver};"
        cmd="bamCoverage --bam "+input.bam+" -o "+output.outbw+" --binSize 25 --smoothLength 75 --ignoreForNormalization "+" ".join(excludedchrs)+" --numberOfProcessors 32 --normalizeUsing RPGC --effectiveGenomeSize "+params.effectivegenomesize
        if pe=="yes":
            cmd+=" --centerReads"
        else:
            file = list(map(lambda z:z.strip().split(),open(input.ppqt,'r').readlines()))
            extender = [ ppqt[1] for ppqt in file if ppqt[0] == wildcards.name ] 
            cmd+=" -e "+extender[0]
        shell(commoncmd+cmd)

rule deeptools_prep:
    input:
        bw=expand(join(workpath,bw_dir,"{name}.{ext}.bw"),name=samples,ext=extensions),
        bam=expand(join(workpath,bam_dir,"{name}.{ext}.bam"),name=samples,ext=extensions2),
        bw2=expand(join(workpath,bw_dir,"{name}.Q5DD.RPGC.inputnorm.bw"),name=sampleswinput)
    output:
        temp(expand(join(workpath,bw_dir,"{ext}.deeptools_prep"),ext=extensions)),
        temp(expand(join(workpath,bam_dir,"{group}.{ext}.deeptools_prep"),group=groups,ext=extensions2)),
        temp(expand(join(workpath,bw_dir,"{group2}.{ext2}.deeptools_prep"), zip, group2=deepgroups,ext2=deepexts)),
    params:
        rname="pl:deeptools_prep",
        batch="--mem=10g --time=1:00:00",
    threads: 1
    run:
        for x in extensions2:
            bws=list(filter(lambda z:z.endswith(x+".RPGC.bw"),input.bw))
            bams=list(filter(lambda z:z.endswith(x+".bam"),input.bam))
            labels=list(map(lambda z:re.sub("."+x+".RPGC.bw","",z),
                list(map(lambda z:os.path.basename(z),bws))))
            bws2=list(filter(lambda z:z.endswith(x+".RPGC.inputnorm.bw"),input.bw2))
            labels2=list(map(lambda z:re.sub("."+x+".RPGC.inputnorm.bw","",z),
                list(map(lambda z:os.path.basename(z),bws2))))
            o=open(join(workpath,bw_dir,x+".RPGC.deeptools_prep"),'w')
            o.write("%s\n"%(x+".RPGC"))
            o.write("%s\n"%(" ".join(bws)))
            o.write("%s\n"%(" ".join(labels)))
            o.close()
            if len(bws2) > 0:
                for i in ["","prot."]:
                    o4=open(join(workpath,bw_dir,"InputNorm."+i+x+".RPGC.deeptools_prep"),'w')
                    o4.write("%s\n"%(x+".RPGC.inputnorm"))
                    o4.write("%s\n"%(" ".join(bws2)))
                    o4.write("%s\n"%(" ".join(labels2)))
                    o4.close()
            for group in groups:
                iter = [ i for i in range(len(labels)) if labels[i] in groupdatawinput[group] ]
                bws3 = [ bws[i] for i in iter ]
                labels3 = [ labels[i] for i in iter ]
                bams3 = [ bams[i] for i in iter ]
                for i in ["","prot."]:
                    o2=open(join(workpath,bam_dir,group+"."+i+x+".deeptools_prep"),'w')
                    o2.write("%s\n"%(x))
                    o2.write("%s\n"%(" ".join(bams3)))
                    o2.write("%s\n"%(" ".join(labels3)))
                    o2.close()
                for i in ["","prot."]:
                    o3=open(join(workpath,bw_dir,group+"."+i+x+".RPGC.deeptools_prep"),'w')
                    o3.write("%s\n"%(x+".RPGC"))
                    o3.write("%s\n"%(" ".join(bws3)))
                    o3.write("%s\n"%(" ".join(labels3)))
                    o3.close()

rule deeptools_QC:
    input:
        join(workpath,bw_dir,"{ext}.deeptools_prep"),
    output:
        heatmap=join(workpath,deeptools_dir,"spearman_heatmap.{ext}.pdf"),
        scatter=join(workpath,deeptools_dir,"spearman_scatterplot.{ext}.pdf"),
        pca=join(workpath,deeptools_dir,"pca.{ext}.pdf"),
	npz=temp(join(workpath,deeptools_dir,"{ext}.npz")),
    params:
        rname="pl:deeptools_QC",
        deeptoolsver=config['bin'][pfamily]['tool_versions']['DEEPTOOLSVER'],
        png=lambda w: join(workpath,deeptools_dir,"spearman_heatmap." + w.ext + "_mqc.png")
    run:
        import re
        commoncmd="module load {params.deeptoolsver}; module load python/2.7;"
        listfile=list(map(lambda z:z.strip().split(),open(input[0],'r').readlines()))
        ext=listfile[0][0]
        bws=listfile[1]
        labels=listfile[2]
        cmd1="multiBigwigSummary bins -b "+" ".join(bws)+" -l "+" ".join(labels)+" -out "+output.npz
        cmd2="plotCorrelation -in "+output.npz+" -o "+output.heatmap+" -c 'spearman' -p 'heatmap' --skipZeros --removeOutliers --plotNumbers"
	cmd3="plotCorrelation -in "+output.npz+" -o "+output.scatter+" -c 'spearman' -p 'scatterplot' --skipZeros --removeOutliers"
        cmd4="plotPCA -in "+output.npz+" -o "+output.pca
        cmd5="plotCorrelation -in "+output.npz+" -o "+ params.png +" -c 'spearman' -p 'heatmap' --skipZeros --removeOutliers --plotNumbers"
	shell(commoncmd+cmd1)
        shell(commoncmd+cmd2)
	shell(commoncmd+cmd3)
	shell(commoncmd+cmd4)
        if "Q5DD" in wildcards.ext:
            shell(commoncmd+cmd5)

rule deeptools_fingerprint:
    input:
        join(workpath,bam_dir,"{group}.sorted.deeptools_prep"),
    output:
        image=join(workpath,deeptools_dir,"{group}.fingerprint.sorted.pdf"),
        metrics=join(workpath,extra_fingerprint_dir,"{group}.fingerprint.metrics.sorted.tsv"),
    params:
        rname="pl:deeptools_fingerprint",
        deeptoolsver=config['bin'][pfamily]['tool_versions']['DEEPTOOLSVER'],
	nthreads="8"
    run:
        import re
        commoncmd="module load {params.deeptoolsver}; module load python/2.7;"
        listfile=list(map(lambda z:z.strip().split(),open(input[0],'r').readlines()))
        ext=listfile[0][0]
        bams=listfile[1]
        labels=listfile[2]
        cmd="plotFingerprint -b "+" ".join(bams)+" --labels "+" ".join(labels)+" -p "+params.nthreads+" --skipZeros --outQualityMetrics "+output.metrics+" --plotFile "+output.image 
        if se == "yes":
            cmd+=" -e 200"
        shell(commoncmd+cmd)

rule deeptools_fingerprint_Q5DD:
    input:
        join(workpath,bam_dir,"{group}.Q5DD.deeptools_prep")
    output:
        image=join(workpath,deeptools_dir,"{group}.fingerprint.Q5DD.pdf"),
        raw=temp(join(workpath,deeptools_dir,"{group}.fingerprint.raw.Q5DD.tab")),
        metrics=join(workpath,deeptools_dir,"{group}.fingerprint.metrics.Q5DD.tsv"),
    params:
        rname="pl:deeptools_fingerprint_Q5DD",
        deeptoolsver=config['bin'][pfamily]['tool_versions']['DEEPTOOLSVER'],
	nthreads="8"
    run:
        import re
        commoncmd="module load {params.deeptoolsver}; module load python/2.7;"
        listfile=list(map(lambda z:z.strip().split(),open(input[0],'r').readlines()))
        ext=listfile[0][0]
        bams=listfile[1]
        labels=listfile[2]
        cmd="plotFingerprint -b "+" ".join(bams)+" --labels "+" ".join(labels)+" -p "+params.nthreads+" --skipZeros --outQualityMetrics "+output.metrics+" --plotFile "+output.image+" --outRawCounts "+output.raw
        if se == "yes":
            cmd+=" -e 200"
        shell(commoncmd+cmd)

rule deeptools_genes:
    input:
        join(workpath,bw_dir,"{group}.{ext}.deeptools_prep")
    output:
        metaheat=join(workpath,deeptools_dir,"{group}.metagene_heatmap.{ext}.pdf"),
        TSSheat=join(workpath,deeptools_dir,"{group}.TSS_heatmap.{ext}.pdf"),
        metaline=join(workpath,deeptools_dir,"{group}.metagene_profile.{ext}.pdf"),
        TSSline=join(workpath,deeptools_dir,"{group}.TSS_profile.{ext}.pdf"),
        metamat=temp(join(workpath,deeptools_dir,"{group}.metagene.{ext}.mat.gz")),
        TSSmat=temp(join(workpath,deeptools_dir,"{group}.TSS.{ext}.mat.gz")),
        bed=temp(join(workpath,deeptools_dir,"{group}.geneinfo.{ext}.bed")),
    params:
        rname="pl:deeptools_genes",
        deeptoolsver=config['bin'][pfamily]['tool_versions']['DEEPTOOLSVER'],
        prebed=config['references'][pfamily]['GENEINFO'],
        nthreads="16"
    run:
        import re
        commoncmd="module load {params.deeptoolsver}; module load python/2.7;"
        listfile=list(map(lambda z:z.strip().split(),open(input[0],'r').readlines()))
        ext=listfile[0][0]
        bws=listfile[1]
        labels=listfile[2]
        if "prot" in wildcards.group:
            cmd1="grep --line-buffered 'protein_coding' "+ params.prebed  +" | awk -v OFS='\t' -F'\t' '{{print $1, $2, $3, $5, \".\", $4}}' > "+output.bed
        else:
            cmd1="awk -v OFS='\t' -F'\t' '{{print $1, $2, $3, $5, \".\", $4}}' "+params.prebed+" > "+output.bed
        cmd2="computeMatrix scale-regions -S "+" ".join(bws)+" -R "+output.bed+" -p "+params.nthreads+" --upstream 1000 --regionBodyLength 2000 --downstream 1000 --skipZeros -o "+output.metamat+" --samplesLabel "+" ".join(labels)
        cmd3="computeMatrix reference-point -S "+" ".join(bws)+" -R "+output.bed+" -p "+params.nthreads+" --referencePoint TSS --upstream 3000 --downstream 3000 --skipZeros -o "+output.TSSmat+" --samplesLabel "+" ".join(labels)
        cmd4="plotHeatmap -m "+output.metamat+" -out "+output.metaheat+" --colorMap 'PuOr_r' --yAxisLabel 'average RPGC' --regionsLabel 'genes' --legendLocation 'none'"
        cmd5="plotHeatmap -m "+output.TSSmat+" -out "+output.TSSheat+" --colorMap 'RdBu_r' --yAxisLabel 'average RPGC' --regionsLabel 'genes' --legendLocation 'none'"
        cmd6="plotProfile -m "+output.metamat+" -out "+output.metaline+" --plotHeight 15 --plotWidth 15 --perGroup --yAxisLabel 'average RPGC' --plotType 'se' --legendLocation upper-right"
        cmd7="plotProfile -m "+output.TSSmat+" -out "+output.TSSline+" --plotHeight 15 --plotWidth 15 --perGroup --yAxisLabel 'average RPGC' --plotType 'se' --legendLocation upper-left"
        shell(commoncmd+cmd1)
        shell(commoncmd+cmd2)
        shell(commoncmd+cmd3)
        shell(commoncmd+cmd4)
        shell(commoncmd+cmd5)
        shell(commoncmd+cmd6)
        shell(commoncmd+cmd7)

rule inputnorm:
    input:
        chip = join(workpath,bw_dir,"{name}.Q5DD.RPGC.bw"),
        ctrl = lambda w : join(workpath,bw_dir,chip2input[w.name] + ".Q5DD.RPGC.bw")
    output:
        join(workpath,bw_dir,"{name}.Q5DD.RPGC.inputnorm.bw")
    params:
        rname="pl:inputnorm",
        deeptoolsver=config['bin'][pfamily]['tool_versions']['DEEPTOOLSVER'],
    shell: """
module load {params.deeptoolsver};
bigwigCompare --binSize 25 --outFileName {output} --outFileFormat 'bigwig' --bigwig1 {input.chip} --bigwig2 {input.ctrl} --operation 'subtract' --skipNonCoveredRegions -p 32;
	"""

rule preseq:
    params:
        rname = "pl:preseq",
        preseqver=config['bin'][pfamily]['tool_versions']['PRESEQVER'],
    input:
        bam = join(workpath,bam_dir,"{name}.sorted.bam"),
    output:
        ccurve = join(workpath,qc_dir,"{name}.ccurve"),
    shell:
        """
module load {params.preseqver};
preseq c_curve -B -o {output.ccurve} {input.bam}            
        """
	
rule NRF:
    input:
        bam=join(workpath,bam_dir,"{name}.sorted.bam"),
    params:
        rname='pl:NRF',
        samtoolsver=config['bin'][pfamily]['tool_versions']['SAMTOOLSVER'],
        rver=config['bin'][pfamily]['tool_versions']['RVER'],
        preseqver=config['bin'][pfamily]['tool_versions']['PRESEQVER'],
        nrfscript=join(workpath,"Scripts","atac_nrf.py "),            
    output:
        preseq=join(workpath,qc_dir,"{name}.preseq.dat"),
        preseqlog=join(workpath,qc_dir,"{name}.preseq.log"),
        nrf=temp(join(workpath,qc_dir,"{name}.nrf")),
    threads: 16
    shell: """
module load {params.preseqver};
preseq lc_extrap -P -B -D -o {output.preseq} {input.bam} -seed 12345 -v -l 100000000000 2> {output.preseqlog}
python {params.nrfscript} {output.preseqlog} > {output.nrf}
        """

rule rawfastqc:
    input:
        expand(join(workpath,"{name}.R1.fastq.gz"), name=samples) if \
            se == "yes" else \
            expand(join(workpath,"{name}.R{rn}.fastq.gz"), name=samples,rn=[1,2])
    output:
        folder=join(workpath,'rawfastQC'),
        html=expand(join(workpath,'rawfastQC',"{name}.R1_fastqc.html"),name=samples)
    params:
        rname='pl:rawfastqc',
        batch='--cpus-per-task=32 --mem=100g --time=48:00:00',
        fastqcver=config['bin'][pfamily]['tool_versions']['FASTQCVER']
    threads: 32
    shell: """
if [ ! -d {output} ]; then
mkdir {output.folder};
fi
module load {params.fastqcver};
fastqc {input} -t {threads} -o {output.folder}
            """

rule fastqc:
    params:
        rname='pl:fastqc',
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        fastqcver=config['bin'][pfamily]['tool_versions']['FASTQCVER']
    input:
        expand(join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),name=samples) if \
            se == "yes" else \
            expand(join(workpath,trim_dir,"{name}.R{rn}.trim.fastq.gz"), name=samples,rn=[1,2])
    output:
        folder=join(workpath,'fastQC'),
        html=expand(join(workpath,'fastQC',"{name}.R1.trim_fastqc.html"),name=samples)
    threads: 32
    shell: """
mkdir -p {output.folder};
module load {params.fastqcver};
fastqc {input} -t {threads} -o {output.folder}
            """

rule fastq_screen:
    input:
        join(workpath,trim_dir,"{name}.R1.trim.fastq.gz") if se == "yes" else \
            expand(join(workpath,trim_dir,"{name}.R{rn}.trim.fastq.gz"),name=samples,rn=[1,2])
    output:
        join(workpath,"FQscreen","{name}.R1.trim_screen.txt") if se == "yes" else \
            expand(join(workpath,"FQscreen","{name}.R{rn}.trim_screen.txt"),name=samples,rn=[1,2]),
        join(workpath,"FQscreen","{name}.R1.trim_screen.png") if se == "yes" else \
            expand(join(workpath,"FQscreen","{name}.R{rn}.trim_screen.png"),name=samples,rn=[1,2]),
        join(workpath,"FQscreen2","{name}.R1.trim_screen.txt") if se == "yes" else \
            expand(join(workpath,"FQscreen2","{name}.R{rn}.trim_screen.txt"),name=samples,rn=[1,2]),
        join(workpath,"FQscreen2","{name}.R1.trim_screen.png") if se == "yes" else \
            expand(join(workpath,"FQscreen2","{name}.R{rn}.trim_screen.png"),name=samples,rn=[1,2]),
    params:
        rname='pl:fqscreen',
        bowtie2ver=config['bin'][pfamily]['tool_versions']['BOWTIE2VER'],
        perlver=config['bin'][pfamily]['tool_versions']['PERLVER'],
        fastq_screen=config['bin'][pfamily]['FASTQ_SCREEN'],
        fastq_screen_config=config['references'][pfamily]['FASTQ_SCREEN_CONFIG'],
        fastq_screen_config2=config['references'][pfamily]['FASTQ_SCREEN_CONFIG2'],
        outdir = "FQscreen",
        outdir2 = "FQscreen2",
    threads: 24
    shell: """
module load {params.bowtie2ver} ;
module load {params.perlver};
{params.fastq_screen} --conf {params.fastq_screen_config} \
    --outdir {params.outdir} --subset 1000000 \
    --aligner bowtie2 --force {input}
{params.fastq_screen} --conf {params.fastq_screen_config2} \
    --outdir {params.outdir2} --subset 1000000 \
    --aligner bowtie2 --force {input}
            """

rule ngsqc:
    input:
        tagAlign=join(workpath,bam_dir,"{name}.Q5DD.tagAlign.gz") if \
          se == "yes" else \
          join(workpath,bam_dir,"{name}.Q5DD.bam")
    output:
        file=join(workpath,qc_dir,"{name}.Q5DD.NGSQC_report.txt"),
    params:
        rname="pl:ngsqc",
        bedtoolsver=config['bin'][pfamily]['tool_versions']['BEDTOOLSVER'],
        ngsqc=config['bin'][pfamily]['NGSQC'],
        genomefile=config['references'][pfamily]['REFLEN']
    shell: """
module load {params.bedtoolsver};
if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID; fi
cd /lscratch/$SLURM_JOBID;
# if se == "yes"
if [ {se} == "yes" ]; then
    cp {input.tagAlign} tmptagAlign.gz
    gzip -d tmptagAlign.gz
else
    bamToBed -i {input.tagAlign} > tmptagAlign
fi
{params.ngsqc} -v -o tmpOut tmptagAlign {params.genomefile}
mv tmpOut/NGSQC_report.txt {output.file}
"""

rule ngsqc_plot:
    input:
        ngsqc=expand(join(workpath,"QC","{name}.Q5DD.NGSQC_report.txt"),name=samples),
    output:
        png=expand(join(workpath,"QC","{group}.NGSQC.Q5DD.png"),group=groups),
    params:
        rname="pl:ngsqc_plot",
        script=join(workpath,"Scripts","ngsqc_plot.py"),
    run:
        commoncmd1 = "if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID ;fi"
        commoncmd2 = "cd /lscratch/$SLURM_JOBID;"
        cmd1 = ""
        cmd2 = ""
        cmd3 = ""
        cmd4 = ""
        sample = []
        for group in groups:
            labels = [ sample for sample in samples if sample in groupdatawinput[group] ]
            cmd1 = cmd1 + "mkdir " + group + "; "
            for label in labels:
                cmd1 = cmd1 + "cp " + workpath + "/QC/" + label + ".Q5DD.NGSQC_report.txt " + group + "; " 
                if label not in sample:
                    cmd4 = cmd4 + "mv " + group + "/" + label + ".Q5DD.NGSQC.txt " + workpath + "/QC; "
                    sample.append(label)
            cmd2 = cmd2 + "python " + params.script + " -d '" + group + "' -e 'Q5DD' -g '" + group + "'; "
            cmd3 = cmd3 + "mv " + group + "/" + group + ".NGSQC.Q5DD.png " + workpath + "/QC/" + group + ".NGSQC.Q5DD.png" + "; "
        shell(commoncmd1)
        shell(commoncmd2 + cmd1 + cmd2 + cmd3 + cmd4)

rule QCstats:
    input:
        flagstat=join(workpath,bam_dir,"{name}.sorted.bam.flagstat"),
        infq=join(workpath,"{name}.R1.fastq.gz"),	
        ddflagstat=join(workpath,bam_dir,"{name}.Q5DD.bam.flagstat"),
        nrf=join(workpath,qc_dir,"{name}.nrf"),
        ppqt=join(workpath,bam_dir,"{name}.Q5DD.ppqt"),
        ppqt2=join(workpath,bam_dir,"Q5DD.ppqt.txt"),
    params:
        rname='pl:QCstats',
        filterCollate=join(workpath,"Scripts","filterMetrics"),   
    output:
        sampleQCfile=temp(join(workpath,qc_dir,"{name}.qcmetrics")),
    threads: 16
    shell: """
# Number of reads
#grep 'in total' {input.flagstat} | awk '{{print $1,$3}}' | {params.filterCollate} {wildcards.name} tnreads > {output.sampleQCfile}
zcat {input.infq} | wc -l | {params.filterCollate} {wildcards.name} tnreads > {output.sampleQCfile}
# Number of mapped reads
grep 'mapped (' {input.flagstat} | awk '{{print $1,$3}}' | {params.filterCollate} {wildcards.name} mnreads >> {output.sampleQCfile}
# Number of uniquely mapped reads
grep 'mapped (' {input.ddflagstat} | awk '{{print $1,$3}}' | {params.filterCollate} {wildcards.name} unreads >> {output.sampleQCfile}
# NRF, PCB1, PCB2
cat {input.nrf} | {params.filterCollate} {wildcards.name} nrf >> {output.sampleQCfile}
# NSC, RSC, Qtag
awk '{{print $(NF-2),$(NF-1),$NF}}' {input.ppqt} | {params.filterCollate} {wildcards.name} ppqt >> {output.sampleQCfile}
# Fragment Length
fragLen=`grep {wildcards.name} {input.ppqt2} | cut -f 2`
echo "{wildcards.name}\tFragmentLength\t$fragLen" >> {output.sampleQCfile}
            """

rule QCTable:
    input:
        expand(join(workpath,qc_dir,"{name}.qcmetrics"), name=samples),
    params:
        rname='pl:QCTable',
        inputstring=" ".join(expand(join(workpath,qc_dir,"{name}.qcmetrics"), name=samples)),
        filterCollate=join(workpath,"Scripts","createtable"),
    output:
        qctable=join(workpath,qc_dir,"QCTable.txt"),
    threads: 16
    shell: """
cat {params.inputstring} | {params.filterCollate} > {output.qctable}
            """

rule multiqc:
    input: 
        expand(join(workpath,"FQscreen","{name}.R1.trim_screen.txt"),name=samples),
        expand(join(workpath,"FQscreen2","{name}.R1.trim_screen.txt"),name=samples),
        expand(join(workpath,qc_dir,"{name}.ccurve"), name=samples),
        expand(join(workpath,bam_dir,"{name}.Q5DD.bam.flagstat"), name=samples),
        expand(join(workpath,bam_dir,"{name}.Q5.bam.flagstat"), name=samples),
        join(workpath,qc_dir,"QCTable.txt"),
        expand(join(workpath,'rawfastQC',"{name}.R1_fastqc.html"),name=samples),
        expand(join(workpath,'fastQC',"{name}.R1.trim_fastqc.html"),name=samples),
        expand(join(workpath,"QC","{group}.NGSQC.Q5DD.png"),group=groups),
        expand(join(workpath,deeptools_dir,"{group}.fingerprint.raw.Q5DD.tab"),group=groups),
        join(workpath,deeptools_dir,"spearman_heatmap.Q5DD.RPGC.pdf"),
    output:
        join(workpath,"Reports","multiqc_report.html")
    params:
        rname="pl:multiqc",
        multiqc=config['bin'][pfamily]['tool_versions']['MULTIQCVER'],
	qcconfig=config['bin'][pfamily]['CONFMULTIQC'],
	dir=join("..",extra_fingerprint_dir)
    shell: """
module load {params.multiqc}
cd Reports && multiqc -f -c {params.qcconfig} --interactive -e cutadapt --ignore {params.dir} -d ../
"""
