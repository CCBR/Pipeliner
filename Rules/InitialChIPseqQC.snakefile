from snakemake.utils import R
from os.path import join

configfile: "run.json"
    
workpath = config['project']['workpath']    
filetype = config['project']['filetype']
readtype = config['project']['readtype']

trim_dir='trim'

#print(samples)

if readtype == 'Single' :
    rule InitialChIPseqQC:
        params: 
            batch='--time=168:00:00'
        input: 
            "Reports/multiqc_report.html",
            "rawQC",
            "QC",
            expand("FQscreen/{name}.R1_screen.txt",name=samples),
            expand("FQscreen/{name}.R1_screen.png",name=samples),
            expand("{name}.sorted.rmdup.bam.bai", name=samples),
            expand("{name}.sorted.rmdup.bam", name=samples),
            expand("{name}.shifts", name=samples),
            expand("{name}.rmdup.shifts", name=samples),
            expand(join(trim_dir,'{name}.R1.trim.fastq.gz'), name=samples)

                   

    rule fastq_screen:
        input: 
            expand("{name}.R1.fastq.gz", name=samples)
        output:
            "FQscreen/{name}.R1_screen.txt",
            "FQscreen/{name}.R1_screen.png",
        params: 
            rname='pl:fqscreen',
            fastq_screen=config['bin'][pfamily]['FASTQ_SCREEN'],
            outdir = "FQscreen",
            config=config['references'][pfamily]['FASTQ_SCREEN_CONFIG'] 
        threads: 24
        shell:
            """
            module load bowtie/2-2.2.9 ;
            module load perl/5.18.2; 
            {params.fastq_screen} --conf {params.config} \
                --outdir {params.outdir} --subset 1000000 \
                --aligner bowtie2 --force {input}
            """

    rule rawfastqc:
        input: 
            expand("{name}.R1.fastq.gz", name=samples) 
        output: 
            'rawQC'
        priority: 2
        params: 
            rname='pl:rawfastqc',
            batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
            fastqcver=config['bin'][pfamily]['FASTQCVER']
        threads: 32
        shell: 
            """
            mkdir -p {output};
            module load {params.fastqcver}; 
            fastqc {input} -t {threads} -o {output}
            """

    # rule trimgalore:
    #     input:
    #         "{name}.R1.fastq.gz"
    #     params:
    #         rname='pl:trimgalore',
    #         d = trim_dir,
    #         of = '{d}/{name}.R1_trimmed.fq.gz'
    #     output:
    #         of1 = '{d}/{name}.R1.trim.fastq.gz',
    #         of2 = '{d}/{name}.R1.fastq.gz_trimming_report.txt',
    #     shell:
    #         """
    #         module load trimgalore; 
    #         trim_galore -o {params.d} {input};
    #         mv {params.of} {output.of1}
    #         """

    rule trim:
        input:
            if1 = "{name}.R1.fastq.gz"
        params:
            rname='pl:trimmomatic_se',
            d=trim_dir,
            batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
            trimmomaticver=config['bin'][pfamily]['TRIMMOMATICVER'],
            fastawithadaptersetc=config['references'][pfamily]['FASTAWITHADAPTERSETC'],
            seedmismatches=config['bin'][pfamily]['SEEDMISMATCHES'],
            palindromeclipthreshold=config['bin'][pfamily]['PALINDROMECLIPTHRESHOLD'],
            simpleclipthreshold=config['bin'][pfamily]['SIMPLECLIPTHRESHOLD'],
            leadingquality=config['bin'][pfamily]['LEADINGQUALITY'],
            trailingquality=config['bin'][pfamily]['TRAILINGQUALITY'],
            minlen=config['bin'][pfamily]['MINLEN'],
        output:
            of1 = '{d}/{name}.R1.trim.fastq.gz',
            of2 = '{d}/{name}.R1.fastq.gz_trimmomatic.err',
        threads: 32
        shell:
            """
            module load {params.trimmomaticver};
            java -classpath $TRIMMOJAR org.usadellab.trimmomatic.TrimmomaticSE -threads {threads} {input.if1} {output.of1} ILLUMINACLIP:{param.fastawithadaptersetc}:{param.seedmismatches}:{param.palindromeclipthreshold}:{param.simpleclipthreshold} LEADING:{param.leadingquality} TRAILING:{param.trailingquality} MINLEN:{param.minlen} 2> {output.of2}
            """            

    rule fastqc:  
        input:
            expand("trim/{name}.R1.trim.fastq.gz", name=samples),
        output: "QC"
        priority: 2
        params: 
            rname='pl:fastqc',
            batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
            fastqcver=config['bin'][pfamily]['FASTQCVER']
        threads: 32
        shell: 
            """
            mkdir -p {output};
            module load {params.fastqcver}; 
            fastqc {input} -t {threads} -o {output}
            """

    rule BWA:
        input:
            "trim/{name}.R1.trim.fastq.gz"
        output:
            "{name}.sorted.bam"
        params:
            rname='pl:bwa',
            reference= config['references'][pfamily]['BWA'],
        threads: 32
        shell: 
            """
            module load bwa;
            module load samtools;
            bwa mem -t {threads} {params.reference} {input} | \
                samtools sort -o {output}    
            """
                
        '''
        shell: 
            """
            module load bwa;
            module load samtools;
            bwa mem -t  {threads} {params.reference} {input} | \
            awk -F '\t' \
            '/^@SQ/{{\
                OFS="\t"; \
                split($2,a,":"); \
                if(a[2]=="MT"){{ \
                    $2 = a[1]":chrM"; }} \
                else {{ $2=a[1]":chr"a[2]}}    \
                print; next;}} \
            /^@/{{print; next}} \
            {{OFS="\t"; \
                if($3=="*"){{}} \
                else if($3=="MT"){{$3="chrM";}}\
                else if($3!~/chr/){{$3="chr"$3;}}; \
                if($7=="*"){{}} \
                else if($7=="MT"){{$7="chrM";}}\
                else if($7!~/chr/){{$7="chr"$7;}}; \
            print}}' | samtools sort -o {output}
            """  
        '''

    rule picard:
        input: 
            file1= "{name}.sorted.bam"
        output:
            outstar1=temp("{name}.bwa_rg_added.sorted.bam"), 
            outstar2="{name}.bwa_rg_added.sorted.dmark.bam",
            outstar3="{name}.bwa.duplic" 
        params:
            rname='pl:picard',
            batch='--mem=24g --time=10:00:00 --gres=lscratch:800',
            picardver=config['bin'][pfamily]['PICARDVER']
        shell: 
             """
             module load {params.picardver}; 
             java -Xmx10g \
                  -jar $PICARDJARPATH/AddOrReplaceReadGroups.jar \
                  INPUT={input.file1} \
                  OUTPUT={output.outstar1} \
                  TMP_DIR=/lscratch/$SLURM_JOBID \
                  RGID=id \
                  RGLB=library \
                  RGPL=illumina \
                  RGPU=machine \
                  RGSM=sample; 
             java -Xmx10g \
                  -jar $PICARDJARPATH/MarkDuplicates.jar \
                  INPUT={output.outstar1} \
                  OUTPUT={output.outstar2} \
                  TMP_DIR=/lscratch/$SLURM_JOBID \
                  CREATE_INDEX=true \
                  VALIDATION_STRINGENCY=SILENT \
                  METRICS_FILE={output.outstar3}
             """
                
    rule rmdup:
        input:
            file1= "{name}.sorted.bam"
        output:
            rmdup_bam="{name}.sorted.rmdup.bam", 
            rmdup_bai="{name}.sorted.rmdup.bam.bai" 
        params:
            rname='pl:rmdup'
        shell:
            """
            module load samtools;
            samtools rmdup {input.file1} {output.rmdup_bam};
            samtools index {output.rmdup_bam};
            """
    
    rule shiftstats:
        input: 
            if1 = "{name}.sorted.bam",
            if2 = "{name}.sorted.rmdup.bam" 
        output:
            of1 = "{name}.shifts",
            of2 = "{name}.rmdup.shifts"
        params:
            rname='pl:shiftstats',
            batch='--mem=24g --time=10:00:00 --gres=lscratch:800'
        shell: 
             """
             touch {output.of1}
             touch {output.of2}
             """

    rule stats:
        input:
            file1= "{name}.bwa_rg_added.sorted.dmark.bam"
        output:
            outstar2="{name}.flagstat.concord.txt" 
        params: 
            rname='pl:stats',
            batch='--mem=24g --time=10:00:00 --gres=lscratch:800',
            picardver=config['bin'][pfamily]['PICARDVER'],
        shell:
            """
            module load samtools; 
            samtools flagstat {input.file1} > {output.outstar2}; 
            echo 0 >> {output.outstar2};
            echo 0 >> {output.outstar2};
            #samtools view -f 0x2 {input.file1} | wc -l >>{output.outstar2}; 
            #samtools view {input.file1} | grep -w -c NH:i:1  >>{output.outstar2}
            """
            
    rule multiqc:
        input: 
            expand("FQscreen/{name}.R1_screen.png",name=samples),
            expand("{name}.flagstat.concord.txt",name=samples),
            #expand("{name}.RnaSeqMetrics.txt",name=samples),
            #rules.picard.output,
            expand("{name}.bwa.duplic", name=samples),
            #rules.fastqc.output,
            expand(join(trim_dir,'{name}.R1.trim.fastq.gz'), name=samples),
            #rules.trimgalore.output,
            expand("FQscreen/{name}.R1_screen.txt",name=samples),
            #rules.fastq_screen.output,
            #rules.stats.output,

            
        output:
            "Reports/multiqc_report.html"
        params:
            rname="pl:multiqc",
            multiqc=config['bin'][pfamily]['MULTIQC'],
        threads: 1
        shell:  """
                module load {params.multiqc}
                cd Reports && multiqc -f ../
                """

elif readtype == 'Paired' :
    rule all:
        params: 
            batch='--time=168:00:00'
        input: 
            config['project']['id']+"_"+config['project']['flowcellid']+".xlsx",
            "Reports/multiqc_report.html",
            expand("{name}.RnaSeqMetrics.txt",name=samples),
            "rawQC",
            "QC",
            expand("FQscreen/{name}.R1_screen.txt",name=samples),
            expand("FQscreen/{name}.R1_screen.png",name=samples),
            expand("FQscreen/{name}.R2_screen.txt",name=samples),
            expand("FQscreen/{name}.R2_screen.png",name=samples),
            expand("{name}.InsertSizeMetrics.txt",name=samples) 


    rule fastq_screen:
        input: 
            "{name}.R1.fastq.gz",
            "{name}.R2.fastq.gz"
        output:
            "FQscreen/{name}.R1_screen.txt",
            "FQscreen/{name}.R1_screen.png",
            "FQscreen/{name}.R2_screen.txt",
            "FQscreen/{name}.R2_screen.png" 
        params: 
            rname='pl:fqscreen',
            fastq_screen=config['bin'][pfamily]['FASTQ_SCREEN'],
            outdir = "FQscreen",
            config=config['references'][pfamily]['FASTQ_SCREEN_CONFIG'] 
        threads: 24
        shell:
            """
            module load bowtie ; 
            {params.fastq_screen} --conf {params.config} \
                --outdir {params.outdir} --subset 1000000 \
                --aligner bowtie2 --force {input}
            """

    rule rawfastqc:
        input: 
            expand("{name}.R1.fastq.gz", name=samples), 
            expand("{name}.R2.fastq.gz", name=samples)
        output: 
            "rawQC"
        priority: 2
        params: 
            rname='pl:rawfastqc',
            batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
            fastqcver=config['bin'][pfamily]['FASTQCVER']
        threads: 32
        shell: 
            """
            mkdir -p {output};
            module load {params.fastqcver}; 
            fastqc {input} -t {threads} -o {output}
            """


    rule trimmomatic_pe:
        input: 
            file1= join(workpath, "{name}.R1."+filetype),
            file2= join(workpath, "{name}.R2."+filetype),
        output:
            out11="trim/{name}_R1_001_trim_paired.fastq.gz",
            out12="trim/{name}_R1_001_trim_unpaired.fastq.gz",
            out21="trim/{name}_R2_001_trim_paired.fastq.gz",
            out22="trim/{name}_R2_001_trim_unpaired.fastq.gz",
            err="QC/{name}_run_trimmomatic.err"
        params:
            rname='pl:trimmomatic_pe',
            batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
            trimmomaticver=config['bin'][pfamily]['TRIMMOMATICVER'],
            fastawithadaptersetc=config['references'][pfamily]['FASTAWITHADAPTERSETC'],
            seedmismatches=config['bin'][pfamily]['SEEDMISMATCHES'],
            palindromeclipthreshold=config['bin'][pfamily]['PALINDROMECLIPTHRESHOLD'],
            simpleclipthreshold=config['bin'][pfamily]['SIMPLECLIPTHRESHOLD'],
            leadingquality=config['bin'][pfamily]['LEADINGQUALITY'],
            trailingquality=config['bin'][pfamily]['TRAILINGQUALITY'],
            windowsize=config['bin'][pfamily]['WINDOWSIZE'],
            windowquality=config['bin'][pfamily]['WINDOWQUALITY'],
            targetlength=config['bin'][pfamily]['TARGETLENGTH'],
            strictness=config['bin'][pfamily]['STRICTNESS'],
            minlen=config['bin'][pfamily]['MINLEN'],
            headcroplength=config['bin'][pfamily]['HEADCROPLENGTH']
        threads:32
        shell:
            """
            module load {params.trimmomaticver}; 
            java -classpath $TRIMMOJAR   org.usadellab.trimmomatic.TrimmomaticPE \
                 -threads {threads} {input.file1} {input.file2} \
                          {output.out11} {output.out12} {output.out21} {output.out22} \
                 ILLUMINACLIP:{params.fastawithadaptersetc}:{params.seedmismatches}:{params.palindromeclipthreshold}:{params.simpleclipthreshold} 2> {output.err}
            """

    rule fastqc:  
        input:
            expand("trim/{name}_R1_001_trim_paired.fastq.gz", name=samples),
            expand("trim/{name}_R2_001_trim_paired.fastq.gz", name=samples)  
        output: "QC"
        priority: 2
        params: 
            rname='pl:fastqc',
            batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
            fastqcver=config['bin'][pfamily]['FASTQCVER']
        threads: 32
        shell: 
            """
            mkdir -p {output};
            module load {params.fastqcver}; 
            fastqc {input} -t {threads} -o {output}
            """

    rule bwa:
        input:
            file1= "trim/{name}_R1_001_trim_paired.fastq.gz",
            file2="trim/{name}_R2_001_trim_paired.fastq.gz"
        output:
            out= "{name}.p2.Aligned.sortedByCoord.out.bam"
        params:
            rname='pl:bwa',
            prefix="{name}",
            batch='--cpus-per-task=32 --mem=32g --time=48:00:00',
            reference= config['references'][pfamily]['BWA'],
        threads: 32
        shell: 
            """
            module load bwa;
            bwa mem -t  {threads} {params.reference} {file1} {file2} | samtools sort -o {out}
            """  

    rule picard:
        input: 
            file1= "{name}.p2.Aligned.sortedByCoord.out.bam"
        output:
            outstar1=temp("{name}.star_rg_added.sorted.bam"), 
            outstar2="{name}.star_rg_added.sorted.dmark.bam",
            outstar3="{name}.star.duplic" 
        params:
            rname='pl:picard',
            batch='--mem=24g --time=10:00:00 --gres=lscratch:800',
            picardver=config['bin'][pfamily]['PICARDVER']
        shell: 
             """
             module load {params.picardver}; 
             java -Xmx10g \
                  -jar $PICARDJARPATH/AddOrReplaceReadGroups.jar \
                  INPUT={input.file1} \
                  OUTPUT={output.outstar1} \
                  TMP_DIR=/lscratch/$SLURM_JOBID \
                  RGID=id \
                  RGLB=library \
                  RGPL=illumina \
                  RGPU=machine \
                  RGSM=sample; 
             java -Xmx10g \
                  -jar $PICARDJARPATH/MarkDuplicates.jar \
                  INPUT={output.outstar1} \
                  OUTPUT={output.outstar2} \
                  TMP_DIR=/lscratch/$SLURM_JOBID \
                  CREATE_INDEX=true \
                  VALIDATION_STRINGENCY=SILENT \
                  METRICS_FILE={output.outstar3}
             """

    rule stats:
        input:
            file1= "{name}.star_rg_added.sorted.dmark.bam"
        output:
            outstar1="{name}.RnaSeqMetrics.txt",
            outstar2="{name}.flagstat.concord.txt", 
            outstar3="{name}.InsertSizeMetrics.txt", 
            outstar4="{name}.InsertSizeHisto.pdf"
        params: 
            rname='pl:stats',
            batch='--mem=24g --time=10:00:00 --gres=lscratch:800',
            picardver=config['bin'][pfamily]['PICARDVER'],
        shell:
            """
            module load R;
            module load {params.picardver};
            java -Xmx10g \
                 -jar $PICARDJARPATH/CollectRnaSeqMetrics.jar \
                 INPUT={input.file1} \
                 OUTPUT={output.outstar1} \
                 TMP_DIR=/lscratch/$SLURM_JOBID  \
                 VALIDATION_STRINGENCY=SILENT ; 
            java -Xmx10g \
                 -jar $PICARDJARPATH/CollectInsertSizeMetrics.jar \
                 INPUT={input.file1} \
                 OUTPUT={output.outstar3} \
                 HISTOGRAM_FILE={output.outstar4} \
                 MINIMUM_PCT=0.5 \
                 TMP_DIR=/lscratch/$SLURM_JOBID ;

            module load samtools; 
            samtools flagstat {input.file1} > {output.outstar2}; 
            samtools view -f 0x2 {input.file1} | wc -l >>{output.outstar2}; 
            samtools view {input.file1} | grep -w -c NH:i:1  >>{output.outstar2}
            """


    rule rnaseq_multiqc:
        input: 
            expand("{name}.Rdist.info",name=samples),
            expand("FQscreen/{name}.R1_screen.png",name=samples),
            expand("FQscreen/{name}.R2_screen.png",name=samples),
            expand("{name}.flagstat.concord.txt",name=samples),
            expand("{name}.RnaSeqMetrics.txt",name=samples),
            expand("{name}.InsertSizeMetrics.txt",name=samples)
        output:
            "Reports/multiqc_report.html"
        params:
            rname="pl:multiqc",
            multiqc=config['bin'][pfamily]['MULTIQC'],
        threads: 1
        shell:  """
                module load {params.multiqc}
                cd Reports && multiqc -f  ../
                """

    rule RNAseq_generate_QC_table:
        input:
            expand("QC/{name}_run_trimmomatic.err",name=samples),
            expand("{name}.star.duplic",name=samples),
            expand("{name}.p2.Log.final.out",name=samples),
            expand("{name}.RnaSeqMetrics.txt",name=samples)
        output:
            config['project']['id']+"_"+config['project']['flowcellid']+".xlsx"
        params:
            workpath,
            project=config['project']['id'],
            flowcell=config['project']['flowcellid'],
            rname="pl:QC_table"
        shell: 
            """
            perl Scripts/CollectPipelineStats2Tab_v2.3.pl \
                -p {params.project}\
                -f {params.flowcell}\
                -d {params.workpath}\
                -r 5\
                -e 2;
            perl Scripts/Tab2Excel_v2.3.pl -i {params.project}_{params.flowcell} -r 5
            """


