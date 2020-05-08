rule pindel_conf:
  params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindelconf"
  input: expand("{x}.recal.bam", x=samples)
  output: "pindel_out/pindel_config",
  shell: "mkdir -p pindel_out; module load perl/5.18.4; perl Scripts/make_pindel_config.pl"

if config['project']['annotation'] == "hg19":
  rule pindel_germline_1:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr1_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 1 -f {params.genome} -o pindel_out/pindel_calls_chr1 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr1 -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr1 -d 01022009 -G -c 1"

  rule pindel_germline_2:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr2_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 2 -f {params.genome} -o pindel_out/pindel_calls_chr2 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr2.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr2 -d 01022009 -G -c 2"

  rule pindel_germline_3:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr3_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 3 -f {params.genome} -o pindel_out/pindel_calls_chr3 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr3.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr3 -d 01022009 -G -c 3"

  rule pindel_germline_4:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr4_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 4 -f {params.genome} -o pindel_out/pindel_calls_chr4 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr4.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr4 -d 01022009 -G -c 4"

  rule pindel_germline_5:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr5_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 5 -f {params.genome} -o pindel_out/pindel_calls_chr5 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr5.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr5 -d 01022009 -G -c 5"

  rule pindel_germline_6:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr6_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 6 -f {params.genome} -o pindel_out/pindel_calls_chr6 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr6.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr6 -d 01022009 -G -c 6"

  rule pindel_germline_7:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr7_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 7 -f {params.genome} -o pindel_out/pindel_calls_chr7 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr7.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr7 -d 01022009 -G -c 7"

  rule pindel_germline_8:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr8_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 8 -f {params.genome} -o pindel_out/pindel_calls_chr8 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr8.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr8 -d 01022009 -G -c 8"

  rule pindel_germline_9:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr9_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 9 -f {params.genome} -o pindel_out/pindel_calls_chr9 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr9.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr9 -d 01022009 -G -c 9"

  rule pindel_germline_10:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr10_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 10 -f {params.genome} -o pindel_out/pindel_calls_chr10 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr10.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr10 -d 01022009 -G -c 10"

  rule pindel_germline_11:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr11_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 11 -f {params.genome} -o pindel_out/pindel_calls_chr11 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr11.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr11 -d 01022009 -G -c 11"

  rule pindel_germline_12:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr12_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 12 -f {params.genome} -o pindel_out/pindel_calls_chr12 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr12.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr12 -d 01022009 -G -c 12"

  rule pindel_germline_13:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr13_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 13 -f {params.genome} -o pindel_out/pindel_calls_chr13 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr13.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr13 -d 01022009 -G -c 13"

  rule pindel_germline_14:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr14_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 14 -f {params.genome} -o pindel_out/pindel_calls_chr14 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr14.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr14 -d 01022009 -G -c 14"

  rule pindel_germline_15:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr15_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 15 -f {params.genome} -o pindel_out/pindel_calls_chr15 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr15.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr15 -d 01022009 -G -c 15"

  rule pindel_germline_16:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr16_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 16 -f {params.genome} -o pindel_out/pindel_calls_chr16 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr16.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr16 -d 01022009 -G -c 16"

  rule pindel_germline_17:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr17_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 17 -f {params.genome} -o pindel_out/pindel_calls_chr17 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr17.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr17 -d 01022009 -G -c 17"

  rule pindel_germline_18:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr18_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 18 -f {params.genome} -o pindel_out/pindel_calls_chr18 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr18.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr18 -d 01022009 -G -c 18"

  rule pindel_germline_19:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr19_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 19 -f {params.genome} -o pindel_out/pindel_calls_chr19 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr19.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr19 -d 01022009 -G -c 19"

  rule pindel_germline_20:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr20_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 20 -f {params.genome} -o pindel_out/pindel_calls_chr20 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr20.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr20 -d 01022009 -G -c 20"

  rule pindel_germline_21:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr21_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 21 -f {params.genome} -o pindel_out/pindel_calls_chr21 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr21.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr21 -d 01022009 -G -c 21"

  rule pindel_germline_22:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr22_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c 22 -f {params.genome} -o pindel_out/pindel_calls_chr22 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr22.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chr22 -d 01022009 -G -c 22"

  rule pindel_germline_X:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chrX_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c X -f {params.genome} -o pindel_out/pindel_calls_chrX --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chrX.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chrX -d 01022009 -G -c X"

  rule pindel_germline_Y:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chrY_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c Y -f {params.genome} -o pindel_out/pindel_calls_chrY --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chrY.vcf -r {params.genome} -R hs37d5 -P pindel_out/pindel_calls_chrY -d 01022009 -G -c Y"

else:

  rule pindel_germline_1:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr1_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr1 -f {params.genome} -o pindel_out/pindel_calls_chr1 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr1.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr1 -d 01122012 -G -c chr1"

  rule pindel_germline_2:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr2_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr2 -f {params.genome} -o pindel_out/pindel_calls_chr2 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr2.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr2 -d 01122012 -G -c chr2"

  rule pindel_germline_3:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr3_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr3 -f {params.genome} -o pindel_out/pindel_calls_chr3 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr3.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr3 -d 01122012 -G -c chr3"

  rule pindel_germline_4:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr4_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr4 -f {params.genome} -o pindel_out/pindel_calls_chr4 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr4.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr4 -d 01122012 -G -c chr4"

  rule pindel_germline_5:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr5_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr5 -f {params.genome} -o pindel_out/pindel_calls_chr5 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr5.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr5 -d 01122012 -G -c chr5"

  rule pindel_germline_6:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr6_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr6 -f {params.genome} -o pindel_out/pindel_calls_chr6 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr6.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr6 -d 01122012 -G -c chr6"

  rule pindel_germline_7:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr7_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr7 -f {params.genome} -o pindel_out/pindel_calls_chr7 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr7.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr7 -d 01122012 -G -c chr7"

  rule pindel_germline_8:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr8_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr8 -f {params.genome} -o pindel_out/pindel_calls_chr8 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr8.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr8 -d 01122012 -G -c chr8"

  rule pindel_germline_9:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr9_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr9 -f {params.genome} -o pindel_out/pindel_calls_chr9 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr9.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr9 -d 01122012 -G -c chr9"

  rule pindel_germline_10:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr10_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr10 -f {params.genome} -o pindel_out/pindel_calls_chr10 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr10.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr10 -d 01122012 -G -c chr10"

  rule pindel_germline_11:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr11_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr11 -f {params.genome} -o pindel_out/pindel_calls_chr11 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr11.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr11 -d 01122012 -G -c chr11"

  rule pindel_germline_12:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr12_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr12 -f {params.genome} -o pindel_out/pindel_calls_chr12 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr12.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr12 -d 01122012 -G -c chr12"

  rule pindel_germline_13:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr13_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr13 -f {params.genome} -o pindel_out/pindel_calls_chr13 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr13.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr13 -d 01122012 -G -c chr13"

  rule pindel_germline_14:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr14_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr14 -f {params.genome} -o pindel_out/pindel_calls_chr14 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr14.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr14 -d 01122012 -G -c chr14"

  rule pindel_germline_15:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr15_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr15 -f {params.genome} -o pindel_out/pindel_calls_chr15 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr15.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr15 -d 01122012 -G -c chr15"

  rule pindel_germline_16:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr16_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr16 -f {params.genome} -o pindel_out/pindel_calls_chr16 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr16.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr16 -d 01122012 -G -c chr16"

  rule pindel_germline_17:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr17_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr17 -f {params.genome} -o pindel_out/pindel_calls_chr17 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr17.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr17 -d 01122012 -G -c chr17"

  rule pindel_germline_18:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr18_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr18 -f {params.genome} -o pindel_out/pindel_calls_chr18 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr18.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr18 -d 01122012 -G -c chr18"

  rule pindel_germline_19:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chr19_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chr19 -f {params.genome} -o pindel_out/pindel_calls_chr19 --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chr19.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chr19 -d 01122012 -G -c chr19"

  rule pindel_germline_X:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chrX_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chrX -f {params.genome} -o pindel_out/pindel_calls_chrX --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chrX.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chrX -d 01122012 -G -c chrX"

  rule pindel_germline_Y:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pindel"
    input: expand("{x}.recal.bam", x=samples),"pindel_out/pindel_config",
    output: "pindel_out/pindel_calls_chrY_INT_final",
    threads: 4
    shell: "mkdir -p pindel_out; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -c chrY -f {params.genome} -o pindel_out/pindel_calls_chrY --number_of_threads {threads}; pindel2vcf -v pindel_out/pindel_calls_chrY.vcf -r {params.genome} -R mm10 -P pindel_out/pindel_calls_chrY -d 01122012 -G -c chrY"