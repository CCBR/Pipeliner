ext='R1_fastqc.html'
prefix='QC/'
rule all_custom:
    input:  expand("{p}{s}.{e}",e=ext,s=samples,p=prefix)
    output: 
