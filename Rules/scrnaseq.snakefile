from snakemake.utils import R

configfile: "run.json"

if config['project']['pipeline'] == "cellranger":
  rule all:
     params: batch='--time=168:00:00',crid=config['project']['CRID']
     input: "{params.crid}/outs/web_summary.html"
            
        
elif config['project']['pipeline'] == "scrnaseqinit":
  rule all:
     params: batch='--time=168:00:00'
     input: "scrna_initial.html", "scrna_jackstraw.html"

elif config['project']['pipeline'] == "scrnaseqcluster":
  rule all:
     params: batch='--time=168:00:00'
     input: "scrna_cluster.html"

rule cellranger: 
   params: rname='pl:cellranger',batch='--cpus-per-task=40 --mem=110g --time=48:00:00',crid=config['project']['CRID'],refer=config['project']['annotation'],dir=config['project']['workpath'],expected=config['project']['EXPECTED'],projectId=config['project']['id']
   output: "{crid}/outs/web_summary.html".format(crid=config['project']['CRID'])
   shell: """
          module load cellranger;
          cellranger count --id={params.projectId} \
                   --transcriptome=$CELLRANGER_REF/{params.refer} \
                   --fastqs={params.dir} \
                   --sample={params.crid} \
                   --localcores=32 \
                   --localmem=100 \
                   --jobmode=slurm --maxjobs=30 > cellranger_{params.projectId}.log
          """

rule scrna_initial: 
   params: rname='pl:scrnainit',batch='--cpus-per-task=40 --mem=110g --time=48:00:00',countspath=config['project']['COUNTSPATH'],refer=config['project']['annotation'],dir=config['project']['workpath'],projDesc=config['project']['description'],projectId=config['project']['id']
   output: "scrna_initial.html", "{params.projectId}_seurat_object.rds"
   shell: "module load R/3.4.0_gcc-6.2.0; Rscript Scripts/scrna_initial_call.R '{params.dir}' '{params.refer}' '{params.projectId}' '{params.projDesc}'"

rule scrna_jackstraw: 
   params: rname='pl:scrnajackstraw',batch='--cpus-per-task=40 --mem=110g --time=48:00:00',dir=config['project']['workpath'],projDesc=config['project']['description'],projectId=config['project']['id']
   input: so="{projectId}_seurat_object.rds".format(projectId=config['project']['id'])
   output: "scrna_jackstraw.html"
   shell: "module load R/3.4.0_gcc-6.2.0; Rscript Scripts/scrna_jackstraw_call.R '{params.dir}' '{input.so}' '{params.projectId}' '{params.projDesc}'"

rule scrna_cluster: 
   params: rname='pl:scrnacluster',batch='--cpus-per-task=40 --mem=110g --time=48:00:00',dir=config['project']['workpath'],pcs=config['project']['PCS'],resolution=config['project']['RESOLUTION'],projDesc=config['project']['description'],projectId=config['project']['id']
   input: so="{projectId}_seurat_object.rds".format(projectId=config['project']['id'])
   output: "scrna_cluster.html"
   shell: "module load R/3.4.0_gcc-6.2.0; Rscript Scripts/scrna_cluster_call.R '{params.dir}' '{input.so}' '{params.pcs}' '{params.resolution}' '{params.projectId}' '{params.projDesc}'"