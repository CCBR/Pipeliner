#!/usr/bin/env python
import sys,os,re
import json
from subprocess import Popen, PIPE, STDOUT

whereiam=os.popen("pwd").read().strip()

C=open("project.json","r")
C=eval(C.read())
pipeline=C['project']['pipeline']
pfamily=C['project']['pfamily']
WP=C['project']['workpath']

# jobscripts for PBS and slurm

#F={"qsub":"/pipeline_ctrl.sh","submit":"/submit.sh","slurm":"/pipeline_ctrls.sh","submit_slurm":"/submit_slurm.sh"}

F={"slurm":"/pipeline_ctrl.sh","submit_slurm":"/submit_slurm.sh"}

parameters=" ".join(C['project']['smparams'])

for S in F.keys():

    Q=open(S+".template","r")
    T=Q.read()
    Q.close()
    T=re.sub("R=/.+","R="+WP,T)
    T=re.sub("D=/.+","D="+whereiam,T)
    T=re.sub("Stats",pipeline+".stats",T)
    T=re.sub("snakemake ","snakemake %s "%parameters,T)    
    Q=open(WP+F[S],"w")
    Q.write(T)
    Q.close()
    p = os.popen("chmod ug+x "+WP+F[S])

annotation=C['project']['annotation']
binset=C['project']['binset']

D=open("rules.json","r")
E=open(annotation+".json","r")
F=open(binset+".json","r")
#H=open("rnaseq.json","r")



D=eval(D.read())
E=eval(E.read())
F=eval(F.read())
#H=eval(H.read())

if pipeline=="custom":
    for i in C['project']['custom']:
        D['rules'][i].append("custom")
    if C['project']['bysample']=="yes":
        Q=open(whereiam+"/Rules/all-custom.template","r")
    else:
        Q=open(whereiam+"/Rules/all-custom.combined.template","r")    
    allrule=Q.read()
    Q.close()
    allrule=re.sub("ext='bam'","ext='%s'"%C['project']['efiletype'],allrule)
    if C['project']['efiletype']=="R1_fastqc.html" or C['project']['efiletype']=="qualimapReport":
        prefix="QC/"
        allrule=re.sub("prefix=''","prefix='%s'"%prefix,allrule)
        
    Q=open(whereiam+"/Rules/all-custom.rl","w")
    Q.write(allrule)
    Q.close()

config=dict(list(C.items())+list(D.items())+list(E.items())+list(F.items()))

PipelineType=config['project']['pipeline']

cluster= config['project']['cluster']
p = os.popen("cp {0}/cluster_{1}.json {2}/cluster.json".format(whereiam,cluster,WP))


with open(C['project']['workpath']+'/run.json', 'w') as F:
#    json.dump(config, F)
    json.dump(config, F, sort_keys = True, indent = 4,ensure_ascii=False)

F.close()

G=open("Rules/"+config['final'][PipelineType],"r")
All=G.read()
T=open("snakefile.template","r")
Template=T.read()

BS=C['project']['batchsize']
workpath=C['project']['workpath']
Template=re.sub("BS=50","BS="+BS,Template)

F=open(workpath+"/Snakefile","w")
F.write(Template+"\n")
F.write(All+"\n")

for R in sorted(config['rules'].keys()):
    if PipelineType in config['rules'][R]:
        F.write("include: \""+whereiam+"/Rules/"+R+".rl\"\n")

F.close()

os.popen("cd "+workpath+"/ && snakemake  --rulegraph|dot -Tpng -o "+workpath+"/Reports/"+PipelineType+".png")
os.popen("cd "+workpath+"/ && convert "+workpath+"/Reports/"+PipelineType+".png "+workpath+"/Reports/"+PipelineType+".gif")

#Make clickable SVG

os.popen("cd "+workpath+"/ && snakemake  --rulegraph > "+workpath+"/Reports/"+PipelineType+".dot")

F=open(workpath+"/Reports/"+PipelineType+".dot","r")
f=F.read()
F.close()

f2=re.sub(r"label = \"(.+?)\",", r"label = '\1' URL='../Rules/\1.rl'", f)



F=open(workpath+"/Reports/"+PipelineType+".dot2","w")
F.write(f2)
F.close()

F=open(workpath+"/Reports/"+PipelineType+".dot2","r")
f=F.read()
f2=re.sub(r'\'', r'"', f)
f2=re.sub(r'_', r'.', f2)
f2=re.sub(r'.dag', r'_dag', f2)

F.close()

#print(f2)

F=open(workpath+"/Reports/"+PipelineType+".dot2","w")
F.write(f2)
F.close()

os.popen("cd "+workpath+"/Reports && dot -Tsvg "+PipelineType+".dot2 > "+PipelineType+".svg")

os.popen("cd "+workpath+"/ && snakemake --detailed-summary > Reports/"+PipelineType+".summary")

