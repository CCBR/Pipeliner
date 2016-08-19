#!/usr/bin/env python3
import json

config=eval(open("../run.json","r").read())

samples=sorted(list(config['project']['units'].keys()))
pairs=sorted(list(config['project']['pairs'].keys()))
pfamily=config['project']['pfamily']
pipeline=config['project']['pipeline']
if pfamily=="rnaseq":
    qcfolder="QC"

D=dict()

D['ProjectSummary']=dict()
L=list()
L.append("../run.json")
D['ProjectSummary']['include']=L
D['QCWorkflow']=dict()
L=list()
L.append("../Reports/initialqcrnaseq.png")
D['QCWorkflow']['include']=L
D['ExomeWorkflow']=dict()
L=list()
L.append("../Reports/rnaseq.png")
D['ExomeWorkflow']['include']=L
D['MultiQC']=dict()
L=list()
L.append("../Reports/multiqc_report.html")
D['MultiQC']['include']=L


L=list()
for s in samples:
    L.append("../rawQC/"+s+".R1_fastqc.html")
    L.append("../rawQC/"+s+".R2_fastqc.html")    

D['RawFastQC']=dict()
D['RawFastQC']['include']=L

with open(config['project']['workpath']+'/Reports/report.json', 'w') as F:
    json.dump(D, F, sort_keys = True, indent = 4,ensure_ascii=False)

F.close()


