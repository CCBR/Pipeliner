#!/usr/bin/env python

import sys
import os
import re






date = os.popen('date').read()

V = os.popen('snakemake -v')
v = V.read()

    
F = open('stats.css', 'r')
CSS = F.read()
F.close()

F = open('toc.js', 'r')
TOC = F.read()
F.close()

F = open('sorttable.js', 'r')
SRT = F.read()
F.close()


F = open('../run.json', 'r')
P = eval(F.read())
F.close()
pid = P['project']['id']
org = P['project']['organism']
user = P['project']['analyst']
pipeline = P['project']['pipeline']
pipeline_ver = P['project']['version']


if pipeline=="initialqc":
    includeFastqcSummary=True
    includeNgsqcSummary=True
    includeQualimapSummary=True

includeFastqcSummary=True
includeNgsqcSummary=True
includeQualimapSummary=True
useStats=True


out=os.popen("cd ../ && snakemake --detailed-summary --dryrun --rerun-incomplete> Reports/"+pipeline+".summary",'r',1).read()


try:
    F = open(pipeline+'.stats', 'r')
    f=F.read()
    D=eval(f)
    F.close()
except:
    print("No "+pipeline+".stats file found. No run stats will be included.\n")
    useStats=False
F.close()

title = 'Pipeline '+pipeline+' Stats'
page="<html><head><style>"+CSS+"</style>"
page=page+"<title>"+title+"</title>"

page=page+"<script type='text/javascript' src='sorttable.js'></script>"
page=page+"<script type='text/javascript' src='toc.js'></script>"
page=page+"<script type='text/javascript' src='pipeliner.js'></script>"


page=page+"</head>"

page=page+"<body onload=generateTOC(document.getElementById('toc'));>"


nsamp = len(P['project']['units'])

page=page+'<h1>Project: <i>{}</i> with {} samples</h1>'.format(pid, nsamp)
page=page+'<h5>Date: <i>{}</i></h5>'.format(date)
page=page+'<h5>Organism: <i>{}</i></h5>'.format(org)
page=page+'<h5>Analyst: <i>{}</i></h5>'.format(user)
page=page+'<h5>Pipeline: <i>{}</i></h5>'.format(pipeline)
page=page+'<h5>Pipeline Version: <i>{}</i></h5>'.format(pipeline_ver)
page=page+'<h5>Snakemake Version: <i>{}</i></h5>'.format(v)

page=page+'<h2>Contents</h2>'
page=page+"<div id='toc'></div>"


page=page+"<div id='main'>"
page=page+'<h2>Pipeline Flow Diagram</h2>'
#page=page+"<img src='"+pipeline+".gif'>"

page=page+"<object type='image/svg+xml' data="+pipeline+".svg>Your browser does not support SVG</object>"

#page=page+"<iframe src="+pipeline+".svg>Your browser does not support SVG</iframe>"

if useStats == True:
    page=page+'<h3>Total Computing Time for Pipeline: {}'.format(D['total_runtime'])
    page=page+'<h2>Execution Times for Rules</h2>'
    page=page+"<table cols=2 class='sortable'>"
    page=page+'<th>Rule</th><th>Average Duration</th><th>Maximum Duration</th><th>Minimum Duration</th>'
    for k in D['rules'].keys():
        page=page+'<ul><tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>'.format(k,D['rules'][k]['mean-runtime'], D['rules'][k]['max-runtime'], D['rules'][k]['min-runtime'])
    page=page+"</table>"
    page=page+'<h2>Starts, Stops, and Processing Times for Files'
    page=page+"<table cols=4 class='sortable'>"
    page=page+'<tr>'
    page=page+'<th>File</th><th>Start</th><th>Stop</th><th>Duration</th>'
    page=page+'</tr>'

    for k in D['files'].keys():
        page=page+'<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>'.format(k,D['files'][k]['start-time'], D['files'][k]['stop-time'],D['files'][k]['duration'])

    page=page+"</table>"



page=page+'Contents of <i>project.json</i>'
page=page+"<table cols=2 class='sortable'>"

page=page+'<tr>'

page=page+'<th>Parameter</th><th>Value</th>'

page=page+'</tr>'

for k in P.keys():
    T = '<table>'
    try:
        for k2 in P[k].keys():
            T2 = '<table>'
            try:
                for k3 in P[k][k2].keys():
                    T2 = T2 \
                        + '<tr><td>{}</td><td>{}</td></tr>'.format(k3,
                            P[k][k2][k3])
                T2 = T2 + '</table>'
            except:
                T2 = P[k][k2]

            T = T + '<tr><td>{}</td><td>{}</td></tr>'.format(k2, T2)
        T = T + '</table>'
    except:
        T = P[k]

    page=page+'<tr><td>{}</td><td>{}</td></tr>'.format(k, T)

page=page+"</table>"


page=page+'<h2>FastQC Summary For All Samples</h2>'

if includeFastqcSummary==True:

    page=page+"<object type='text/html' name='fastqcsum' data='aggregate_fastqc_report.html' width='800' height='800' scrolling='yes'></object>"

page=page+'<h2>NGSQC Summary For All Samples</h2>'

if includeNgsqcSummary==True:

    page=page+"<object type='text/html' name='ngsqcsum' data='ngser.html' width='800' height='800' scrolling='yes'></object>"


if includeQualimapSummary==True:

    page=page+"<object type='text/html' name='qualimapsum data='QUALIMAP_SUMMARY_REPORT.html' width='800' height='800' scrolling='yes'></object>"


page=page+"<form name='fastqcform'>"

page=page+'<h2>FastQC Reports</h2>'


page=page+"<h5> Select Sample ID</h5><select name='fastqcselect' onchange='fastqc()'>"

first=sorted(P['project']['units'].keys())[0]
             
for unit in sorted(P['project']['units'].keys()):
    page=page+"<option value='%s'>%s</option>" % (unit,unit)
page=page+"</select>"


page=page+"<object type='text/html' name='fastqcr1html' data='../QC/%s.R1_fastqc.html' width='800' height='800' scrolling='yes'></object>" % first

page=page+"<object type='text/html' name='fastqcr2html' data='../QC/%s.R2_fastqc.html' width='800' height='800' scrolling='yes'></object>" % first


page=page+'<h2>FastQC Reports for Trimmed Reads</h2>'


page=page+"<h5> Select Sample ID</h5><select name='fastqcselecttrimmed' onchange='fastqctrimmed()'>"


for unit in sorted(P['project']['units'].keys()):
    page=page+"<option value='%s'>%s</option>" % (unit,unit)
page=page+"</select>"


page=page+"<object type='text/html' name='fastqcr1trimmedhtml' data='../QC/%s.R1.trimmed_fastqc.html' width='800' height='800' scrolling='yes'></object>" % first

page=page+"<object type='text/html' name='fastqcr2trimmedhtml' data='../QC/%s.R2.trimmed_fastqc.html' width='800' height='800' scrolling='yes'></object>" % first





page=page+'<h2>QualiMap Reports</h2>'


page=page+"<h5> Select Sample ID</h5><select name='qualimapselect' onchange='qualimap()'>"
for unit in sorted(P['project']['units'].keys()):
    page=page+"<option value='%s'>%s</option>" % (unit,unit)
page=page+"</select>"


page=page+"<object type='text/html' name='qualimaphtml' data='../QC/%s.qualimapReport/qualimapReport.html' width='800' height='400' scrolling='yes'></object>" % first


page=page+"<h2>Snakemake Files Summary</h2>"


page=page+"<object type='text/html' name='summaryhtml' data='../Reports/%s.summary' width='800' height='400' scrolling='yes'></object>" % pipeline

page=page+"</form>"




page=page+"</body></html>"

G = open(pipeline+'.html', 'w')

G.write(str(page))
G.close()
