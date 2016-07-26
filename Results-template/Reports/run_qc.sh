#!/bin/bash


qcdir=../QC
qualimap_pat=*quali*html
qualimap_dir=

/ngser.pl ../QC ./

python2 summarize_qualimap_reports.py *quali*html ../QC/qualimap/ ./

python2 summarize_fastqc_reports.py fastqc.zip ../QC ./
