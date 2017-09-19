#!/bin/sh

mkdir all
mkdir all/merged
mkdir all/fc
mkdir all/sf
mkdir all/of
mkdir summary

ls starfusion/oncofuse/ > mylist.txt
cat mylist.txt | while read line
do
    python Scripts/fusingFusion.py fusioncatcher/$line/final-list_candidate-fusion-genes.txt starfusion/$line/star-fusion.fusion_candidates.final starfusion/oncofuse/$line/$line.oncofuse.output fusioncatcher/oncofuse/$line/$line.oncofuse.output starfusion/fusioninspector/$line/$line.fusion_predictions.final fusioncatcher/fusioninspector/$line/$line.fusion_predictions.final all/merged/$line.txt all/fc/$line.txt all/sf/$line.txt all/of/$line.txt
done


rm mylist.txt 

cd Scripts/

python combine.py ../all/merged/
cut -d$'\t' -f 3- ../all/easySum.txt > ../summary/allSamples.tsv
python combine.py ../all/fc/
cut -d$'\t' -f 3- ../all/easySum.txt > ../summary/fcSamples.tsv
python combine.py ../all/sf/
cut -d$'\t' -f 3- ../all/easySum.txt > ../summary/sfSamples.tsv
python combine.py ../all/of/
cut -d$'\t' -f 3- ../all/easySum.txt > ../summary/ofSamples.tsv

python fusionCounter.py ../summary/ofSamples.tsv


rm ../summary/allSamples.tsv ../summary/ofSamples.tsv

rm -rf ../all/

