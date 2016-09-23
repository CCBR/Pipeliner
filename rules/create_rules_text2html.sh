#!/usr/bin/sh
# Developer: Carl McIntosh

## Usage example
# python create_rules_text2html.py \
# --rule_file test/annovar.rl \
# --rule_name annovar \
# --out_file test/rule_annovar.html

# python create_rules_text2html.py \
# -i test/annovar.rl \
# -l annovar \
# -o test/rule_annovar.html

RULES_DIR=/Users/mcintoshc/Desktop/Pipeline-test/Rules
RULES2HTML_DIR=/Users/mcintoshc/GitHub/ExomePipeline/Rules
SCRIPT_DIR=$RULES2HTML_DIR
 
## These were fullishly done by hand
#

clear;

 echo "<html>" > $RULES2HTML_DIR/MenuRules.html
 echo "<head>" >> $RULES2HTML_DIR/MenuRules.html
 echo "<title></title>" >> $RULES2HTML_DIR/MenuRules.html
 echo "</head>" >> $RULES2HTML_DIR/MenuRules.html
 echo "<img src=\"../fk_logo_thumbnail_id.jpg\" alt=\"CCBR Logo\" style=\"width:184px;height:46px\">" >> $RULES2HTML_DIR/MenuRules.html
for RULE_NAME in all-exomeseq-germline all-exomeseq-pairs all-initialqc annovar apply_variant_recal batchgvcf batchvcf callvariant checkqc combine coverage_qc fastq fastqc fastqc2 flagstats flagstats_dedup gatk_genotyping gatk_htc gatk-apply-vqsr gatk-genotyping gatk-htc gatk-variant-recal genotypegvcfs genotypegvcfs2 headers index_ref initialqc map_pe markdups mpileup ngsqc novoalign qualimap realign recal sam2bam sort trimmomatic variant_annotator variant_annotator_alt variant_annotator2 variant_recal varscan
do
  echo "Working on $RULES2HTML_DIR/rule_$RULE_NAME.html"
  python create_rules_text2html.py \
  --rule_file $RULES_DIR/$RULE_NAME.rl \
  --rule_name $RULE_NAME \
  --out_file $RULES2HTML_DIR/rule_$RULE_NAME.html


  echo "<!-- #### $RULE_NAME #### -->" >> $RULES2HTML_DIR/MenuRules.html
  echo "<br><a href=\"rule_$RULE_NAME.html\" target=\"view\">$RULE_NAME.rl</a>" >> $RULES2HTML_DIR/MenuRules.html
  echo "" >> $RULES2HTML_DIR/MenuRules.html
done

echo "<body>" >> $RULES2HTML_DIR/MenuRules.html
echo "</body>" >> $RULES2HTML_DIR/MenuRules.html
echo "<html>" >> $RULES2HTML_DIR/MenuRules.html
