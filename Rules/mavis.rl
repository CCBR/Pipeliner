rule rule mavis:
     input:  star="starfusion/fusioninspector/{x}/{x}.fusion_predictions.final",fusioncatcher="fusioncatcher/fusioninspector/{x}/{x}.fusion_predictions.final",
     output: ???
     params: rname='mavis',sample="{x}",starlib=config['references'][pfamily]['STARFUSIONLIB']
     shell: "source /data/CCBR_Pipeliner/db/PipeDB/bin/conda/etc/profile.d/conda.sh
             conda activate base
             mkdir -p mavis
             mkdir -p mavis/{params.sample}
             /data/CCBR_Pipeliner/db/PipeDB/bin/mavis/bin/mavis convert --file_type starfusion --outputfile {params.sample}_starfusion --inputs starfusion/fusioninspector/{params.sample}/{params.sample}.fusion_predictions.final.abridged --strand_specific TRUE
             /data/CCBR_Pipeliner/db/PipeDB/bin/mavis/bin/mavis convert --file_type starfusion --outputfile {params.sample}_fusioncatcher --inputs fusioncatcher/fusioninspector/{params.sample}/{params.sample}.fusion_predictions.final.abridged --strand_specific TRUE
             mkdir -p {params.sample}_merged
             /data/CCBR_Pipeliner/db/PipeDB/bin/mavis/bin/mavis cluster --disease_status diseased --library {params.sample} --protocol transcriptome --strand_specific TRUE --inputs {params.sample}_starfusion {params.sample}_fusioncatcher -o {params.sample}_merged
             /data/CCBR_Pipeliner/db/PipeDB/bin/mavis/bin/mavis annotate --annotations /data/CCBR_Pipeliner/db/PipeDB/lib/mavis_resources/hg19/reference_inputs/ensembl69_hg19_annotations.json --library {params.sample}_starfusion --protocol transcriptome --reference_genome /data/CCBR_Pipeliner/db/PipeDB/lib/mavis_resources/hg19/reference_inputs/hg19.fa -n test_input -o test_annotate
             mkdir -p {params.sample}_pairings
             /data/CCBR_Pipeliner/db/PipeDB/bin/mavis/bin/mavis pairing --annotations /data/CCBR_Pipeliner/db/PipeDB/lib/mavis_resources/hg19/reference_inputs/ensembl69_hg19_annotations.json --inputs {params.sample}_annotate/annotations.tab -o {params.sample}_pairings
             mkdir -p {params.sample}_summary
             /data/CCBR_Pipeliner/db/PipeDB/bin/mavis/bin/mavis summary --annotations /data/CCBR_Pipeliner/db/PipeDB/lib/mavis_resources/hg19/reference_inputs/ensembl69_hg19_annotations.json --inputs {params.sample}_pairings/mavis_paired.tab -o {params.sample}_summary --contig_call_distance 100 --flanking_call_distance 100 --spanning_call_distance 100 --split_call_distance 100"