rule admixture_wgs:
    input: "combined.strictFilter.vcf"
    output: vcf=temp("admixture_out/samples_noINDEL_nomissing.recode.vcf"),
            mergedvcf=temp("admixture_out/samples_and_knowns.vcf"),
            ped=temp("admixture_out/samples_and_knowns_filtered.ped"),
            map=temp("admixture_out/samples_and_knowns_filtered.map"),
            recodeped="admixture_out/samples_and_knowns_filtered_recode.ped",
            admix="admixture_out/samples_and_knowns_filtered_recode.Q",
            table="admixture_out/admixture_table.tsv",
            plot="admixture_out/admixture_mqc.png"
    params: gatk=config['bin'][pfamily]['GATK'],ref=config['project']['annotation'],genome=config['references'][pfamily]['GENOME'],key=config['references'][pfamily]['ADMIXTUREKEY'],refcount=config['references'][pfamily]['ADMIXTUREREFS'],knowns=config['references'][pfamily]['KNOWNANCESTRY'],rname="admixture"
    threads: 8
    shell: """
         mkdir -p admixture_out; module load vcftools; vcftools --vcf {input} --remove-indels --max-missing 1 --recode --recode-INFO-all --out admixture_out/samples_noINDEL_nomissing; module load GATK/3.5-0; GATK -m 48G CombineVariants -R {params.genome} --genotypemergeoption UNSORTED -o {output.mergedvcf} --variant {params.knowns} --variant admixture_out/samples_noINDEL_nomissing.recode.vcf --minimumN 2 -nt 4; vcftools --vcf {output.mergedvcf} --maf 0.05 --remove-indels --plink --out admixture_out/samples_and_knowns_filtered; module load plink/1.9.0-beta4.4; plink --noweb --recode12 --out admixture_out/samples_and_knowns_filtered_recode --file admixture_out/samples_and_knowns_filtered; module load perl/5.18.4; perl Scripts/admixture_prep.pl {params.key} admixture_out/samples_and_knowns_filtered_recode.pop admixture_out/samples_and_knowns_filtered_recode.ped; /data/CCBR_Pipeliner/db/PipeDB/bin/admixture_linux-1.3.0/admixture admixture_out/samples_and_knowns_filtered_recode.ped {params.refcount} --supervised -j32; mv samples_and_knowns_filtered_recode.{params.refcount}.P admixture_out/samples_and_knowns_filtered_recode.P; mv samples_and_knowns_filtered_recode.{params.refcount}.Q admixture_out/samples_and_knowns_filtered_recode.Q; perl Scripts/admixture_post.pl {params.key} {output.table} {output.admix} {params.ref} {output.recodeped}; module load R/3.5; Rscript Scripts/admixplot.R

           """