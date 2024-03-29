Welcome to **Pipeliner** - an open-source and scalable solution to NGS analysis powered by the NIH's [Biowulf](https://hpc.nih.gov/) cluster. 

Pipeliner provides access to a set of best-practices NGS pipelines developed, tested, and benchmarked by experts at [CCBR](https://ccbr.ccr.cancer.gov/) and [NCBR](https://ncbr.ncifcrf.gov/). This wiki is the main source of documentation for users and developers working with or contributing to [`Pipeliner`](https://github.com/CCBR/Pipeliner). If you are unfamiliar with Pipeliner or this is the first time hearing about this project, please check out the links below.

## Quick Navigation

| <img src="https://github.com/CCBR/Pipeliner/wiki/Images/assets/dna_icon.png" width="64px" height="64px"/> |<img src="https://github.com/CCBR/Pipeliner/wiki/Images/assets/group_icon.png" width="64px" height="64px"/> | <img src="https://github.com/CCBR/Pipeliner/wiki/Images/assets/question_icon.png" width="64px" height="64px"/> | <img src="https://github.com/CCBR/Pipeliner/wiki/Images/assets/doc_icon.png" width="64px" height="64px"/> |
|:----------------------------:|:---------------------------------:|:-------------------------------:|:---------------------------:|
| **[Getting Started](https://github.com/CCBR/Pipeliner/wiki/1.-Introduction-to-Pipeliner#introduction-to-pipeliner)** | **[Project & Community](https://github.com/CCBR/Pipeliner/wiki/Community-and-contribution)**       | **[FAQs](https://github.com/CCBR/Pipeliner/wiki/7.-FAQ#faq)** | **[Pipeline Documentation](https://github.com/CCBR/Pipeliner/wiki/Pipeline-Documentation)** |
| Introducing Pipeliner, prerequisites and how to get started | About our community, how to contribute and make feature requests | Frequently asked questions and solutions to common problems | Detailed technical documentation for users and developers about each of our pipelines |

## Questions or need help?
Please check out our [FAQ](https://github.com/CCBR/Pipeliner/wiki/7.-FAQ#faq) or [contact](https://github.com/CCBR/Pipeliner/wiki/Contact-us) page for different ways of getting in touch with the team.





## Getting started
1. [Pre-requisites](https://github.com/CCBR/Pipeliner/wiki/1.-Introduction-to-Pipeliner#pre-requisites)  
2. [How to launch the Pipeliner GUI](https://github.com/CCBR/Pipeliner/wiki/1.-Introduction-to-Pipeliner#how-to-launch-the-pipeliner-gui)


## RNA-seq  
#### I. [Gene and isoform expression pipeline](https://github.com/CCBR/Pipeliner/wiki/Gene-and-isoform-expression-pipeline)
1. [Introduction](https://github.com/CCBR/Pipeliner/wiki/Gene-and-isoform-expression-pipeline#1-introduction)    
    1.1 [Theory and practical guide to RNA-seq](https://github.com/CCBR/Pipeliner/wiki/Theory-and-practical-guide-for-RNA-seq)    
    1.2 [Supported Reference Genomes](https://github.com/CCBR/Pipeliner/wiki/Differential-expression-pipeline-tools-and-versions#reference-genomes)  
    1.3 [Tools and Versions](https://github.com/CCBR/Pipeliner/wiki/Differential-expression-pipeline-tools-and-versions#tools-and-versions)  
    &emsp;1.3.1 [Quantification and QC pipeline](https://github.com/CCBR/Pipeliner/wiki/Differential-expression-pipeline-tools-and-versions#quality-control-pipeline)  
    &emsp;1.3.2 [Differential Expression pipeline](https://github.com/CCBR/Pipeliner/wiki/Differential-expression-pipeline-tools-and-versions#differential-expression-pipeline)   
2. [Overview](https://github.com/CCBR/Pipeliner/wiki/Gene-and-isoform-expression-pipeline#2-overview)  
    2.1 [Quantification and QC pipeline](https://github.com/CCBR/Pipeliner/wiki/Gene-and-isoform-expression-pipeline#21-quantification-and-quality-control-pipeline)  
    2.2 [Differential Expression pipeline](https://github.com/CCBR/Pipeliner/wiki/Gene-and-isoform-expression-pipeline#22-differential-expression-pipeline)  
3. [Tutorial](https://github.com/CCBR/Pipeliner/wiki/Gene-and-isoform-expression-pipeline#3-tutorial)   
    3.1 [Quantification and QC pipeline](https://github.com/CCBR/Pipeliner/wiki/Gene-and-isoform-expression-pipeline#31-getting-started-with-the-quantification-and-quality-control-pipeline)  
4. [TLDR](https://github.com/CCBR/Pipeliner/wiki/TLDR-RNA-seq)  
    4.1 [Launch Pipeliner](https://github.com/CCBR/Pipeliner/wiki/TLDR-RNA-seq#launch-pipeliner)  
    4.1 [Run the Quantification and QC pipeline](https://github.com/CCBR/Pipeliner/wiki/TLDR-RNA-seq#run-quantification-and-quality-control-pipeline)  
5. [References](https://github.com/CCBR/Pipeliner/wiki/Differential-expression-pipeline-tools-and-versions#references) 

#### II. [Gene-fusions pipeline](https://github.com/CCBR/Pipeliner/wiki/4.-RNASeq-Fusion-Detection)
1. [Getting Started](https://github.com/CCBR/Pipeliner/wiki/4.-RNASeq-Fusion-Detection#getting-started)
2. [Output from Quality-control](https://github.com/CCBR/Pipeliner/wiki/4.-RNASeq-Fusion-Detection#quality-control)
3. [Output from Gene-fusions analysis](https://github.com/CCBR/Pipeliner/wiki/4.-RNASeq-Fusion-Detection#gene-fusion-outputs)
4. [Additional Outputs](https://github.com/CCBR/Pipeliner/wiki/4.-RNASeq-Fusion-Detection#additional-outputs)


#### III. [RNA variant calling pipeline](https://github.com/CCBR/Pipeliner/wiki/6.-RNAseq-Variant-Calling)
1. [Getting Started](https://github.com/CCBR/Pipeliner/wiki/6.-RNAseq-Variant-Calling#getting-started)
2. [Outputs from Quality-control](https://github.com/CCBR/Pipeliner/wiki/6.-RNAseq-Variant-Calling#quality-control)
3. [Outputs from Variant Calling](https://github.com/CCBR/Pipeliner/wiki/6.-RNAseq-Variant-Calling#variant-calling-outputs)
4. [Additional Outputs](https://github.com/CCBR/Pipeliner/wiki/6.-RNAseq-Variant-Calling#additional-outputs)


## ChIP-seq
#### I. [Quality-control pipeline](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq)
1. [Getting started with the Quality-control pipeline](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#getting-started-with-the-quality-control-pipeline)  
2. [About the Demo Dataset](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#about-the-demo-dataset)  
3. [Tutorial](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#phase-1-initial-qc)  
  3.1 [Step 0. Fill out the `Project Information` section](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#step-0-fill-out-the-project-information-section)  
  3.2 [Step 1. Setting `Pipeline`](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#step-1-setting-pipeline)  
  3.3 [Step 2. Select `Genome`](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#step-2-select-genome)  
  3.4 [Step 3. Select your Data Directory](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#step-3-select-your-data-directory)  
  3.5 [Step 4. Select your Working Directory](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#step-4-select-your-working-directory)  
  3.6 [Step 5. Initialize your Working Directory](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#step-5-initialize-your-working-directory)  
  3.7 [Step 6. Set Peak information](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#step-6-set-peak-information)  
  3.8 [Step 7. Perform Dry Run](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#step-7-perform-dry-run)  
  3.9 [Step 8. Run](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#step-8-run)  
4. [Check progress](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#check-progress)  
5. [Confirm successful completion of Phase 1](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#confirm-successful-completion-of-phase-1)  
6. [Output files and directory structure](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq-phase1-output-files)
7. [Tools and Versions](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq-tools)  
  7.1 [Quality control assessment tools](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq-tools#quality-control-assessment-tools)  
  7.2 [Data processing tools](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq-tools#data-processing-tools)  
 

#### II. [Peak calling, differential peak calling, annotations, and motif pipeline](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#phase-2-peak-calling-differential-peak-calling-annotations-and-motif-searches)  
1. [Phase 2: Peak calling, differential peak calling, annotations, and motif searches](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#phase-2-peak-calling-differential-peak-calling-annotations-and-motif-searches)  
  1.1 [Step 1. Initial Setup](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#step-1-initial-setup)  
  1.2 [Step 2. Select Pipeline Options](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#step-2-select-pipeline-options)  
  1.3 [Step 3. Set contrasts (Optional)](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#step-3-set-contrasts-optional)  
  1.4 [Step 4. Perform Dry Run](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#step-4-perform-dry-run)  
  1.5 [Step 5. Run](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#step-5-run)  
2. [Check progress](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq#check-progress-1)  
3. [Output files and directory structure](https://github.com/CCBR/Pipeliner/wiki/ChIP-seq-phase2-output-files)  

## DNA-seq
#### [I. Whole exome sequencing pipeline](https://github.com/CCBR/Pipeliner/wiki/3.-Whole-Exome-Sequencing)
1. [Getting Started](https://github.com/CCBR/Pipeliner/wiki/3.-Whole-Exome-Sequencing#getting-started)
2. [Quality-control](https://github.com/CCBR/Pipeliner/wiki/3.-Whole-Exome-Sequencing#initialqc-workflow-and-outputs)
3. [Variant Detection Workflows](https://github.com/CCBR/Pipeliner/wiki/3.-Whole-Exome-Sequencing#wes-variant-detection-workflows)  
    3.1 [Germline variant calling pipeline](https://github.com/CCBR/Pipeliner/wiki/3.-Whole-Exome-Sequencing#germline-variant-detection-pipeline)  
    3.2 [Tumor-Normal somatic variant calling pipeline](https://github.com/CCBR/Pipeliner/wiki/3.-Whole-Exome-Sequencing#tumor-normal-somatic-variant-detection-pipeline)    
    3.3 [Tumor-only somatic variant calling pipeline](https://github.com/CCBR/Pipeliner/wiki/3.-Whole-Exome-Sequencing#tumor-only-somatic-variant-detection-pipeline)  

#### [II. Whole genome sequencing pipeline](https://github.com/CCBR/Pipeliner/wiki/5.-Whole-Genome-Sequencing)
1. [Getting Started](https://github.com/CCBR/Pipeliner/wiki/5.-Whole-Genome-Sequencing#getting-started)
2. [Quality-control](https://github.com/CCBR/Pipeliner/wiki/5.-Whole-Genome-Sequencing#initialqc-workflow-and-outputs)
3. [Variant Detection Workflows](https://github.com/CCBR/Pipeliner/wiki/5.-Whole-Genome-Sequencing#wgs-variant-detection-workflows)   
    3.1 [Germline variant calling pipeline](https://github.com/CCBR/Pipeliner/wiki/5.-Whole-Genome-Sequencing#germline-variant-detection-pipeline)  
    3.2 [Tumor-Normal somatic variant calling pipeline](https://github.com/CCBR/Pipeliner/wiki/5.-Whole-Genome-Sequencing#tumor-normal-somatic-variant-detection-pipeline)  
    3.3 [Tumor-only somatic variant calling pipeline](https://github.com/CCBR/Pipeliner/wiki/5.-Whole-Genome-Sequencing#tumor-only-somatic-variant-detection-pipeline)  
