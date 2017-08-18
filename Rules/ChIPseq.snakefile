from snakemake.utils import R
from os.path import join, abspath, basename
from os import environ as env
from io import StringIO
from tempfile import TemporaryFile

from pysam import Samfile, FastaFile
from collections import Counter

#pipehome = '/scratch/kimb8/Pipeliner/'
pipehome = '/home/kopardevn/Pipeliner/'
bam_dir='bam'

def normalize_bam_file_chromosomes(
    bamfns,
    obamfns=[],
    suffix='.common_chrom.bam' ) :
    
    counts = []
    for bamfn1 in bamfns :
        bam1 = Samfile(bamfn1)

        cnt1 = Counter()
        for aread1 in bam1 :
            if aread1.is_unmapped :
                continue
            cnt1[ aread1.reference_name ] += 1 
        bam1.close()
        counts.append(cnt1)
    
    common = None
    for cnt in counts :
        if common != None :
            common = common.intersection( cnt )
        else :
            common = set( cnt )
    
    #################
    # check if common chromosomes exists!
    assert common
    
    for i, bamfn1 in enumerate(bamfns) :
        bam1 = Samfile(bamfn1)
        
        if obamfns :
            obamfn1 = obamfns[i] #list of output BAM filenames
        else :
            obamfn1 = bamfn1.replace( '.bam', '.common_chrom.bam' )
        obam1 = Samfile( obamfn1, 'wb', template=bam1 )
        
        for aread1 in bam1:
            if aread1.is_unmapped :
                continue
                
            if aread1.reference_name in common and len(aread1.reference_name) < 6 :
                obam1.write(aread1)
                
        obam1.close()
        bam1.close()
    

configfile: "run.json"
include: join( pipehome, "Rules", "InitialChIPseqQC.snakefile" )
#include: join( "Rules", "InitialChIPseqQC.snakefile" )
    
workpath = config['project']['workpath']    
filetype = config['project']['filetype']
readtype = config['project']['readtype']

bam_suffix = '.sorted.bam'
#####################
#samples to process
#####################
samples = config['project']['peaks']['chips']
chip2input = config['project']['peaks']['inputs']

inputs = [v for v in chip2input.values() if v] #need to be removed!
uniq_inputs = list(sorted(set([v for v in chip2input.values() if v])))

groups = config['project']['groups'] #group2chips
chip2group = {}
for group, chips in groups.items() :
    for chip in chips :
        chip2group[chip]=group

#ordered group!!
ogroups = [g for g in [chip2group[c] for c in samples]]

contrast = config['project']['contrast']
clabel_sep = '_vs_'
clabels = [g1 + clabel_sep + g2 for g1,g2 in contrast]

#debugging
print( 'groups:',groups, file=sys.stderr)
print( 'samples:', samples, file=sys.stderr)
print( 'ogroups:', ogroups, file=sys.stderr)

#####################
#Peak caller Output info
#####################
macsn_dir = 'macs_narrow'
macsb_dir = 'macs_broad'
sicer_dir = "sicer"
peak_dirs = [macsn_dir, macsb_dir, sicer_dir]

macsn_suffix = '_peaks.narrowPeak'
macsb_suffix = '_peaks.broadPeak'
sicer_suffix = '_broadpeaks.bed'
peak_suffixes = [macsn_suffix, macsb_suffix, sicer_suffix]

callers = ['macsn', 'macsb', 'sicer']
caller2dir = {k:v for k,v in zip(callers, peak_dirs)}
dir2caller = {v:k for k,v in zip(callers, peak_dirs)}
suffix2caller = {v:k for k,v in zip(callers, peak_suffixes)}
caller2suffix = {k:v for k,v in zip(callers, peak_suffixes)}

genome = config['project']['annotation']
if genome == 'hg19' :
    ngsplot_regions = ["tss", "tes", "genebody", 'cgi', "dhs", "enhancer", "exon"]
else :
    ngsplot_regions = ["tss", "tes", "genebody", 'cgi', "exon"]


#####################
# directories
#####################
#sicer_rmdup_dir = "sicer_rmdup"
memechip_dir = "MEMEChIP"
idr_dir = 'IDR'
pepr_dir = 'PePr'
ngsplot_dir = "Reports"
ceas_dir = 'CEAS'
chipseeker_dir = 'ChIPseeker'

#now it is deprecated
#because I put the chr from BWA process
#in InitialChIPseqQC pipeline.
def fix_chrom( fn ) :
    changed = 0
    fp = TemporaryFile('w+')
    for l in open( fn ) :
        line =[e.strip() for e in l.split('\t')]
        if not line[0].startswith( 'chr' ) :
            changed = 1
            if line[0] == 'MT' :
                line[0] = 'chrM'
            else :
                line[0] = 'chr'+line[0]
        print('\t'.join(line), file=fp, sep='\t')
        
    fp.seek(0)
    if changed :
        ofp = open(fn, 'w')
        for l in fp :
            ofp.write(l)
        ofp.close()
            

#####################
# Local programs
#####################
#chipseeker_rmd = join( pipehome, 'ChIPseeker.rmd' )
#chipseeker_r = join( pipehome, 'runChIPseeker.R' )
chipseeker_rmd = join( 'Scripts' , 'ChIPseeker.rmd' )
chipseeker_r = join( 'Scripts' , 'runChIPseeker.R' )
    
if genome == 'hg19' :
    rule ChIPseq:
        params: 
            batch='--time=168:00:00'
        input: 
            expand(
                join('ChIPQC','{caller}','ChIPQCreport.html'),
                caller=callers
            ),
            [ join(macsn_dir, "{g}", "{n}"+macsn_suffix ).format( n=n, g=g )
             for n,g in zip(samples, ogroups) ],
            #expand(
            #    join(macsn_dir,"{group}","{name}"+macsn_suffix),
            #    zip, name=samples,
            #         group=ogroups
            #),
            [join(macsb_dir,"{g}","{n}"+macsb_suffix).format(n=n,g=g) for n,g in zip(samples, ogroups) ],

            [join(sicer_dir,"{g}","{n}"+sicer_suffix).format(n=n,g=g) for n,g in zip(samples, ogroups) ],
            #expand( 
            #    join(sicer_dir,"{group}","{name}"+sicer_suffix), 
            #    zip, name=samples, 
            #         group=ogroups
            #),

            [join(memechip_dir,'{g}','{n}_memechip').format(n=n,g=g) for n,g in zip(samples, ogroups) ],
            #expand( 
            #    join(memechip_dir,"{group}","{name}_memechip"), 
            #    zip, name=samples, 
            #         group=ogroups
            #),
            #expand( join("{dir}","{group}.idrValue.txt"), 
            #       group = groups, 
            #       dir = [idr_mn_dir, idr_mb_dir, idr_sc_dir]),
            expand( join(pepr_dir,'{clabel}__PePr_chip1_peaks.bed'),
                   clabel=clabels ),
            expand( join(pepr_dir,'{clabel}__PePr_chip2_peaks.bed'),
                   clabel=clabels ),

            expand( join(ngsplot_dir, "{region}.heatmap.pdf"),
                   region=ngsplot_regions ),
            expand( join(idr_dir,'{caller}','{group}.idrValue.txt'), 
                   caller=callers, group=groups ),
            expand( expand( join(ceas_dir, '{caller}', '{{group}}', '{{name}}.pdf'), 
                   caller=callers), zip, group=ogroups, name=samples ),
            expand(expand(join(chipseeker_dir, '{caller}', '{{group}}', '{{name}}.html'),
                caller=callers), zip, group=ogroups, name=samples),
            #"Reports/multiqc_report.html",
else :
    rule ChIPseq:
        params: 
            batch='--time=168:00:00'
        input: 
            expand(
                join('ChIPQC','{caller}','ChIPQCreport.html'),
                caller=callers
            ),
            [ join(macsn_dir, "{g}", "{n}"+macsn_suffix ).format( n=n, g=g )
             for n,g in zip(samples, ogroups) ],
            #expand(
            #    join(macsn_dir,"{group}","{name}"+macsn_suffix),
            #    zip, name=samples,
            #         group=ogroups
            #),
            [join(macsb_dir,"{g}","{n}"+macsb_suffix).format(n=n,g=g) for n,g in zip(samples, ogroups) ],

            [join(sicer_dir,"{g}","{n}"+sicer_suffix).format(n=n,g=g) for n,g in zip(samples, ogroups) ],
            #expand( 
            #    join(sicer_dir,"{group}","{name}"+sicer_suffix), 
            #    zip, name=samples, 
            #         group=ogroups
            #),

            [join(memechip_dir,'{g}','{n}_memechip').format(n=n,g=g) for n,g in zip(samples, ogroups) ],
            #expand( 
            #    join(memechip_dir,"{group}","{name}_memechip"), 
            #    zip, name=samples, 
            #         group=ogroups
            #),
            #expand( join("{dir}","{group}.idrValue.txt"), 
            #       group = groups, 
            #       dir = [idr_mn_dir, idr_mb_dir, idr_sc_dir]),
            expand( join(pepr_dir,'{clabel}__PePr_chip1_peaks.bed'),
                   clabel=clabels ),
            expand( join(pepr_dir,'{clabel}__PePr_chip2_peaks.bed'),
                   clabel=clabels ),

            expand( join(ngsplot_dir, "{region}.heatmap.pdf"),
                   region=ngsplot_regions ),
            expand( join(idr_dir,'{caller}','{group}.idrValue.txt'), 
                   caller=callers, group=groups ),
            #expand( expand( join(ceas_dir, '{caller}', '{{group}}', '{{name}}.pdf'), 
            #       caller=callers), zip, group=ogroups, name=samples ),
            expand(expand(join(chipseeker_dir, '{caller}', '{{group}}', '{{name}}.html'),
                caller=callers), zip, group=ogroups, name=samples),

rule ChIPseeker :
    params:
        rname="pl:ChIPseeker:{group}:{name}",
        genome = config['project']['annotation']
    input:
        join( macsn_dir, '{group}', '{name}'+macsn_suffix ),
        join( macsb_dir, '{group}', '{name}'+macsb_suffix ),
        join( sicer_dir, '{group}', '{name}'+sicer_suffix ),
        
    output:
        expand(
            join(chipseeker_dir, '{caller}', '{{group}}', '{{name}}.html'),
            caller=callers
        )
    shadow: 'shallow'
        
    run:
        for in_fn, out_fn in zip(input, output) :  
            in_fn = abspath(in_fn)
            out_fn = abspath(out_fn)
            
            #fix_chrom(in_fn)
            
            shell( '''
            module load R;
            chipseeker=`basename {chipseeker_rmd}`
            #ln -s {chipseeker_rmd}
            Rscript {chipseeker_r} $chipseeker {in_fn} {params.genome} {out_fn}
            ''' )
            
rule CEAS :
    params:
        rname="pl:CEAS:{group}:{name}",
        outdirs = lambda w: [join(ceas_dir, caller, w.group)
                             for caller in callers ],
        outnames = lambda w: [join(ceas_dir, caller, w.group, w.name) 
                              for caller in callers ],
        genome = config['project']['annotation'],
    input:
        join( macsn_dir, '{group}', '{name}'+macsn_suffix ),
        join( macsb_dir, '{group}', '{name}'+macsb_suffix ),
        join( sicer_dir, '{group}', '{name}'+sicer_suffix ),
    output:
        expand(
            join(ceas_dir, '{caller}', '{{group}}', '{{name}}.pdf'),
            caller=callers
        )
        
    run:
        #print( 'outdirs:', params.outdirs )
        #print( 'outnames:', params.outnames )
        
        for in_fn, odir, oname in zip(input, params.outdirs, params.outnames) :    
            shell( '''
            mkdir -p {odir};
            ''')
            
            #take the first 5 column and make 
            #cleaned up file in the ceas directory
            in_fn2 = join(odir,basename(in_fn)) #cleaned input
            in_fp = open(in_fn2, 'w')
            for l in open(in_fn) :
                line = l.split()
                if len(line[0]) < 6 :
                    print( "\t".join(line[:5]), file=in_fp )
            in_fp.close()
            
            shell('''
            module load ceas;
            module load R;
            ceas -g /fdb/CEAS/{params.genome}.refGene \
                 -b {in_fn2} --name {oname};
            ''' )
        
rule IDR:
    params:
        rname="pl:IDR:{caller}:{group}",
        outdir=join(idr_dir,'{caller}'),
    
    input:
        lambda w: [
            join(caller2dir[w.caller], w.group, c + caller2suffix[w.caller]) 
            for c in groups[w.group]
        ]
        
    output:
        join( idr_dir, '{caller}', '{group}.idrValue.txt' )
    run:
        if wildcards.caller == 'macsn' :
            shell(
            '''
            module load idr
            mkdir -p {params.outdir}
            idr -s {input} -o {output} --input-file-type narrowPeak --plot 
            ''')
        elif wildcards.caller == 'macsb' :
            shell(
            '''
            module load idr
            mkdir -p {params.outdir}
            idr -s {input} -o {output} --input-file-type broadPeak --plot 
            ''')
        elif wildcards.caller == 'sicer' :
            input2 = [fn+'.2' for fn in input]
            shell(
            '''
            module load idr
            
            mkdir -p {params.outdir}
            
            for fn in {input}
            do
                cat $fn | awk 'BEGIN{{OFS="\\t"}} {{print $0, -log($NF)}}' >  $fn.2 ;
            done
            
            idr -s {input2} -o {output} --input-file-type bed --rank 7 --plot 
            ''')
        
rule ngsplot :
    params:
        rname="pl:NGSplot:{region}",
        genome = config['project']['annotation'],
        outbase = join( '{ngsplot_dir}', '{region}' ),
        config = '{region}.config',
        batch='--cpus-per-task=16 --mem=32g --time=24:00:00',
        name = [*samples, *inputs],
    input:
        expand( join(bam_dir,"{name}.sorted.dedup.bam"), name=[*samples,*inputs] )
    output:
        join("{ngsplot_dir}", "{region}.heatmap.pdf")
    run:
        sfp = StringIO()
        for n, fn in zip( params.name, input ) :
            print( fn, -1, n, file=sfp, sep='\t' )
            
        ofp = open(params.config, 'w' )
        ofp.write(sfp.getvalue())
        ofp.close()
        
        shell( 
        """
        module load ngsplot ;
        ngs.plot.r -G {params.genome} -R {wildcards.region} -C {params.config} -O {params.outbase}
        """)
        
rule PePr :
    params:
        rname="pl:PePr:{clabel}",
    input:
        chips1 = lambda w: [join(bam_dir,chip+".sorted.bam")
                   for chip in
                   groups[w.clabel.split(clabel_sep)[0]] ],
        
        inputs1 = lambda w: [join(bam_dir,ctrl+".sorted.bam")
                   for ctrl in
                   [ chip2input[c] for c in groups[w.clabel.split(clabel_sep)[0]] ]  
                   if ctrl ],
        
        chips2 = lambda w: [join(bam_dir,chip+".sorted.bam")
                   for chip in
                   groups[w.clabel.split(clabel_sep)[1]] ],
        
        inputs2 = lambda w: [join(bam_dir,ctrl+".sorted.bam")
                   for ctrl in
                   [ chip2input[c] for c in groups[w.clabel.split(clabel_sep)[1]] ]
                   if ctrl ],
    output:
        join(pepr_dir,'{clabel}__PePr_chip1_peaks.bed'),
        join(pepr_dir,'{clabel}__PePr_chip2_peaks.bed'),
        
    run:
        group1, group2 = wildcards.clabel.split( clabel_sep )
        
        chips1 = groups[group1]
        chips2 = groups[group2]
        
        inputs1 = list(set([ chip2input[c] for c in chips1 ]))
        inputs2 = list(set([ chip2input[c] for c in chips2 ]))
        
        c1 = ','.join([join(bam_dir,c+'.sorted.bam') for c in chips1 ])
        i1 = ','.join([join(bam_dir,c+'.sorted.bam') for c in inputs1])
        c2 = ','.join([join(bam_dir,c+'.sorted.bam') for c in chips2 ])
        i2 = ','.join([join(bam_dir,c+'.sorted.bam') for c in inputs2])
        
        '''
        bams = [ *chips1, *chips2, *inputs1, *inputs2 ]
        bams2com_chroms = { b:b+'.sorted.common_chrom.bam' for b in bams }
        
        ibamfns = [k+'.sorted.bam' for k,v in bams2com_chroms.items()]
        obamfns = [v for k,v in bams2com_chroms.items()]
        
        normalize_bam_file_chromosomes( ibamfns, obamfns )
        
        #replace the filenames
        chips1 = [ bams2com_chroms[b] for b in chips1 ]
        chips2 = [ bams2com_chroms[b] for b in chips2 ]
        inputs1 = [ bams2com_chroms[b] for b in inputs1 ]
        inputs2 = [ bams2com_chroms[b] for b in inputs2 ]
        
        c1 = ','.join(chips1)
        i1 = ','.join(inputs1)
        c2 = ','.join(chips2)
        i2 = ','.join(inputs2)
        '''
            
        if all(inputs1) and all(inputs2) :
            shell( """
            module load PePr;
            PePr --diff -f bam -c {c1} --chip2 {c2} -i {i1} --input2 {i2} -n {wildcards.clabel} --output-directory {pepr_dir}
            """ )
        else :
            shell( """
            module load PePr;
            PePr --diff -f bam -c {c1} --chip2 {c2} -n {wildcards.clabel} --output-directory {pepr_dir}
            """ )
        
rule ChIPQC:
    params:
        rname='pl:ChIPQC:{caller}',
        genome = config['project']['annotation'],
        
    input:
        peakinfo_fn = 'peakcall.tab',
        peak_files = lambda w: [ join(caller2dir[w.caller], 
                         g, 
                         c+caller2suffix[w.caller]
                        ) for g in groups 
                              for c in samples ]
        
    output:
        join('ChIPQC','{caller}','ChIPQCreport.html')
        
    run:
        shell( 'mkdir -p ChIPQC/{caller}'.format(caller=wildcards.caller) )
        chipqcin_fn = 'chipqctab_{caller}.in'.format( caller=wildcards.caller )
        
        #to avoid problems we need to correct names
        cnames = {} #this is set of corrected names
        for s in config['project']['groups'] : #group2chips
            if s[0].isalpha :
                cnames[s] = s
            else :
                cnames[s] = 'v'+s
        for s in samples :
            if s[0].isalpha :
                cnames[s] = s
            else :
                cnames[s] = 'v'+s
        for s in inputs :
            if s[0].isalpha :
                cnames[s] = s
            else :
                cnames[s] = 'v'+s
               
        suffix = caller2suffix[wildcards.caller]
        peak_dir = caller2dir[wildcards.caller]
        #chipqcin = 'ChIPQC{c}.in'.format( c=wildcards.caller )

        sfp = StringIO()
        groups = {}
        for l in open(input.peakinfo_fn) :
            c, i, g = l.split('\t')
            g=g.strip()
            if g in groups :
                groups[g].append([c,i])
            else :
                groups[g]=[[c,i]]

        print(
            'SampleID', #group id 
            'Condition', #sample id
            'Replicate',
            'bamReads',
            'ControlID',
            'bamControl',
            'Peaks',
            file=sfp, sep='\t' 
        )

        for g, v in groups.items() :
            for i,(chp,inp) in enumerate(v) :
                if inp :
                    print(
                        cnames[chp], 
                        cnames[g],
                        i+1, 
                        join(bam_dir,chp+'.sorted.bam'),
                        cnames[inp],
                        join(bam_dir,inp+'.sorted.bam'),
                        join(peak_dir, g, chp+suffix),
                        file=sfp, sep="\t"
                    )
                else :
                    print(
                        cnames[chp], 
                        cnames[g],
                        i+1, 
                        join(bam_dir,chp+'.sorted.bam'),
                        '',
                        '',
                        join(peak_dir, g, chp+suffix),
                        file=sfp, sep="\t"
                    )

        ofp = open( chipqcin_fn, 'w' )
        ofp.write( sfp.getvalue() )
        ofp.close()
            
        R_string = """
        library(ChIPQC)
        samples <- read.table('{tab}', header=1, sep="\t")
        result <- ChIPQC(samples, annotation='{genome}')
        ChIPQCreport( result, reportName="ChIPQCreport", reportFolder="ChIPQC/{caller}" )
        """.format( tab=chipqcin_fn, caller=wildcards.caller, genome=params.genome )
        
        rscript_fn = "chipqc_run_{c}.R".format( c=wildcards.caller )
        of=open(rscript_fn,'w')
        of.write(R_string)
        of.close()
        
        shell( '''
            module load R
            Rscript {r}
            '''.format(r=rscript_fn)
        )
        

rule MEMEChIP:
    input:
        join(macsn_dir,'{group}',"{name}"+macsn_suffix)
        #join(macsb_dir,"{name}_peaks.xls"),
        #join(sicer_dir,"{name}_broadpeaks.bed"),
        #join(sicer_rmdup_dir,"{name}_broadpeaks.bed"),
    params:
        rname='pl:MEMEChIP:{group}:{name}',
        sorted_peak = "{name}.sorted",
        sorted_fa = "{name}.sorted.fa",
        genome_fasta= config['references']['ChIPseq']['GENOME'],
    threads: 16
    output:
        join(memechip_dir,'{group}',"{name}_memechip")
        
    run :
        #shell("""
        #sort -k9,9gr {i} | head -n 500 > {sorted_peak} ;
        #""".format(i=input, sorted_peak=params.sorted_peak )
        #)
        
        content = []
        for l in open(input[0]):
            line = l.split()
            if line :
                content.append( (float(line[8]), l) )
        
        fp = open( params.sorted_peak, 'w' )
        for i, (_,l) in enumerate(sorted(content)) :
            if i < 500 :
                print( *(l.split()[:4]), sep="\t", file=fp )
            else :
                break
        fp.close()
        
        shell("""
        module load bedtools ;
        bedtools getfasta -fi {genome_fa} -bed {sorted_peak} -fo {sorted_fa} ;
        """.format(genome_fa = params.genome_fasta,
                   sorted_peak = params.sorted_peak,
                   sorted_fa = params.sorted_fa,
            )
        )
        
        shell("""
        module load meme ;
        meme-chip --oc {o} -dna {sorted_fa} ;
        """.format(o = output, 
                   sorted_fa = params.sorted_fa,
              )
        )
        
rule MACS2_narrow:
    output:
        join( macsn_dir, '{group}', '{name}'+macsn_suffix )
        
    input: 
        join(bam_dir,"{name}"+bam_suffix)
        
    params: 
        rname='pl:MACS2_narrow:{group}:{name}',
        gsize=config['project']['gsize'],
        ctrl = lambda w : join(bam_dir,chip2input[w.name] + bam_suffix),
        genome=config['project']['annotation'],
        
    shell:
        """
        module load macs

        if [ -n "{params.ctrl}" ] 
        then
            macs2 callpeak -t {input} -c {params.ctrl} -f BAM -g {params.gsize} -n {wildcards.name} --outdir {macsn_dir}/{wildcards.group} -B -q 0.01 ;
        else                
            macs2 callpeak -t {input} -f BAM -g {params.gsize} -n {wildcards.name} --outdir {macsn_dir}/{wildcards.group} -B -q 0.01;
        fi

        #module load ceas
        #ceas -g /fdb/CEAS/{params.genome}.refGene -b {output}
        """

rule MACS2_broad:
    input: 
        join( bam_dir, "{name}.sorted.bam")
    output:
        join( macsb_dir, '{group}', '{name}'+macsb_suffix )
    params: 
        rname='pl:MACS2_broad:{group}:{name}',
        gsize=config['project']['gsize'],
        ctrl = lambda w : config['project']['peaks']['inputs'][w.name],
        genome=config['project']['annotation'],
    shell:
        """
        module load macs
        ctrl="{params.ctrl}"

        if [ -n "$ctrl" ] 
        then
            macs2 callpeak -t {input} -c {bam_dir}/${{ctrl}}.sorted.bam -f BAM -g {params.gsize} -n {wildcards.name} --outdir {macsb_dir}/{wildcards.group} --broad --broad-cutoff 0.1 ;
        else 
            macs2 callpeak -t {input} -f BAM -g {params.gsize} -n {wildcards.name} --outdir {macsb_dir}/{wildcards.group} --broad --broad-cutoff 0.1;
        fi

        #module load ceas
        #ceas -g /fdb/CEAS/{params.genome}.refGene -b {output}
        """

rule SICER :       
    input:
        join(bam_dir,"{name}.sorted.bam")
    output:
        #lambda w: join(sicer,groups[w.name],w.name+sicer_suffix),
        join( sicer_dir, '{group}', '{name}'+sicer_suffix )
    params:
        rname = "pl:SICER:{group}:{name}",
        ctrl = lambda w : config['project']['peaks']['inputs'][w.name],
        SICERDIR = "/usr/local/apps/sicer/1.1",
        genome = config['project']['annotation'],
    run:
        if params.ctrl: 
            shell( 
            """
            module load bedtools
            #module load ceas

            mkdir -p {sicer_dir}/{wildcards.group}/{wildcards.name}
            cd {sicer_dir}/{wildcards.group}/{wildcards.name}

            bamToBed -i ../../../{input} > {wildcards.name}.bed
            bamToBed -i ../../../{bam_dir}/{params.ctrl}.sorted.bam > {params.ctrl}.bed
            """)
            
            #fix_chrom( join(sicer_dir,
            #               wildcards.group,
            #               wildcards.name,
            #               wildcards.name+'.bed') )
            
            #fix_chrom( join(sicer_dir, 
            #               wildcards.group, 
            #               wildcards.name, 
            #               params.ctrl+'.bed') )
            
            shell("""
            cd {sicer_dir}/{wildcards.group}/{wildcards.name}
            module load sicer
            bash {params.SICERDIR}/SICER.sh ./ {wildcards.name}.bed {params.ctrl}.bed ./ {params.genome} 1 300 300 0.75 600 1E-2
            #ceas -g /fdb/CEAS/{params.genome}.refGene -b {wildcards.name}-W300-G600-FDR1E-2-island.bed -w {wildcards.name}-W300-normalized.wig
            ln {wildcards.name}-W300-G600-islands-summary-FDR1E-2 ../{wildcards.name}_broadpeaks.bed
            """)

        else : 
            shell( """
            module load sicer
            module load bedtools
            #module load ceas

            mkdir -p {sicer_dir}/{wildcards.group}/{wildcards.name}
            cd {sicer_dir}/{wildcards.group}/{wildcards.name}

            bamToBed -i ../../../{input} > {wildcards.name}.bed
            """)
            
            #fix_chrom( join(sicer_dir, 
            #               wildcards.group, 
            #               wildcards.name, 
            #               wildcards.name+'.bed') )
            
            shell ("""
            cd {sicer_dir}/{wildcards.group}/{wildcards.name}
            bash {params.SICERDIR}/SICER-rb.sh ./ {wildcards.name}.bed ./ {params.genome} 1 300 300 0.75 600 100
            #ceas -g /fdb/CEAS/{params.genome}.refGene -b {wildcards.name}-W300-G600-E100.scoreisland
            ln {wildcards.name}-W300-G600-E100.scoreisland ../{wildcards.name}_broadpeaks.bed
            """)
            
            

                        
