### Prerequisits:

- perl
- python
- bedtools
- [UCSC genome browser utilities](http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads): liftOver, bigWigToBedGraph, bedToBigBed, towBitToFa
- [PHAST Utilities](http://compgen.cshl.edu/phast/): maf_parse
- MEME suite
- R

# Generate annotations

### Retrieve and parse repeat annotations

	cd annotations/
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
	zcat rmsk.txt.gz | cut -f 6-8 >rmsk.bed
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz
	zcat simpleRepeat.txt.gz | cut -f 2-4 >simpleRepeat.bed
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz
	zcat genomicSuperDups.txt.gz | cut -f 2-4 >genomicSuperDups.bed
	cat rmsk.bed simpleRepeat.bed genomicSuperDups.bed | sort -k1,1 -k 2,2n | bedtools merge -i - >repeats_merge.bed; gzip repeats_merge.bed

	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsSelf/hg38.hg38.all.chain.gz
	cat hg38.hg38.all.chain.gz | grep -P "^chain " | perl -lane 'next if($F[5] == $F[10] && $F[6] == $F[11]); print "$F[2]\t$F[5]\t$F[6]\n$F[7]\t$F[10]\t$F[11]";' >hg38.hg38.all.chain.bed
	sort -k 1,1 -k 2,2n hg38.hg38.all.chain.bed | uniq | bedtools merge -i - >hg38_hg38_all_chain_merge.bed; gzip hg38_hg38_all_chain_merge.bed

Pseudogene annotations

	wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.2wayconspseudos.gtf.gz

### Generate exon and promoter annotations

	wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz
	zcat gencode.v40.annotation.gtf.gz | grep -P "\texon\t" >gencode_exon.v40.annotation.gtf; gzip gencode_exon.v40.annotation.gtf
	../scripts/extract-promotor.pl -i -c hg38.chrom.sizes <(zcat gencode.v40.annotation.gtf.gz) >gencode_promotor-5kup1kdn.v40.annotation_ENS.bed; gzip gencode_promotor-5kup1kdn.v40.annotation_ENS.bed
	../scripts/extract-promotor.pl -i -c hg38.chrom.sizes -5 1000 -3 0 <(zcat gencode.v40.annotation.gtf.gz) >gencode_promotor-1kup0kdn.v40.annotation_ENS.bed; gzip gencode_promotor-1kup0kdn.v40.annotation_ENS.bed

### TFBS motifs

JASPAR

	mkdir -p meme-database; cd meme-database
	wget wget https://jaspar.elixir.no/download/data/2020/CORE/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.zip; unzip JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.zip

HOCOMOCO

	wget https://hocomoco11.autosome.org/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme
	wget https://hocomoco11.autosome.org/final_bundle/hocomoco11/core/MOUSE/mono/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme
	cd ..

Motifs from Jolma et al. Cell 2013, extracted from Table S3, are included in GitHub data

### Miscellaneous annotation files

liftOver files

	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm10.over.chain.gz
	wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz
	wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg38.over.chain.gz

phastCons scores

	wget --no-check-certificate https://genome.senckenberg.de//data/hg38/phastConsEline_hg38_multiz120Mammals.bw
	bigWigToBedGraph phastConsEline_hg38_multiz120Mammals.bw phastConsEline_hg38_multiz120Mammals.bdg

Branch list for R, and a second version to match with Zoonomia species (Fig. 3C)

	Rscript ../scripts/branch-list.R
	Rscript ../scripts/zoonomia-branch-list.R


gnomAD data

	wget https://datasetgnomad.blob.core.windows.net/dataset/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
	gunzip -c gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz >gnomad.v2.1.1.lof_metrics.by_gene.txt; gzip gnomad.v2.1.1.lof_metrics.by_gene.txt
	cd ..

# Generate PRE dataset

	mkdir -p data; cd data/
	wget https://bds.mpi-cbg.de/hillerlab/120MammalAlignment/Human120way/data/conservation/phastConsElements_hg38_multiz120Mammals.bed.gz
	bedtools intersect -v -a phastConsElements_hg38_multiz120Mammals.bed.gz -b ../annotations/gencode_exon.v40.annotation.gtf.gz | bedtools intersect -v -a - -b ../annotations/gencode_promotor-1kup0kdn.v40.annotation_ENS.bed.gz | bedtools merge -i - -d 9 | perl -lane 'print if($F[2]-$F[1]>50);' | grep -P "^chr[0-9XY]+\t" >phastConsElements_hg38_multiz120Mammals_noExon_merge10_l49.bed; gzip phastConsElements_hg38_multiz120Mammals_noExon_merge10_l49.bed
	bedtools intersect -v -a phastConsElements_hg38_multiz120Mammals_noExon_merge10_l49.bed.gz -b <(zcat ../annotations/hg38_hg38_all_chain_merge.bed.gz ../annotations/repeats_merge.bed.gz) | bedtools intersect -v -a - -b ../annotations/gencode.v40.2wayconspseudos.gtf.gz >phastConsElements_hg38_multiz120Mammals_filtered.bed; gzip phastConsElements_hg38_multiz120Mammals_filtered.bed

	wget https://api.wenglab.org/screen_v13/fdownloads/GRCh38-cCREs.bed
	bedtools intersect -v -a GRCh38-cCREs.bed -b ../annotations/gencode_exon.v40.annotation.gtf.gz | bedtools intersect -v -a - -b ../annotations/gencode_promotor-1kup0kdn.v40.annotation_ENS.bed.gz | bedtools intersect -v -a - -b <(zcat ../annotations/hg38_hg38_all_chain_merge.bed.gz ../annotations/repeats_merge.bed.gz) | bedtools intersect -v -a - -b ../annotations/gencode.v40.2wayconspseudos.gtf.gz | perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[0]:$F[1]-$F[2]";' >GRCh38-cCREs_filtered_named.bed

Mouse ChIP data

	mkdir -p mouse_chip; cd mouse_chip/
	wget https://api.wenglab.org/screen_v13/fdownloads/mm10-ccREs.bed
	liftOver -minMatch=.7 mm10-ccREs.bed ../../annotations/mm10ToHg38.over.chain.gz mm10TOhg38-ccREs.bed unMapped
	for i in $(cat files.txt); do wget $i; done
	for i in $(ls *_segments.bed.gz | perl -lne '$_ =~ s/.bed.gz//; print;'); do zcat $i\.bed.gz | grep -vP 'U1[0-9]\t' | cut -f 1-3 | sort -k1,1 -k2,2n | bedtools merge -i - >$i\_merge.bed; done # filter out inactive chromatin regions (U10-U15)
	zcat *_segments_merge.bed | sort -k1,1 -k2,2n | bedtools merge -i - >all_chromHMM.bed
	cat all_chromHMM.bed <(cut -f 1-3 mm10-ccREs.bed) | sort -k1,1 -k2,2n | bedtools merge -i - >all_mouse.bed

chromHMM regions are larger than cCREs and therefore lift worse, so liftOver cCREs to mouse instead:

	liftOver -minMatch=.7 ../GRCh38-cCREs_filtered_named.bed ../../annotations/hg38ToMm10.over.chain.gz GRCh38TOmm10-cCREs_filtered_named.bed unMapped

Intersect lifted human cCREs with mouse ChIP-type data and revert into hg38 coordinates using the name field:

	bedtools intersect -u -a GRCh38TOmm10-cCREs_filtered_named.bed -b all_mouse.bed | perl -lane '@a = split /:|-/, $F[3]; print "$a[0]\t$a[1]\t$a[2]";' >../GRCh38-cCREs_mouse-cCRE-chromHMM-intersect.bed
	cd ../; gzip GRCh38-cCREs_mouse-cCRE-chromHMM-intersect.bed

### Annotate PRE list with PRE type, sequence space type, phastCons score

PRE type

	cat <(zcat phastConsElements_hg38_multiz120Mammals_filtered.bed.gz | perl -lne 'print "$_\tpCE";') <(zcat GRCh38-cCREs_mouse-cCRE-chromHMM-intersect.bed.gz | perl -lane '$" = "\t"; print "@F[0..2]\tcCRE";') | sort -k1,1 -k2,2n >type.bed

Promoter-proximal

	bedtools intersect -u -a type.bed -b ../annotations/gencode_promotor-5kup1kdn.v40.annotation_ENS.bed.gz | perl -lane '$" = "\t"; print "@F[0..3]\tpromoter-proximal";' >promoter

Intronic

	zcat ../annotations/gencode.v40.annotation.gtf.gz | grep transcript | bedtools intersect -u -a type.bed -b - | bedtools intersect -v -a - -b ../annotations/gencode_promotor-5kup1kdn.v40.annotation_ENS.bed.gz | perl -lane '$" = "\t"; print "@F[0..3]\tintronic";' >intronic

Intergenic

	zcat ../annotations/gencode.v40.annotation.gtf.gz | grep transcript | bedtools intersect -v -a type.bed -b - | bedtools intersect -v -a - -b ../annotations/gencode_promotor-5kup1kdn.v40.annotation_ENS.bed.gz | perl -lane '$" = "\t"; print "@F[0..3]\tintergenic";' >intergenic

	cat promoter intronic intergenic | sort -k1,1 -k2,2n >sequence.bed

phastCons

	bedtools map -a sequence.bed -b ../annotations/phastConsEline_hg38_multiz120Mammals.bdg -c 4 -o mean >PREs.bed; gzip PREs.bed

### Download and parse MAF files

	mkdir -p MAFs; cd MAFs/
	perl ../../scripts/prepRun_extractMafs.pl
	cd ../

This will generate a job file with one job per line. Run jobs as appropriate for your computing cluster. Scripts here are written for use with Slurm and Dead Simple Queue (https://github.com/ycrc/dsq).

### Detect accelerated sequence substitution on individual branches using phyloP

	mkdir -p phyloP_out; cd phyloP_out/
	perl ../../scripts/prepRun_phyPacc.pl # Run jobs
	cd ../

### Branch-wise filters against GC-biased gene conversion and indels

Run phastBias to detect likely GC-biased gene conversion

	mkdir -p phastBias_out; cd phastBias_out/
	../../scripts/prepRun_phastBias.pl # Run jobs
	cd ../

Individual sequence-based indel filter

	cd MAFs/
	../../scripts/prepRun_qualFilter.pl # Run jobs
	cd ../

Then run 

	mkdir -p qualFilter; cd qualFilter/
	cat *_use.bed | sort -u | sort -k1,1 -k2,2n >qualFilter_useTips.bed
	gzip qualFilter_useTips.bed
	python ../../scripts/qualFilter-step2.py qualFilter_useTips.bed.gz ../../annotations/branchlist.pyData qualFilter_useBranches.tsv; gzip qualFilter_useBranches.tsv
	zcat qualFilter_useBranches.tsv.gz | perl -lane '@p = split /:|-/, $F[0]; $" = "\t"; print "@p\t$F[1]";' | sort -k1,1 -k2,2n >qualFilter_useBranches.bed; gzip qualFilter_useBranches.bed
	cd ../

### Collect accelerated PRE dataset

Combine pCEs and cCREs (per branch)

	cd phyloP_out/
	for i in $(ls *Es_phyP_ACC.out | cut -f 2 -d '_' | sort -u); do echo $i; cat ../../annotations/header >$i\_phyloP_ACC.tsv; cat chr*_$i\_*Es_phyP_ACC.out | grep -vP "^#chr\t" | perl -lane '$" = "\t"; print "@F[0..2]\t@F[4..8]\t'$i'";' | sort -k1,1 -k2,2n >>$i\_phyloP_ACC.tsv; gzip -f $i\_phyloP_ACC.tsv; done

Filter out results in phastBias tracks

	for i in $(ls *_phyloP_ACC.tsv.gz | cut -f 2 -d '_' | sort -u); do echo $i; zcat $i\_phyloP_ACC.tsv.gz | head -1 >$i\_phyloP-phastBias.tsv; zcat $i\_phyloP_ACC.tsv.gz | grep -vP "^chr\t" | bedtools intersect -v -a - -b ../phastBias_out/$i\_phastBias.gff.gz >>$i\_phyloP-phastBias.tsv; gzip -f $i\_phyloP-phastBias.tsv; done
	cd ../../

# Gene-PRE connectivity data

### Jung et al. Nature Genetics 2019 (original data hg19)
NPC_P-O-contacts.bed extracted from Table S3

	mkdir -p data/PRE-gene-dataset/Jung/; cd data/PRE-gene-dataset/Jung/
	liftOver NPC_P-O-contacts.bed ../../../annotations/hg19ToHg38.over.chain NPC_P-O-contacts_hg38.bed unMapped
	perl -lane 'if($_ =~ /;/){@s = split /;/, $F[3]; $" = "\t"; foreach(@s){print "@F[0..2]\t$_";}}else{print;}' <(zcat NPC_P-O-contacts_hg38.bed.gz) >NPC_P-O-contacts_hg38_1perline.bed; gzip NPC_P-O-contacts_hg38_1perline.bed
	Rscript ../../../scripts/Jung-pcHiC.R
	bedtools intersect -wo -a ../../PREs.bed.gz -b NPC_P-O-contacts_hg38_ENSEMBL.bed.gz | perl -lane '$" = "\t"; print "@F[0..2]\t$F[9]\t$F[3]\tNPC-pcHiC";' | sort -u >../connection-collection.bed

### Song et al. Nature 2020 (original data hg38)
Tables download from NeMO Archive (controlled access)

	mkdir -p ../Song/; cd ../Song/
	zcat RG.MAPS.peaks.txt.gz | tail -n +2 | bedtools intersect -wo -a - -b ../../../annotations/gencode_promotor-5kup1kdn.v40.annotation_ENS.bed.gz | cut -f 4-6,13 | sort -u >tmp; zcat RG.MAPS.peaks.txt.gz | tail -n +2 | perl -lane '$" = "\t"; print "@F[3..5]\t@F[0..2]";' | bedtools intersect -wo -a - -b ../../../annotations/gencode_promotor-5kup1kdn.v40.annotation_ENS.bed.gz | cut -f 4-6,10 | sort -u >>tmp; sort -k1,1 -k2,2n tmp >RG_gencode-v40-promotor_intersect.bed
	zcat IPC.MAPS.peaks.txt.gz | tail -n +2 | bedtools intersect -wo -a - -b ../../../annotations/gencode_promotor-5kup1kdn.v40.annotation_ENS.bed.gz | cut -f 4-6,13 | sort -u >tmp; zcat IPC.MAPS.peaks.txt.gz | tail -n +2 | perl -lane '$" = "\t"; print "@F[3..5]\t@F[0..2]";' | bedtools intersect -wo -a - -b ../../../annotations/gencode_promotor-5kup1kdn.v40.annotation_ENS.bed.gz | cut -f 4-6,10 | sort -u >>tmp; sort -k1,1 -k2,2n tmp >IPC_gencode-v40-promotor_intersect.bed
	zcat eN.MAPS.peaks.txt.gz | tail -n +2 | bedtools intersect -wo -a - -b ../../../annotations/gencode_promotor-5kup1kdn.v40.annotation_ENS.bed.gz | cut -f 4-6,13 | sort -u >tmp; zcat eN.MAPS.peaks.txt.gz | tail -n +2 | perl -lane '$" = "\t"; print "@F[3..5]\t@F[0..2]";' | bedtools intersect -wo -a - -b ../../../annotations/gencode_promotor-5kup1kdn.v40.annotation_ENS.bed.gz | cut -f 4-6,10 | sort -u >>tmp; sort -k1,1 -k2,2n tmp >eN_gencode-v40-promotor_intersect.bed
	zcat iN.MAPS.peaks.txt.gz | tail -n +2 | bedtools intersect -wo -a - -b ../../../annotations/gencode_promotor-5kup1kdn.v40.annotation_ENS.bed.gz | cut -f 4-6,13 | sort -u >tmp; zcat iN.MAPS.peaks.txt.gz | tail -n +2 | perl -lane '$" = "\t"; print "@F[3..5]\t@F[0..2]";' | bedtools intersect -wo -a - -b ../../../annotations/gencode_promotor-5kup1kdn.v40.annotation_ENS.bed.gz | cut -f 4-6,10 | sort -u >>tmp; sort -k1,1 -k2,2n tmp >iN_gencode-v40-promotor_intersect.bed
	gzip -f *_gencode-v40-promotor_intersect.bed

	bedtools intersect -wo -a ../../PREs.bed.gz -b RG_gencode-v40-promotor_intersect.bed.gz | perl -lane '$" = "\t"; print "@F[0..2]\t$F[9]\t$F[3]\tRG-PLAC-seq";' | sort -u >>../connection-collection.bed
	bedtools intersect -wo -a ../../PREs.bed.gz -b IPC_gencode-v40-promotor_intersect.bed.gz | perl -lane '$" = "\t"; print "@F[0..2]\t$F[9]\t$F[3]\tIPC-PLAC-seq";' | sort -u >>../connection-collection.bed
	bedtools intersect -wo -a ../../PREs.bed.gz -b eN_gencode-v40-promotor_intersect.bed.gz | perl -lane '$" = "\t"; print "@F[0..2]\t$F[9]\t$F[3]\teN-PLAC-seq";' | sort -u >>../connection-collection.bed
	bedtools intersect -wo -a ../../PREs.bed.gz -b iN_gencode-v40-promotor_intersect.bed.gz | perl -lane '$" = "\t"; print "@F[0..2]\t$F[9]\t$F[3]\tiN-PLAC-seq";' | sort -u >>../connection-collection.bed

### Won et al. Nature 2016 (original data hg19)
Data extracted from Tables S22,S23 (saved as .tsv from .xlsx)

	mkdir -p ../Won/; cd ../Won/
	grep -vP "^chr\t" Won_TSS_CP.tsv | cut -f 1,4,5,7 | sort -k1,1 -k2,2n >Won_TSS_CP_hg19.bed
	liftOver Won_TSS_CP_hg19.bed ~/project/genomes/hg19ToHg38.over.chain Won_TSS_CP_hg38.bed unMapped
	grep -vP "^chr\t" Won_TSS_GZ.tsv | cut -f 1,4,5,7 | sort -k1,1 -k2,2n >Won_TSS_GZ_hg19.bed
	liftOver Won_TSS_GZ_hg19.bed ~/project/genomes/hg19ToHg38.over.chain Won_TSS_GZ_hg38.bed unMapped

	bedtools intersect -wo -a ../../PREs.bed.gz -b /Won_TSS_CP_hg38.bed.gz | perl -lane '$" = "\t"; print "@F[0..2]\t$F[9]\t$F[3]\tCP-HiC";' | sort -u >>../connection-collection.bed
	bedtools intersect -wo -a ../../PREs.bed.gz -b Won_TSS_GZ_hg38.bed.gz | perl -lane '$" = "\t"; print "@F[0..2]\t$F[9]\t$F[3]\tGZ-HiC";' | sort -u >>../connection-collection.bed

### Rajarajan et al. (PsychENCODE) Science 2018 (original data hg19)
clean_master_hiccups_loops_nonsubsampled.txt extracted from Supplemental tables xls

	mkdir -p ../PsychENCODE/; cd ../PsychENCODE/
	for i in $(cut -f 21 clean_master_hiccups_loops_nonsubsampled.txt  | sort -u); do echo $i; grep -P "\t$i\t" clean_master_hiccups_loops_nonsubsampled.txt | bedtools intersect -wo -a - -b ../../../annotations/gencode_promotor-5kup1kdn.v40lift37.annotation_ENS.bed.gz | cut -f 4-6,29 | sort -u >tmp; grep -P "\t$i\t" clean_master_hiccups_loops_nonsubsampled.txt | perl -lane '$" = "\t"; print "@F[3..5]\t@F[0..2]";' | bedtools intersect -wo -a - -b ../../../annotations/gencode_promotor-5kup1kdn.v40lift37.annotation_ENS.bed.gz | cut -f 4-6,10 | sort -u >>tmp; sort -k1,1 -k2,2n tmp >$i\_hiccups_loops_nonsubsampled_hg19.bed; liftOver $i\_hiccups_loops_nonsubsampled_hg19.bed ../../../annotations/hg19ToHg38.over.chain $i\_hiccups_loops_nonsubsampled_hg38.bed unMapped_$i; done; gzip -f *.bed

	bedtools intersect -wo -a ../../PREs.bed.gz -b NPC_hiccups_loops_nonsubsampled_hg38.bed.gz | perl -lane '$" = "\t"; print "@F[0..2]\t$F[9]\t$F[3]\tNPC-HiC";' | sort -u >>../connection-collection.bed
	bedtools intersect -wo -a ../../PREs.bed.gz -b Neu_hiccups_loops_nonsubsampled_hg38.bed.gz | perl -lane '$" = "\t"; print "@F[0..2]\t$F[9]\t$F[3]\tNeu-HiC";' | sort -u >>../connection-collection.bed
	bedtools intersect -wo -a ../../PREs.bed.gz -b Astro_hiccups_loops_nonsubsampled_hg38.bed.gz | perl -lane '$" = "\t"; print "@F[0..2]\t$F[9]\t$F[3]\tAstro-HiC";' | sort -u >>../connection-collection.bed

### 4DNucleosome data

	mkdir -p ../4Dnucleosome/; cd ../4Dnucleosome/

Du, Zheng et al. Nature 2017 (original data in mm9)
First, run hiccups (juicertools), which needs a GPU system. Here, using Google Colab (switch runtime to GPU!):

	!wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar
	!curl -O --user [user token] https://data.4dnucleome.org/files-processed/4DNFII3JV8I1/@@download/4DNFII3JV8I1.hic
	!java -Xmx10g -jar juicer_tools_1.22.01.jar hiccups 4DNFII3JV8I1.hic cortex_results --ignore-sparsity

	mv merged_loops.bedpe 4DNESUQT299T_mouse-P28-cortex-mm9_hiccups-merged-loops.bed

add "chr" to chromosome field

	perl -lane 'if(/^#/){print;}else{$F[0] = "chr$F[0]"; $F[3] = "chr$F[3]"; $" = "\t"; print "@F";}' 4DNESUQT299T_mouse-P28-cortex-mm9_hiccups-merged-loops.bed >tmp; mv tmp 4DNESUQT299T_mouse-P28-cortex-mm9_hiccups-merged-loops.bed

liftOver both anchors from mm9 to mm10

	perl -lane '$" = "_"; print "$F[0]\t$F[1]\t$F[2]\t@F[3..23]";' 4DNESUQT299T_mouse-P28-cortex-mm9_hiccups-merged-loops.bed >tmp
	liftOver tmp ../../../annotations/mm9ToMm10.over.chain.gz tmp_mm10 unMapped
	perl -lane '$" = "_"; @u = split /_/, $F[3]; print "$u[0]\t$u[1]\t$u[2]\t@F[0..2]_@u[3..20]";' tmp_mm10 >tmp_mm9_2
	liftOver tmp_mm9_2 ../../../annotations/mm9ToMm10.over.chain.gz 4DNESUQT299T_mouse-P28-cortex_hiccups-merged-loops_mm10.bed unMapped

liftOver mm10 to hg38 and translate ENSM to ENS gene IDs

	bedtools intersect -wo -a 4DNESUQT299T_mouse-P28-cortex_hiccups-merged-loops_mm10.bed -b ../../../annotations/gencode_promotor-5kup1kdn.vM29.annotation_ENS.bed.gz | cut -f 4-6,28 | sort -u >tmp; perl -lane 'next if(/^#/); $" = "\t"; print "@F[3..5]\t@F[0..2]";' 4DNESUQT299T_mouse-P28-cortex_hiccups-merged-loops_mm10.bed | bedtools intersect -wo -a - -b ../../../annotations/gencode_promotor-5kup1kdn.vM29.annotation_ENS.bed.gz | cut -f 4-6,10 | sort -u >>tmp; sort -k1,1 -k2,2n tmp >cortex_hiccups_loops_mm10.bed; liftOver -minMatch=.1 cortex_hiccups_loops_mm10.bed ~/project/genomes/mm/mm10ToHg38.over.chain.gz cortex_hiccups_loops_hg38.bed unMapped
	Rscript ../../../scripts/4Dnucleosome-cortex-ENS.R

	bedtools intersect -wo -a ../../PREs.bed.gz -b cortex_hiccups_loops_hg38_ENS.bed.gz | perl -lane '$" = "\t"; print "@F[0..2]\t$F[9]\t$F[3]\tP28-mouse-cortex-insitu-HiC";' | sort -u >>../connection-collection.bed

Bertero, Fields et al. Nature Communications 2019 (original data hg38)
First, run hiccups (juicertools), which needs a GPU system. Here, using Google Colab:

	!curl -O --user [user token] https://data.4dnucleome.org/files-processed/4DNFIF1EWFUP/@@download/4DNFIF1EWFUP.hic
	!java -Xmx10g -jar juicer_tools_1.22.01.jar hiccups 4DNFIF1EWFUP.hic heart_results --ignore-sparsity

	mv merged_loops.bedpe 4DNESZFHB53P_human-fetal-heart-hg38_hiccups-merged-loops.bed

add "chr" to chromosome field

	perl -lane 'if(/^#/){print;}else{$F[0] = "chr$F[0]"; $F[3] = "chr$F[3]"; $" = "\t"; print "@F";}' 4DNESZFHB53P_human-fetal-heart-hg38_hiccups-merged-loops.bed >tmp; mv tmp 4DNESZFHB53P_human-fetal-heart-hg38_hiccups-merged-loops.bed

	bedtools intersect -wo -a 4DNESZFHB53P_human-fetal-heart-hg38_hiccups-merged-loops.bed -b ../../../annotations/gencode_promotor-5kup1kdn.v40.annotation_ENS.bed.gz | cut -f 4-6,28 | sort -u >tmp; perl -lane 'next if(/^#/); $" = "\t"; print "@F[3..5]\t@F[0..2]";' 4DNESZFHB53P_human-fetal-heart-hg38_hiccups-merged-loops.bed | bedtools intersect -wo -a - -b ../../../annotations/gencode_promotor-5kup1kdn.v40.annotation_ENS.bed.gz | cut -f 4-6,10 | sort -u >>tmp; sort -k1,1 -k2,2n tmp >heart_hiccups_loops_hg38.bed; gzip *.bed

	bedtools intersect -wo -a ../../PREs.bed.gz -b heart_hiccups_loops_hg38.bed.gz | perl -lane '$" = "\t"; print "@F[0..2]\t$F[9]\t$F[3]\tfetal-heart-DNase-HiC";' | sort -u >>../connection-collection.bed

### Promoter overlap

	cd ../
	bedtools intersect -wo -a ../PREs.bed.gz -b ../../annotations/gencode_promotor-5kup1kdn.v40.annotation_ENS.bed.gz | perl -lane '$" = "\t"; print "@F[0..2]\t$F[9]\t$F[3]\tPromoter-overlap";' | sort -u >>connection-collection.bed

### Collect full connectivity dataset

	cat connection-collection.bed | sort -k1,1 -k2,2n >tmp; mv tmp connection-collection.bed; gzip -f connection-collection.bed

# Assemble data and generate outputs

### Collect dataset and generate Figure 1

Fig. 1A was generated by plotting annotations/tree_hg38_120mammals_named.mod in FigTree v1.4.4.

For Fig. 1B,C, S3, run

	Rscript scripts/collect-phyloP.R

All raw output figures were edited in Inkscape

### Figure 2

	Rscript scripts/constraint.R

### Generate and analyze data for TFBS analysis of Fig. 3B

	for i in $(zcat data/PRE-gene-dataset/PREs_qualFiltered.bed.gz | cut -f 1 | sort -u); do echo $i; zcat data/PRE-gene-dataset/PREs_qualFiltered.bed.gz | cut -f 1-3 | grep -P "^$i\\t" >data/MAFs/tmp_mafs/$i\/PREs_$i\.bed; done

Download all 120 genomes, transform into fasta, make Markov models, map motifs to genomes

	perl scripts/prepRun_download-genomes.pl # Run download-prep-genomes.txt
	perl scripts/prepRun_fimoMap.pl # Run fimo_job_*.txt jobs

Generate one bed file with all mapped TFBSs per genome

	perl scripts/prepRun_allbed.pl

Extract one MAF file for each PRE to translate PREs to other genomes

	perl scripts/prepRun_PRE-MAFs.pl # Run parseMAFs-job.txt
	perl scripts/prepRun_translate.pl # Run perl-translate-job.txt

Lists of all PREs per chr, used in following step

	cd data/phyloP_out/; for i in $(zcat phyloP-all.tsv.gz | grep -vP "^chr\t" | cut -f 1 | sort -u); do echo $i; zcat phyloP-all.tsv.gz | grep -P "$i\\t" | perl -lane 'print "$F[0]:$F[1]-$F[2]";' | sort -u >../meme/$i\/phyloP-$i\_all.lst; done

Intersect mapped TFBSs with PREs in their respective genome coordinates

	perl scripts/prep_fimo-intersect.pl # Run chr*_fimo-intersect_job.txt

Calculate distance measures for TFBS sets in accelerated and non-accelerated PRE-branch combinations

	Rscript scripts/dist-sign.R
	Rscript scripts/dist-nsign.R

Analyze and create Fig. 3B

	Rscript scripts/dist-analysis.R

### Integrate and analyze OCR data for Fig. 3C

Download multispecies model OCR data

	mkdir -p data/TACIT; cd data/TACIT
	wget http://daphne.compbio.cs.cmu.edu/files/ikaplow/TACITSupplement/PredictionMatrices/cortex_ocrs_filtered_named.txt
	wget http://daphne.compbio.cs.cmu.edu/files/ikaplow/TACITSupplement/PredictionMatrices/liver_ocrs_filtered_named.txt
	wget http://daphne.compbio.cs.cmu.edu/files/ikaplow/TACITSupplement/PredictionMatrices/pv_ocrs_filtered_named.txt

Transform into BED files for intersecting

	perl -lane 'next if(/Assembly_Chr/); @d = split /\./, $F[0]; @u = split /_/, $d[0]; print "$u[1]\t$d[1]\t$d[2]\t$F[0]";' cortex_ocrs_filtered_named.txt >cortex_ocrs_filtered_named.bed
	perl -lane 'next if(/Assembly_Chr/); @d = split /\./, $F[0]; @u = split /_/, $d[0]; print "$u[1]\t$d[1]\t$d[2]\t$F[0]";' liver_ocrs_filtered_named.txt >liver_ocrs_filtered_named.bed
	perl -lane 'next if(/Assembly_Chr/); @d = split /\./, $F[0]; @u = split /_/, $d[0]; print "$u[1]\t$d[1]\t$d[2]\t$F[0]";' pv_ocrs_filtered_named.txt >pv_ocrs_filtered_named.bed

Lift over from mouse to human

	liftOver -bedPlus=3 -minMatch=.25 cortex_ocrs_filtered_named.bed ../../annotations/mm10ToHg38.over.chain.gz cortex_ocrs_filtered_named_hg38.bed unMapped
	liftOver -bedPlus=3 -minMatch=.25 liver_ocrs_filtered_named.bed ../../annotations/mm10ToHg38.over.chain.gz liver_ocrs_filtered_named_hg38.bed unMapped
	liftOver -bedPlus=3 -minMatch=.25 pv_ocrs_filtered_named.bed ../../annotations/mm10ToHg38.over.chain.gz pv_ocrs_filtered_named_hg38.bed unMapped

Intersect PREs with OCRs (including file headers)

	echo -e "chr\tstart\tend\tAssembly_Chr.Start.End.SummitOffset" >cortex_ocrs_filtered_named_hg38_CRE-intersect.tsv
	bedtools intersect -wo -f .5 -b cortex_ocrs_filtered_named_hg38.bed -a ../PRE-gene-dataset/PREs_qualFiltered.bed.gz | cut -f 1-3,10 >>cortex_ocrs_filtered_named_hg38_PRE-intersect.tsv

	echo -e "chr\tstart\tend\tAssembly_Chr.Start.End.SummitOffset" >liver_ocrs_filtered_named_hg38_CRE-intersect.tsv
	bedtools intersect -wo -f .5 -b liver_ocrs_filtered_named_hg38.bed -a ../PRE-gene-dataset/PREs_qualFiltered.bed.gz | cut -f 1-3,10 >>liver_ocrs_filtered_named_hg38_PRE-intersect.tsv

	echo -e "chr\tstart\tend\tAssembly_Chr.Start.End" >pv_ocrs_filtered_named_hg38_CRE-intersect.tsv
	bedtools intersect -wo -f .5 -b pv_ocrs_filtered_named_hg38.bed -a ../PRE-gene-dataset/PREs_qualFiltered.bed.gz | cut -f 1-3,10 >>pv_ocrs_filtered_named_hg38_PRE-intersect.tsv

Analyze and output Fig. 3C

	Rscript scripts/TACIT-multispecies.R

### Figure 4

	Rscript scripts/hyperGO.R # Fig. 4A-D, S4
	Rscript scripts/pathway-genes.R # Fig. 4E,F

### Figure 6

	Rscript scripts/lacZ.R

Fig. 6D is a subpart of Fig. S5
