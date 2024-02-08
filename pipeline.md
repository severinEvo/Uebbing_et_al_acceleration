### Prerequisits:

- perl
- python
- bedtools
- [UCSC genome browser utilities](http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads): liftOver, bigWigToBedGraph, bedToBigBed
- [PHAST Utilities](http://compgen.cshl.edu/phast/): maf_parse
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

### Additioal annotation files

liftOver files

	wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg38.over.chain.gz
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm10.over.chain.gz

phastCons scores

	wget --no-check-certificate https://genome.senckenberg.de//data/hg38/phastConsEline_hg38_multiz120Mammals.bw
	bigWigToBedGraph phastConsEline_hg38_multiz120Mammals.bw phastConsEline_hg38_multiz120Mammals.bdg

Branch list for R

	R CMD BATCH ../scripts/branch-list.R
	cd ..

# Generate PRE dataset

	mkdir data; cd data/
	wget https://bds.mpi-cbg.de/hillerlab/120MammalAlignment/Human120way/data/conservation/phastConsElements_hg38_multiz120Mammals.bed.gz
	bedtools intersect -v -a phastConsElements_hg38_multiz120Mammals.bed.gz -b ../annotations/gencode_exon.v40.annotation.gtf.gz | bedtools intersect -v -a - -b ../annotations/gencode_promotor-1kup0kdn.v40.annotation_ENS.bed.gz | bedtools merge -i - -d 9 | perl -lane 'print if($F[2]-$F[1]>50);' | grep -P "^chr[0-9XY]+\t" >phastConsElements_hg38_multiz120Mammals_noExon_merge10_l49.bed; gzip phastConsElements_hg38_multiz120Mammals_noExon_merge10_l49.bed
	bedtools intersect -v -a phastConsElements_hg38_multiz120Mammals_noExon_merge10_l49.bed.gz -b <(zcat ../annotations/hg38_hg38_all_chain_merge.bed.gz ../annotations/repeats_merge.bed.gz) | bedtools intersect -v -a - -b ../annotations/gencode.v40.2wayconspseudos.gtf.gz >phastConsElements_hg38_multiz120Mammals_filtered.bed; gzip phastConsElements_hg38_multiz120Mammals_filtered.bed

	wget https://api.wenglab.org/screen_v13/fdownloads/GRCh38-cCREs.bed
	bedtools intersect -v -a GRCh38-cCREs.bed -b ../annotations/gencode_exon.v40.annotation.gtf.gz | bedtools intersect -v -a - -b ../annotations/gencode_promotor-1kup0kdn.v40.annotation_ENS.bed.gz | bedtools intersect -v -a - -b <(zcat ../annotations/hg38_hg38_all_chain_merge.bed.gz ../annotations/repeats_merge.bed.gz) | bedtools intersect -v -a - -b ../annotations/gencode.v40.2wayconspseudos.gtf.gz | perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[0]:$F[1]-$F[2]";' >GRCh38-cCREs_filtered_named.bed

Mouse ChIP data

	mkdir mouse_chip; cd mouse_chip/
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

	mkdir MAFs; cd MAFs/
	perl ../../scripts/prepRun_extractMafs.pl
	cd ../

This will generate a job file with one job per line. Run jobs as appropriate for your computing cluster. Scripts here are written for use with Slurm and Dead Simple Queue (https://github.com/ycrc/dsq).

### Detect accelerated sequence substitution on individual branches using phyloP

	mkdir phyloP_out; cd phyloP_out/
	perl ../../scripts/prepRun_phyPacc.pl # Run jobs
	cd ../

### Branch-wise filters against GC-biased gene conversion and indels

Run phastBias to detect likely GC-biased gene conversion

	mkdir phastBias_out; cd phastBias_out/
	../../scripts/prepRun_phastBias.pl # Run jobs
	cd ../

Individual sequence-based indel filter

	cd MAFs/
	../../scripts/prepRun_qualFilter.pl # Run jobs
	cd ../

Then run 

	mkdir qualFilter; cd qualFilter/
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

# Collect dataset and generate Figure 1

Fig 1A was generated by plotting annotations/tree_hg38_120mammals_named.mod in FigTree v1.4.4

	R CMD BATCH collect-phyloP.R

Raw output figures were edited in Inkscape
