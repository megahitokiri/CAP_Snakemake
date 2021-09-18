#--------------------------------------------------------------------------------
# MAIN CONFIGURATION PARAMETERS (INTERNAL CONFIG FILE)
#Cancel all jobs if fail: squeue -u $USER -h | awk '{print $1}' | xargs scancel
#--------------------------------------------------------------------------------

CHRS = '1 2 3 4 5 6 7 8 9 10 11'.split()
PROJECT = "PK_hap2"
REFERENCE = "PK_hap2.09082021.fasta"
NANOPORE_FASTQ = "FAP02733_Porechop_No_Trim_Final.fastq.gz"
AED_FILTER = '0.6 0.7 0.8'.split()

#--------------------------------------------------------------------------------
# TargetRule FINAL_GFF3
#--------------------------------------------------------------------------------

rule FINAL_GFF3:
	input:
#		"Ref/Project_Main.fasta",
		expand("Ref/scaffold_{Chrs}.fasta",Chrs = CHRS),
#		expand("EDTA_Files/scaffold_{Chrs}.fasta.mod.EDTA.TElib.fa",Chrs = CHRS),
#		expand("EDTA_Files/scaffold_{Chrs}.fasta.mod.EDTA.intact.gff3",Chrs = CHRS),
#		expand("EDTA_Files/scaffold_{Chrs}.masked.fasta",Chrs = CHRS),
#		"Minimap_Aligned/minimap2.sorted.MAPQ20.dedup.bam",
#		"Minimap_Aligned/minimap2.sorted.MAPQ20.dedup.bam.bai",
		expand("Sorted_Chromosomes/Nanopore.chr{Chrs}.fasta",Chrs = CHRS),
#		"Sorted_Chromosomes/Minor_Scaffolds.fasta",
		expand("Maker_Files/Chr{Chrs}/scaffold_{Chrs}.maker.output/scaffold_{Chrs}_master_datastore_index.log",Chrs = CHRS),
		expand("Maker_Files/Chr{Chrs}/MAKER_Filtered.scaffold_{Chrs}.all.AED_{AED_filter}.gff3",Chrs = CHRS,AED_filter=AED_FILTER),
		expand("Post_Maker_Files/MAKER_ORF_Filtered_{Project}.scaffold_{Chrs}.AED_{AED_filter}.gff3",Project=PROJECT,Chrs = CHRS,AED_filter=AED_FILTER),
		expand("FINAL_ANNOTATION/FINAL_{Project}.AED_{AED_filter}.sorted.gff3",Project=PROJECT,AED_filter=AED_FILTER),
		expand("COGNATE/COGNATE_{Project}.AED_{AED_filter}/COGNATE_{Project}.AED_{AED_filter}_00-analyzed_transcripts.fa",Project=PROJECT,AED_filter=AED_FILTER),
		expand("COGNATE/COGNATE_{Project}.AED_{AED_filter}/COGNATE_{Project}.AED_{AED_filter}_01-summary.tsv",Project=PROJECT,AED_filter=AED_FILTER),
		
#--------------------------------------------------------------------------------
# Init: Initializing files and folder
#--------------------------------------------------------------------------------
rule Init:
	input:
		reference={REFERENCE},
		maker_files="Maker_Files.zip",
		Assembly_spliter="Assembly_Chr_splitter.R",
	output:
		"Ref/Project_Main.fasta",
	log:
		"logs/Init.log"
	shell:
		"""
		snakemake --dag | dot -Tsvg > dag.svg
		mkdir logs
		mkdir EDTA_Files
		mkdir Minimap_Aligned
		mkdir Sorted_Chromosomes
		mkdir Nanopore_FASTQ
		mkdir Ref
		unzip -o Maker_Files.zip
		mkdir Post_Maker_Files
		mkdir FINAL_ANNOTATION
		mkdir COGNATE

		cp  {input.reference} {output}

		ml samtools

		#Index FASTA file
		samtools faidx {output} 
		
		ml unload samtools
		echo Init succesfuly run > {log} 2>&1
		"""

#--------------------------------------------------------------------------------
# Chr_splitting: Split the Main Assembly in Chromosomes for easy handling
#--------------------------------------------------------------------------------

rule Chr_splitting:
	input:
		rules.Init.output,
	output:
		"Ref/scaffold_1.fasta",
		"Ref/scaffold_2.fasta",
		"Ref/scaffold_3.fasta",
		"Ref/scaffold_4.fasta",
		"Ref/scaffold_5.fasta",
		"Ref/scaffold_6.fasta",
		"Ref/scaffold_7.fasta",
		"Ref/scaffold_8.fasta",
		"Ref/scaffold_9.fasta",
		"Ref/scaffold_10.fasta",
		"Ref/scaffold_11.fasta",
	log:
		"logs/chr_split.scaffold_{input.Chr_Number}.log",
	shell:
		"""
		echo Assembly_Chr_splitter.R --args -f {input} 
		ml r/4.1.0

		echo Assembly split into Chromosomes
		R --vanilla < Assembly_Chr_splitter.R --args -f {input} &&
		mv *scaffold*.fasta Ref/

		ml unload r/4.1.0
		ml unload samtools

		cp Ref/*.* EDTA_Files
		echo Splitting into chromosomes succesfuly run > {log} 2>&1
		"""

#--------------------------------------------------------------------------------
# EDTA_individual: Look for TE elements on individual Fasta Chr
#--------------------------------------------------------------------------------

rule EDTA_individual:
	input:
		rules.Chr_splitting.output,
	output:
		fa_file="EDTA_Files/scaffold_{Chrs}.fasta.mod.EDTA.TElib.fa",
		gff3_file="EDTA_Files/scaffold_{Chrs}.fasta.mod.EDTA.intact.gff3",
	log:
		"../logs/EDTA_individual_{Chrs}.log"
	shell:
		"""
		cd EDTA_Files
		eval "$(conda shell.bash hook)"

		conda activate EDTA
			echo starting EDTA process on: scaffold_{wildcards.Chrs}.fasta
			EDTA.pl --genome scaffold_{wildcards.Chrs}.fasta  2> {log}
		conda deactivate
		"""
		
#--------------------------------------------------------------------------------
# Masked_FASTA: Create masked fasta for further analysis from EDTA results.
#--------------------------------------------------------------------------------

rule Masked_FASTA:
	input:
		EDTA_repeats_file=rules.EDTA_individual.output.gff3_file,
		reference=rules.Chr_splitting.output,
	output:
		masked_fasta_file="EDTA_Files/scaffold_{Chrs}.masked.fasta",
	log:
		"../logs/Masked_FASTA{Chrs}.log"
	shell:
		"""
		ml bedtools
		cd EDTA_Files
		echo "Creating Masked Reference Genome for scaffold_{wildcards.Chrs}.masked.fasta"
		bedtools maskfasta -fi scaffold_{wildcards.Chrs}.fasta -bed scaffold_{wildcards.Chrs}.fasta.mod.EDTA.intact.gff3 -fo scaffold_{wildcards.Chrs}.masked.fasta 2> {log}
		ml unload bedtools	
		"""
#--------------------------------------------------------------------------------
# MINIMAP2: Align the Nanopore raw data to CHR bams for further processing.
#--------------------------------------------------------------------------------

rule MINIMAP2:
	input:
		Nanopore_File=NANOPORE_FASTQ,
		reference=rules.Init.output,
	output:
		Bam_file="Minimap_Aligned/minimap2.sorted.MAPQ20.dedup.bam",
		Bai_file="Minimap_Aligned/minimap2.sorted.MAPQ20.dedup.bam.bai",
	log:
		"logs/MINIMAP2.log"
	shell:
		"""
		ml minimap2
		echo Minimap will proceed with the alignment.
		echo Minimap Indexing
		minimap2 -ax map-ont -t 12 -2 {input.reference} {input.Nanopore_File} > Minimap_Aligned/minimap2.sam 
		ml unload minimap2

		ml samtools
		echo SAM to BAM convertion
		samtools view -S -b Minimap_Aligned/minimap2.sam > Minimap_Aligned/minimap2.bam

		echo  BAM sorting
		samtools sort Minimap_Aligned/minimap2.bam -o Minimap_Aligned/minimap2.sorted.bam

		echo BAM indexing.....
		samtools index  Minimap_Aligned/minimap2.sorted.bam
		​
		echo Applying MAPQ20 quality filter
		samtools view -bq 20 Minimap_Aligned/minimap2.sorted.bam > Minimap_Aligned/minimap2.sorted.MAPQ20.bam

		echo Removing duplicates from Bam File
		samtools rmdup -s Minimap_Aligned/minimap2.sorted.MAPQ20.bam  Minimap_Aligned/minimap2.sorted.MAPQ20.dedup.bam

		echo BAM indexing SORTED QUALITY FILTERED AND DEDUP.....
		samtools index  Minimap_Aligned/minimap2.sorted.MAPQ20.dedup.bam

		samtools idxstats Minimap_Aligned/minimap2.sorted.MAPQ20.dedup.bam > {log}
		"""	

#--------------------------------------------------------------------------------
# Nanopore_Raw_to_Chr: Split the raw Nanopore data into chromosomes for analysis
#--------------------------------------------------------------------------------

rule Nanopore_Raw_to_Chr:
	input:
		rules.MINIMAP2.output.Bam_file,
	output:
		"Sorted_Chromosomes/Nanopore.chr{Chrs}.fasta",
	log:
		"logs/Sorted_Chromosomes.{Chrs}.log",
	shell:
		"""
		ml samtools
		ml seqtk
		
		echo Extracting chromosome: {wildcards.Chrs}
		samtools view -b Minimap_Aligned/minimap2.sorted.MAPQ20.dedup.bam scaffold_{wildcards.Chrs} > Minimap_Aligned/chr{wildcards.Chrs}.bam
		samtools index Minimap_Aligned/chr{wildcards.Chrs}.bam
		
		echo Transforming Bam to Fasta
		samtools bam2fq Minimap_Aligned/chr{wildcards.Chrs}.bam | seqtk seq -A - > Sorted_Chromosomes/Nanopore.chr{wildcards.Chrs}.fasta 
		
		ml unload seqtk
		ml unload samtools
		"""
#--------------------------------------------------------------------------------
# Minor_Scaffolds_Correction: Add all the fasta lines to Chr11
#--------------------------------------------------------------------------------

rule Minor_Scaffolds_correction:
	input:
		"Sorted_Chromosomes/Nanopore.chr11.fasta",
	output:
		"Sorted_Chromosomes/Minor_Scaffolds.fasta",
	log:
		"logs/Minor_Scaffolds_correction.log",
	shell:
		"""
		ml samtools
		ml seqtk
		
		samtools idxstats Minimap_Aligned/minimap2.sorted.MAPQ20.dedup.bam
		Contig_Number=$( samtools idxstats Minimap_Aligned/minimap2.sorted.MAPQ20.dedup.bam | wc -l )
		Contig_Number=$((Contig_Number - 1))
		echo NUMBER OF CONTIGS: $Contig_Number 2> {log}

		for ((i=12; i<=$Contig_Number; i++))
			do
			samtools view -b Minimap_Aligned/minimap2.sorted.MAPQ20.dedup.bam scaffold_$i > Minimap_Aligned/chr$i.bam
			samtools index Minimap_Aligned/.chr$i.bam
			echo Adding to Chr11 Fasta
			samtools bam2fq Minimap_Aligned/chr$i.bam | seqtk seq -A - >> Sorted_Chromosomes/Nanopore.chr11.fasta
			done
		cp Sorted_Chromosomes/Nanopore.chr11.fasta Sorted_Chromosomes/Minor_Scaffolds.fasta
		ml unload seqtk
		ml unload samtools
		"""				

#--------------------------------------------------------------------------------
# MAKER3: Perform Maker Analysis, using 3 iterations.
#--------------------------------------------------------------------------------

rule MAKER3:
	input:
		EST_File=rules.Nanopore_Raw_to_Chr.output,
		reference=rules.Chr_splitting.output,
		Protein_File="Maker_Files/csa.trans.Protein.10072011.fasta",
		Repeats_File="Maker_Files/PK_Repeat_contigs_RE_filtered_min500bp.fa",
	output:
		"Maker_Files/Chr{Chrs}/scaffold_{Chrs}.maker.output/scaffold_{Chrs}_master_datastore_index.log",
	shell:
		"""
		BASEDIR=$PWD 
		cd Maker_Files
		
		echo Creating Maker Files Chr: {wildcards.Chrs}
		mkdir Chr{wildcards.Chrs}
		cd Chr{wildcards.Chrs}
		maker -CTL
		cd ..
		cp maker_opts.template Chr{wildcards.Chrs}/maker_opts.ctl
		echo "########################" >> Chr{wildcards.Chrs}/maker_opts.ctl
		echo "#-----Custom Parameters (these are always required)" >> Chr{wildcards.Chrs}/maker_opts.ctl
		echo genome=$BASEDIR/Ref/scaffold_{wildcards.Chrs}.fasta >> Chr{wildcards.Chrs}/maker_opts.ctl
		echo "est=$BASEDIR/{input.EST_File}" >> Chr{wildcards.Chrs}/maker_opts.ctl
		echo "protein=$BASEDIR/{input.Protein_File}" >> Chr{wildcards.Chrs}/maker_opts.ctl
		echo "repeat_protein=$BASEDIR/{input.Repeats_File}" >> Chr{wildcards.Chrs}/maker_opts.ctl
		echo "augustus_species=arabidopsis" >> Chr{wildcards.Chrs}/maker_opts.ctl
		
		cd Chr{wildcards.Chrs}
		
		ml postgresql/13.2
		ml openmpi/4.0.3 
		mpiexec --use-hwthread-cpus maker
		ml unload postgresql/13.2
		ml openmpi/4.0.3 
		"""
#--------------------------------------------------------------------------------
# POST_MAKER: Perform Maker Analysis, using 3 differente AED filter values
#--------------------------------------------------------------------------------

rule POST_MAKER:
	input:
		Maker_File=rules.MAKER3.output,
	output:
		"Maker_Files/Chr{Chrs}/MAKER_Filtered.scaffold_{Chrs}.all.AED_{AED_filter}.gff3",
	params:
		project=PROJECT,
		AED_filter=AED_FILTER,
	shell:
		"""
		cd Maker_Files
		cp quality_filter.pl Chr{wildcards.Chrs}/scaffold_{wildcards.Chrs}.maker.output
		
		cd Chr{wildcards.Chrs}/scaffold_{wildcards.Chrs}.maker.output
		fasta_merge -d scaffold_{wildcards.Chrs}_master_datastore_index.log
		gff3_merge -d scaffold_{wildcards.Chrs}_master_datastore_index.log
		maker_map_ids --prefix {params.project}_Chr{wildcards.Chrs} --justify 8 scaffold_{wildcards.Chrs}.all.gff > map
		map_gff_ids map scaffold_{wildcards.Chrs}.all.gff
		#map_fasta_ids map scaffold_{wildcards.Chrs}.all.maker.proteins.fasta
		#map_fasta_ids map scaffold_{wildcards.Chrs}.all.maker.transcripts.fasta
		#####AED filter#####
		
		perl quality_filter.pl -a {wildcards.AED_filter} scaffold_{wildcards.Chrs}.all.gff > scaffold_{wildcards.Chrs}.all.AED_{wildcards.AED_filter}.gff 

		cd ..
		sed -n 1p scaffold_{wildcards.Chrs}.maker.output/scaffold_{wildcards.Chrs}.all.AED_{wildcards.AED_filter}.gff  > MAKER_Filtered.scaffold_{wildcards.Chrs}.all.AED_{wildcards.AED_filter}.gff3
		awk '$2 == "maker"' scaffold_{wildcards.Chrs}.maker.output/scaffold_{wildcards.Chrs}.all.AED_{wildcards.AED_filter}.gff >> MAKER_Filtered.scaffold_{wildcards.Chrs}.all.AED_{wildcards.AED_filter}.gff3
		cd ../../.
		"""

#------------------------------------------------------------------------------------
# ORF_analysis: Uses EDTA maked fasta file and ORF analysis to determine viable genes
#------------------------------------------------------------------------------------

rule ORF_analysis:
	input:
		GFF3_File=rules.POST_MAKER.output,
		Ref_File=rules.Chr_splitting.output,
		Masked_FASTA_File=rules.Masked_FASTA.output,
	output:
		"Post_Maker_Files/MAKER_ORF_Filtered_{Project}.scaffold_{Chrs}.AED_{AED_filter}.gff3",
	params:
		project=PROJECT,
		AED_filter=AED_FILTER,
	shell:
		"""
		BASEDIR=$PWD 
		cp -v Post_Maker_ORF_Analysis_Terminal.R Post_Maker_Files
		cd Post_Maker_Files 

		ml r/4.1.0
		R --vanilla < Post_Maker_ORF_Analysis_Terminal.R --args -g $BASEDIR/Maker_Files/Chr{wildcards.Chrs}/MAKER_Filtered.scaffold_{wildcards.Chrs}.all.AED_{wildcards.AED_filter}.gff3 -a $BASEDIR/Ref/scaffold_{wildcards.Chrs}.fasta -m $BASEDIR/EDTA_Files/scaffold_{wildcards.Chrs}.masked.fasta -o {params.project}.scaffold_{wildcards.Chrs}.AED_{wildcards.AED_filter}
		ml unload r/4.1.0
		cd ..
		"""

#------------------------------------------------------------------------------------
# Chr_merge: Fuse all gff3 individual chromosomes into complete assembly again
#------------------------------------------------------------------------------------

rule Chr_merge:
	input:
		expand("Post_Maker_Files/MAKER_ORF_Filtered_{Project}.scaffold_{Chrs}.AED_{AED_filter}.gff3",Project=PROJECT,Chrs = CHRS,AED_filter=AED_FILTER),
	output:
		"FINAL_ANNOTATION/FINAL_{Project}.AED_{AED_filter}.sorted.gff3",
	params:
		project=PROJECT,
		Chrs=CHRS,
	shell:
		"""
		cat Post_Maker_Files/MAKER_ORF_Filtered_{params.project}.scaffold_1.AED_{wildcards.AED_filter}.gff3 > FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.gff3
		
		for i in {{2..11}}
			do
			tail -n +4 Post_Maker_Files/MAKER_ORF_Filtered_{params.project}.scaffold_$i.AED_{wildcards.AED_filter}.gff3 >> FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.gff3
			done
		
		/home/jmlazaro/github/gff3sort/gff3sort.pl --chr_order original FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.gff3 > FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.sorted.gff3
		"""
#------------------------------------------------------------------------------------
# COGNATE: Generate Protein Fasta file and Analysis of gff3
#------------------------------------------------------------------------------------

rule COGNATE:
	input:
		GFF_file=rules.Chr_merge.output,
		Ref_file=rules.Init.output
	output:
		"COGNATE/COGNATE_{Project}.AED_{AED_filter}/COGNATE_{Project}.AED_{AED_filter}_00-analyzed_transcripts.fa",
		"COGNATE/COGNATE_{Project}.AED_{AED_filter}/COGNATE_{Project}.AED_{AED_filter}_01-summary.tsv",
	params:
		project=PROJECT,
	shell:
		"""
		BASEDIR=$PWD
		cp FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.sorted.gff3 COGNATE
		
		ml perl
		cd COGNATE
		WORKDIR=$PWD
		perl /home/jmlazaro/github/COGNATE/COGNATE_v1.0/COGNATE_v1.0.pl --gff $BASEDIR/COGNATE/FINAL_{params.project}.AED_{wildcards.AED_filter}.sorted.gff3 --fasta $BASEDIR/{input.Ref_file} --name {params.project}.AED_{wildcards.AED_filter} --workingdir $WORKDIR --overwrite
		ml unload perl
		"""


		
#perl /home/jmlazaro/github/COGNATE/COGNATE_v1.0/COGNATE_v1.0.pl --gff /home/jmlazaro/scratch/Sunflower5/COGNATE/HAN412_Eugene_curated_v1_1.gff3 --fasta /home/jmlazaro/scratch/Sunflower5/COGNATE/Ha412HOv2.0-20181130.fasta --name HAN412 	
#(BUSCO) [jose@mrfox XRQv2_BUSCO]$ busco -m genome -i XRQv2.genes.fa -o V2_busco --auto-lineage