#--------------------------------------------------------------------------------
# MAIN CONFIGURATION PARAMETERS (INTERNAL CONFIG FILE)
#--------------------------------------------------------------------------------

CHRS = '1 2 3 4 5 6 7 8 9 10 11'.split()
PROJECT = "PK_hap2"
REFERENCE = "PK_hap2.09082021.fasta"
NANOPORE_FASTQ = "FAP02733_Porechop_No_Trim_Final.fastq.gz"

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
		expand("Maker_Files/chr{Chrs}/scaffold_{Chrs}_master_datastore_index.log",Chrs = CHRS),

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
# MAKER3: Split the raw Nanopore data into chromosomes for analysis
#--------------------------------------------------------------------------------

rule MAKER3:
	input:
		EST_File=rules.Nanopore_Raw_to_Chr.output,
		reference=rules.Chr_splitting.output,
		Protein_File="Maker_Files/csa.trans.Protein.10072011.fasta",
		Repeats_File="Maker_Files/PK_Repeat_contigs_RE_filtered_min500bp.fa",
	output:
		"Maker_Files/chr{Chrs}/scaffold_{Chrs}_master_datastore_index.log",
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
		
