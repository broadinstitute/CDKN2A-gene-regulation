### STARsolo commands used to align the data in the screen
###
# Generate index for mRNA alignment
STAR --runMode genomeGenerate --genomeDir . --genomeFastaFiles ../reference/GRCh38.p13.genome.fa --sjdbGTFfile ../reference/gencode.v34.annotation.gtf

# Generate index for CROP/dialup alignment
# In this FASTA the guide sequences + CROP are added as pseudo-chromosomes
# The files CROP_sgRNA.fa and CROP_sgRNA.gtf are in the Google bucket
# These were created using this Python code: https://github.com/powellgenomicslab/CROP-seq
# (modifed for our sequences)
# The below fa is the above CROP_sgRNA.fa concatenated to the standard references above
# The gtf is this one: gencode.v34.GRCh38.annotation.patched_p16_locus.v2.gtf and the CROP gtf concatenated
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir . --genomeFastaFiles GRCh38.CROP.fa --sjdbGTFfile GRCh38.CROP.gtf

# Align mRNA (per array)
STAR --runThreadN 16 --genomeDir ../reference_original --readFilesCommand zcat --readFilesIn mRNA_array3_S7_R2_001.fastq.gz mRNA_array3_S7_R1_001.fastq.gz --outFileNamePrefix mRNA_array3_orig_ref_ --outSAMtype BAM SortedByCoordinate --soloType CB_UMI_Simple --soloCBwhitelist None --soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 --quantMode TranscriptomeSAM GeneCounts --alignSoftClipAtReferenceEnds Yes --outSAMstrandField intronMotif --outFilterMismatchNoverLmax 0.1 --chimMainSegmentMultNmax 1 --outFilterScoreMin 30 --outSAMattributes CR UR CY UY CB UB --outFilterMatchNminOverLread 0.33 --outSAMattrRGline ID:rg1 SM:sm1 --chimJunctionOverhangMin 15 --alignIntronMax 1000000 --chimOutType Junctions WithinBAM SoftClip --outFilterMultimapNmax 20 --limitSjdbInsertNsj 1200000 --alignSJoverhangMin 8 --outFilterType BySJout --outFilterScoreMinOverLread 0.33 --outFilterIntronMotifs None --chimSegmentMin 15 --alignIntronMin 20 --alignSJDBoverhangMin 1 --alignMatesGapMax 999999 --outFilterMismatchNmax 999

# Align dialup (per array)
STAR --runThreadN 16 --genomeDir ../reference_CROP_stitch/ --readFilesCommand zcat --readFilesIn Dialup-CROP-array3_R2_001.fastq.gz,Dialup-p16-array3_R2_001.fastq.gz Dialup-CROP-array3_R1_001.fastq.gz,Dialup-p16-array3_R1_001.fastq.gz --outFileNamePrefix Dialup-CROP-p16-array3 --outSAMtype BAM SortedByCoordinate --soloType CB_UMI_Simple --soloCBwhitelist None --soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 --quantMode TranscriptomeSAM GeneCounts --alignSoftClipAtReferenceEnds Yes --outSAMstrandField intronMotif --outFilterMismatchNoverLmax 0.1 --chimMainSegmentMultNmax 1 --outFilterScoreMin 30 --outSAMattributes CR UR CY UY CB UB --outFilterMatchNminOverLread 0 --outSAMattrRGline ID:rg1 SM:sm1 --chimJunctionOverhangMin 15 --alignIntronMax 1000000 --chimOutType Junctions WithinBAM SoftClip --outFilterMultimapNmax 20 --limitSjdbInsertNsj 1200000 --alignSJoverhangMin 8 --outFilterType BySJout --outFilterScoreMinOverLread 0 --outFilterIntronMotifs None --chimSegmentMin 15 --alignIntronMin 20 --alignSJDBoverhangMin 1 --alignMatesGapMax 999999 --outFilterMismatchNmax 999 --outFilterMatchNmin 0

