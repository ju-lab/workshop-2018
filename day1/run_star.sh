mkdir -p star/example1
STAR --runMode alignReads \
	--outSAMtype BAM Unsorted \
	--genomeDir /home/users/kjyi/ref/hg19/star_index \
	--outFileNamePrefix star/example1/example1 \
	--chimSegmentMin 20 \
	--runThreadN 20 \
	--chimOutType WithinBAM \
	--readFilesCommand zcat \
	--readFilesIn fastq/example_R1.fastq.gz fastq/example_R2.fastq.gz
