mkdir -p rsem
/usr/local/bin/rsem-calculate-expression \
	--num-threads 2 \
	--no-bam-output \
	--estimate-rspd \
	--bam ./star/NSCLC_01_NTIL.Aligned.toTranscriptome.out.bam \
	/home/users/kjyi/ref/hg19/rsem_reference/rsem_reference \
	rsem/NSCLC_01_NTIL
