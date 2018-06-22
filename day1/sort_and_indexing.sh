#!/bin/bash
samtools sort star/example1/example1Aligned.out.bam > star/example1/example1Aligned.sort.bam
samtools index star/example1/example1Aligned.sort.bam
ls -l star/example
