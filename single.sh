#!/bin/bash

#PBS -l nodes=1:ppn=10
#PBS -l gres=localhd:20
#PBS -l mem=20gb
#PBS -l vmem=20gb
#PBS -l walltime=4:00:00

module load bwa/0.7.8
module load samtools/1.5
module load python/3.7.0

dir=/hpf/projects/nsondheimer/mitoseq/acelik_fixes/parse_pileup
ref=/hpf/projects/nsondheimer/mitoseq/acelik_fixes/parse_pileup/ref/drCRS.fa

mkdir -p $outdir/alignment

if [ $fastq2 = "None" ]
then
	echo "Single end reads"
	bwa mem -t 10 $ref \
  	<(gunzip -d -c $fastq1) \
 	| samtools view -@10 -b - | samtools sort -@ 10 - > $outdir/alignment/$sample.bam 2>$outdir/$sample.align.log
else
	echo "Paired end reads"
	bwa mem -t 10 $ref \
    <(gunzip -d -c $fastq1) \
    <(gunzip -d -c $fastq2) \
  | samtools view -@ 10 -b - | samtools sort -@ 10 - > $outdir/alignment/$sample.bam 2> $outdir/$sample.align.log
fi

mkdir -p $outdir/pileup

samtools mpileup -d 500000 -B -f $ref $outdir/alignment/$sample.bam > $outdir/pileup/$sample.pileup 2> $outdir/$sample.pileup.log

mkdir -p $outdir/processed

python3 $dir/parse_pileup.py -m $outdir/pileup/$sample.pileup \
        -o $outdir/processed/$sample.processed.csv -r $ref -f 2> $outdir/$sample.process.log

