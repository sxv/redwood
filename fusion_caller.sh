#!/bin/bash

WORKING_DIR=/media/lab/drive3/Stanford_LMS_fusions
ALIGN=/media/lab/drive2/RNA_known_fusions
CHIMERASCAN=/media/lab/drive2/chimerascan-0.4.5
DEFUSE=/media/lab/drive2/defuse-0.6.2
FUSIONCATCHER=/media/lab/drive3/fusioncatcher

for path in ''
do
        casename=$(basename $path)
        cd $path
        unpigz -p 40 ./*fastq.gz
        for file1 in ./*1.fastq; do mv $file1 ${file1/1.fastq/1.end1.fastq}; done
        for file2 in ./*2.fastq; do mv $file2 ${file2/2.fastq/2.end2.fastq}; done
        fastq1=*end1.fastq
        fastq2=*end2.fastq
        bash $WORKING_DIR/align.sh $fastq1 $fastq2
        bash $DEFUSE/defuse.sh $fastq1 $fastq2
        bash $CHIMERASCAN/chimerascan.sh $fastq1 $fastq2
        bash $FUSIONCATCHER/fusion_catcher.sh ./
#       bash $FUSIONCATCHER/fusion_catcher_paranoid.sh ./
        ls -lhtr | awk '{print $6,$7,$8,$9,$10}' > timestamp.txt
        cp ./defuse_out/results.tsv ./${casename}.fumble.results.tsv
        cp ./chimerascan_out/chimeras.bedpe ./${casename}.fumble.chimeras.bedpe
        cp ./*catcher*/final-list_candidate-fusion-genes.GRCh37.txt ./${casename}.fumble.final-list_candidate-fusion-genes.GRCh37.txt
        cp ./star-fusion.fusion_candidates.txt ./${casename}.fumble.star-fusion.fusion_candidates.txt
        cd ..
done
