#!/usr/bin/env bash
counter=1

for file in control/*.fastq
do
    filename=$(basename "${file%_1.fastq}")

    if [[ $filename == "SRR16013071" ]]; then
        echo "$filename"

       # fastqc control/"$filename"_1.fastq control/"$filename"_2.fastq -o fastqc/control
       # echo "fastqc done"

        trimmomatic PE -threads 10 -phred33 control/"$filename"_1.fastq control/"$filename"_2.fastq trimctrl/"$filename"_1_trimmed.fastq trimctrl/"$filename"_1_unpaired.fastq trimctrl/"$filename"_2_trimmed.fastq trimctrl/"$filename"_2_unpaired.fastq TRAILING:10
        echo "trimmomatic done"

        hisat2 -q -x genome/myindex -1 trimctrl/"$filename"_1_trimmed.fastq -2 trimctrl/"$filename"_2_trimmed.fastq | samtools sort -o bam/control/"$filename".bam 
        echo "hisat$counter done"

        featureCounts -a genome/Homo_sapiens.GRCh38.106.gtf -o counts/control/"$filename"_featurecounts.txt bam/control/"$filename".bam 
        echo "done$counter"
        
        ((counter++))
    fi
done
