#!/bin/bash 

folder=$1
index=$2
donor=$3
recipient=$4
nb_line=0
ls $folder > fastq.txt
list_fast=()
list_sam=()
tmp="def"
bwa index -p $index $index
while read line; do 
    list_fast+=($line)
    if [ $(expr $nb_line % 2) == "1" ];
        then
        out=$(awk -F'_S' '{print $1}' <<< ${list_fast[0]})
        if [ "$tmp" != "$out" ];
            then
                tmp=$out
                list_sam+=($out)
        fi
        bwa mem $index $folder/${list_fast[0]} $1/${list_fast[1]} >> sam/$out.sam
        list_fast=()
        nb_line=0
        else
        let "nb_line+=1"
    fi
done <fastq.txt

for sam in ${list_sam[@]}
do
    python ecoli.py $sam $recipient $donor
done

rm fastq.txt
rm sam/*.sam
rm index/*.fasta.*
