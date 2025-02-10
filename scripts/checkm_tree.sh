#! /bin/bash

mkdir $1/drep/mags_faa
for file in $1/drep/dereplicated_genomes/*
do
    prodigal -i $1/drep/dereplicated_genomes/"$(basename "${file}")" -a $1/drep/mags_faa/"$(basename "${file%.fa}")".faa -p meta
done

checkm tree $1/drep/mags_faa -t 10 -g -x faa $1/phylo_tree

checkm tree_qa $1/phylo_tree -o 5 -f $1/phylo_tree/checkm_alignment