#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --job-name=fetch_refs
#SBATCH --error=logs/fetch_refs.%J.err
#SBATCH --output=logs/fetch_refs.%J.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200GB
#SBATCH --partition=guest

module purge
module load entrez-direct 
module load seqkit 

cd /work/richlab/aliciarich/marmoset_microbiome

uid_file="dataframes/fetch_references.txt"
accessions_file="dataframes/accessions.txt"
fasta_file="refs/abund_refseqs_new.fasta"
renamed_fasta="refs/abund_refseqs_new_renamed.fasta"

> "$fasta_file"

cut -f1 "$uid_file" > "$accessions_file"

if efetch -db nuccore -input "$accessions_file" -format fasta -email "aliciarich@unomaha.edu" > "$fasta_file"; then
    echo "Sequences fetched successfully."
else
    echo "Error fetching sequences." >&2
    exit 1
fi

seqkit replace $fasta_file -p '^(\S+)(.+?)$' -r '{kv}$2' -k $uid_file -o $renamed_fasta

rm $fasta_file
mv $renamed_fasta $fasta_file
