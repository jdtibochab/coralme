#!/bin/bash

mkdir -p $1/blast
mv $1/org.faa $1/blast/
mv $1/ref.faa $1/blast/
cd $1/blast/

echo "\n################# Creating DB files #################\n"
makeblastdb -in org.faa -dbtype prot # -title {org name?} -out {prefix}
makeblastdb -in ref.faa -dbtype prot # -title {org name?} -out {prefix}

echo "\n################# Blasting #################\n"

blastp -db org.faa -query ref.faa -num_threads 4 -out org_as_db.txt -outfmt 6
blastp -db ref.faa -query org.faa -num_threads 4 -out ref_as_db.txt -outfmt 6

echo "\n################# Processing output files #################\n"
sed -i 's/>//g' org_as_db.txt
sed -i 's/>//g' ref_as_db.txt

BLAST_COLUMNS="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n"
sed -i "1s/^/$BLAST_COLUMNS/" org_as_db.txt
sed -i "1s/^/$BLAST_COLUMNS/" ref_as_db.txt
