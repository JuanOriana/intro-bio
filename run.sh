#!/bin/bash

gene=CNGA3-NM_001298
sequences_amount=25
ej3_input=${gene}_BLAST
ej4_input=${gene}_FASTA

python3 ej1.py GenBank/${gene}.gb
python3 ej2.py FASTA/${gene}.fasta
python3 utils/blast_to_fasta.py FASTA/${gene}.fasta BLAST/${gene}.xml FASTA/${ej3_input}.fasta ${sequences_amount}
if ! command -v mafft &>/dev/null; then
    echo "MAFFT is not installed. Installing..."
    brew install mafft
fi && \
python3 ej3.py FASTA/${ej3_input}.fasta
python3 ej4.py inputs/${ej4_input}.fasta  output.patmatmotifs
#python3 ej5.py FASTA/NUCLEOTIDE-${gene}.fasta
