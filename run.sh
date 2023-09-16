#!/bin/bash

gene=CNGA3-NM_001298
sequences_amount=10
ej3_input=${gene}_BLAST

python ej1.py GenBank/${gene}.gb
# python ej2.py FASTA/${gene}.fasta
python utils/blast_to_fasta.py FASTA/${gene}.fasta BLAST/${gene}.xml FASTA/${ej3_input}.fasta ${sequences_amount}
#if ! command -v mafft &>/dev/null; then
#    echo "MAFFT is not installed. Installing..."
#    brew install mafft
#fi && \
python ej3.py FASTA/${ej3_input}.fasta
