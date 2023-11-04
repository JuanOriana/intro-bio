import sys

# pip install primer3-py
import primer3
from Bio import SeqIO
from Bio.Seq import Seq

from ej1 import read_genbank_file


def get_primers(sequence, num_primers=5, min_len=18, max_len=24, min_gc=50, max_gc=60, max_tm=67):
    design = primer3.bindings.design_primers(
        {
            'SEQUENCE_ID': 'example',
            'SEQUENCE_TEMPLATE': sequence,
        },
        {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_MIN_SIZE': min_len,
            'PRIMER_MAX_SIZE': max_len,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 50.0,
            'PRIMER_MAX_TM': max_tm,
            'PRIMER_MIN_GC': min_gc,
            'PRIMER_MAX_GC': max_gc,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8
        }
    )

    primers = [Seq(design['PRIMER_LEFT_{}_SEQUENCE'.format(i)]) for i in range(num_primers)]
    return primers


def main():
    if len(sys.argv) != 2:
        print("Wrong arguments amount")
        exit(1)

    faa_filename = sys.argv[1]
    in_file = open(faa_filename)
    records = SeqIO.parse(in_file, format="fasta")
    for record in records:
        primers = get_primers(record.seq)
        print(primers)


if __name__ == '__main__':
    main()
