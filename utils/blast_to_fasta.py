from Bio.Blast import NCBIXML
from Bio import Entrez, SeqIO
import sys
import shutil


def blast_to_fasta(blast_results, top_n):
    accession_ids = []
    for blast_record in blast_results:
        for alignment in blast_record.alignments:
            accession_id = alignment.accession
            accession_ids.append(accession_id)
            if len(accession_ids) >= top_n:
                break
        if len(accession_ids) >= top_n:
            break

    sequences = []
    try:
        handle = Entrez.efetch(db="protein", id=accession_ids, rettype="fasta", retmode="text")
        batch_sequences = handle.read()
        handle.close()
        sequences.append(batch_sequences)

    except Exception as e:
        print(f"Failed to download sequences for batch: {str(e)}")

    print(f"Downloaded and saved {len(accession_ids)} sequences.")
    print("-----------------")
    return sequences


def main():
    if len(sys.argv) <= 3:
        print("Wrong amount of arguments, must be at least 3")
        exit(1)

    faa_filename = sys.argv[1]
    blast_filename = sys.argv[2]
    output_filename = sys.argv[3]
    Entrez.email = "pdomingues@itba.edu.ar"
    top_n = int(sys.argv[4]) if len(sys.argv) >= 5 else 10

    if faa_filename.split(".")[-1] != "fasta":
        print("Wrong format for input or output file")
        exit(1)

    if blast_filename.split(".")[-1] != "xml":
        print("Wrong format for input or output file")
        sys.exit(1)

    shutil.copy(faa_filename, output_filename)

    blast_results = NCBIXML.parse(open(blast_filename))
    sequences = blast_to_fasta(blast_results, top_n)

    with open(output_filename, "a") as seq_file:
        seq_file.write('\n')
        for sequence in sequences:
            seq_file.write(sequence)


if __name__ == '__main__':
    main()
