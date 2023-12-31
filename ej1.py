from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os


def read_genbank_file(filename):
    with open(filename, "r") as file:
        records = SeqIO.parse(file, "genbank")
        return list(records)


def translate_mrna_to_aminoacids(to_trans_mrna_sequence, min_seq_len=75, start_codon="ATG"):
    translations = []
    nucleotides = []

    is_complement = False
    for mrna_sequence in to_trans_mrna_sequence, to_trans_mrna_sequence.reverse_complement():
        for i in range(len(mrna_sequence) - 2):
            if mrna_sequence[i] + mrna_sequence[i + 1] + mrna_sequence[i + 2] == start_codon:
                max_seq_length = len(mrna_sequence) - i
                max_seq_length_rounded = max_seq_length - (max_seq_length % 3)
                codon_sequence = mrna_sequence[i:i + max_seq_length_rounded]
                translation = codon_sequence.translate(to_stop=True)
                if len(translation) > min_seq_len:
                    if is_complement:
                        translations.append((len(translation), translation, len(mrna_sequence) - (i + 1),
                                             len(mrna_sequence) - ((i + 1) + len(translation) * 3)))
                    else:
                        translations.append((len(translation), translation, (i + 1), (i + 1) + len(translation) * 3))
                nucleotides.append((len(translation), codon_sequence[: len(translation) * 3 + 3]))
        is_complement = True

    translations.sort(reverse=True)
    nucleotides.sort(reverse=True)

    return translations, nucleotides


def write_fasta_file(output_filename, translation, record_data):
    with open(output_filename, "w") as file:
        record = SeqRecord(translation, id=record_data.name, description=record_data.description)
        SeqIO.write(record, file, "fasta")


def translate_genbank_to_fasta(records):
    translations = {}
    nucleotides = {}
    for record in records:
        translations[record.id], nucleotides[record.id] = translate_mrna_to_aminoacids(record.seq)
    return translations, nucleotides


def main():
    if len(sys.argv) != 2:
        print("Wrong arguments amount")
        exit(1)

    gbk_filename = sys.argv[1]

    if gbk_filename.split(".")[-1] != "gbk" and gbk_filename.split(".")[-1] != "gb":
        print("Wrong format for input file")
        exit(1)

    faa_filename = "./FASTA/" + gbk_filename.split(".")[0].split("/")[-1] + ".fasta"
    nfaa_filename = "./FASTA/NUCLEOTIDE-" + gbk_filename.split(".")[0].split("/")[-1] + ".fasta"

    print(faa_filename)
    print(nfaa_filename)
    os.makedirs(os.path.dirname(faa_filename), exist_ok=True)

    genbank_records = []
    try:
        genbank_records = read_genbank_file(gbk_filename)
    except IOError:
        print("Could not open " + gbk_filename)

    translations, nucleotides = translate_genbank_to_fasta(genbank_records)

    for record in genbank_records:
        # TODO: if a .gb file contains many records, make a fasta file for each
        # faa_filename = "./FASTA/" + record.name + ".fasta"
        # print(faa_filename)

        print(">>> RECORD " + record.id)
        print("-----------------")
        for translation in translations[record.id]:
            print("sequence: ", translation[1], "\nlength: ", translation[0], ", start position: ", translation[2],
                  ", end position: ", translation[3])
            print("-----------------")

        try:
            write_fasta_file(faa_filename, translations[record.id][0][1], record)
            write_fasta_file(nfaa_filename, nucleotides[record.id][0][1], record)
        except IOError:
            print("Could not open " + gbk_filename)

        print(
            f"Las secuencias de aminoácidos del {record.name} ({record.description}) "
            f"se han guardado en {faa_filename}.")


if __name__ == '__main__':
    main()
