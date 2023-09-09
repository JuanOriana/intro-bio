from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os


def read_genbank_file(filename):
    with open(filename, "r") as file:
        records = SeqIO.parse(file, "genbank")
        return list(records)


def translate_mrna_to_aminoacids(to_trans_mrna_sequence, min_seq_len=90, start_codon="ATG"):
    translations = []

    for mrna_sequence in to_trans_mrna_sequence, to_trans_mrna_sequence.complement():
        for i in range(len(mrna_sequence) - 2):
            if mrna_sequence[i] + mrna_sequence[i + 1] + mrna_sequence[i + 2] == start_codon:
                max_seq_length = len(mrna_sequence) - i
                max_seq_length_rounded = max_seq_length - (max_seq_length % 3)
                translation = mrna_sequence[i:max_seq_length_rounded + i].translate(to_stop=True)
                if len(translation) > min_seq_len:
                    translations.append((len(translation), translation, i))

    translations.sort(reverse=True)
    return translations


def write_fasta_file(output_filename, translation, record_data):
    with open(output_filename, "w") as file:
        record = SeqRecord(translation, id=record_data.name, description=record_data.description)
        SeqIO.write(record, file, "fasta")


def translate_genbank_to_fasta(records):
    translations = {}
    for record in records:
        translations[record.id] = translate_mrna_to_aminoacids(record.seq)
    return translations


def main():
    if len(sys.argv) != 2:
        print("Wrong arguments amount")
        exit(1)

    gbk_filename = sys.argv[1]

    if gbk_filename.split(".")[-1] != "gbk" and gbk_filename.split(".")[-1] != "gb":
        print("Wrong format for input file")
        exit(1)

    faa_filename = "./FASTA/" + gbk_filename.split(".")[0].split("/")[-1] + ".fasta"
    print(faa_filename)
    os.makedirs(os.path.dirname(faa_filename), exist_ok=True)

    genbank_records = []
    try:
        genbank_records = read_genbank_file(gbk_filename)
    except IOError:
        print("Could not open " + gbk_filename)

    translations = translate_genbank_to_fasta(genbank_records)

    for record in genbank_records:
        # TODO: if a .gb file contains many records, make a fasta file for each
        # faa_filename = "./FASTA/" + record.name + ".fasta"
        # print(faa_filename)

        print(">>> RECORD " + record.id)
        print("-----------------")
        for translation in translations[record.id]:
            print("sequence: ", translation[1], "length: ", translation[0], " start position: ", translation[2])
            print("-----------------")

        try:
            write_fasta_file(faa_filename, translations[record.id][0][1], record)
        except IOError:
            print("Could not open " + gbk_filename)

        print(
            f"Las secuencias de aminoácidos del {record.name} ({record.description}) "
            f"se han guardado en {faa_filename}.")


if __name__ == '__main__':
    main()
