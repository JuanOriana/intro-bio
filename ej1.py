from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys


def read_genbank_file(filename):
    with open(filename, "r") as file:
        records = SeqIO.parse(file, "genbank")
        return list(records)


def translate_mrna_to_aminoacids(mrna_sequence):
    start_codon = "ATG"  # deberia ser AUG pero es ATG porque tenemos adn no arn
    translations = []

    while start_codon in mrna_sequence:
        start_position = mrna_sequence.find(start_codon)
        sequence_length = len(mrna_sequence) - start_position
        new_length = sequence_length - (sequence_length % 3)
        to_translate = mrna_sequence[start_position:new_length + start_position]

        translation = to_translate.translate(to_stop=True)
        mrna_sequence = mrna_sequence[start_position + len(start_codon):]
        translations.append(translation)
    return translations


def write_fasta_file(output_filename, translation, record_data):
    with open(output_filename, "w") as file:
        record = SeqRecord(translation, id=record_data.name, description=record_data.description)
        SeqIO.write(record, file, "fasta")


def translate_genbank_to_fasta(input_filename, output_filename):
    genbank_records = read_genbank_file(input_filename)

    for record in genbank_records:
        translations = translate_mrna_to_aminoacids(record.seq)
        longest_item = max(translations, key=len)
        write_fasta_file(output_filename, Seq(longest_item), record)
        print(
            f"Las secuencias de aminoácidos del {record.name}({record.description}) "
            f"se han guardado en {output_filename}.")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Wrong arguments amount")
        exit(1)

    gbk_filename = sys.argv[1]

    line_length = 70

    if gbk_filename.split(".")[1] != "gbk" and gbk_filename.split(".")[1] != "gb":
        print("Wrong format for input or output file")
        exit(1)

    faa_filename = "./FASTA/" + gbk_filename.split(".")[0].split("/")[-1] + ".fasta"
    print(faa_filename)

    translate_genbank_to_fasta(gbk_filename, faa_filename)
