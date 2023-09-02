from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys


def read_genbank_file(filename):
    with open(filename, "r") as file:
        records = SeqIO.parse(file, "genbank")
        return list(records)


def translate_mrna_to_aminoacids(mrna_sequence, start_position):
    sequence_length = len(mrna_sequence) - start_position
    new_length = sequence_length - (sequence_length % 3) + start_position
    translations = []
    mrna_sequence = mrna_sequence[start_position:new_length]
    translations.append(mrna_sequence.translate(to_stop=True))
    return translations


def write_fasta_file(output_filename, translations, record_data):
    with open(output_filename, "w") as file:
        for i, translation in enumerate(translations):
            record = SeqRecord(translation, id=record_data.name, description=record_data.description)
            SeqIO.write(record, file, "fasta")


def find_start_codon_position(mrna_sequence):
    start_codon = "ATG"  # deberia ser AUG pero es ATG porque tenemos adn no arn
    return mrna_sequence.find(start_codon)


def translate_genbank_to_fasta(input_filename, output_filename):
    genbank_records = read_genbank_file(input_filename)

    for idx, record in enumerate(genbank_records):
        start_position = find_start_codon_position(record.seq)
        if start_position != -1:
            translations = translate_mrna_to_aminoacids(record.seq, start_position)
            write_fasta_file(output_filename, translations, record)
            # print(record)
            print(
                f"Las secuencias de aminoácidos del {record.name}({record.description}) se han guardado en {output_filename}")
            return
        else:
            print(f"No se encontró el codón de inicio 'AUG' en la secuencia mRNA {record.description}.")


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
    # input_handle = open(gbk_filename, "r")
    # output_handle = open(faa_filename, "w")
    #
    # for seq_record in SeqIO.parse(input_handle, "genbank"):
    #     print("Dealing with GenBank record %s" % seq_record.id)
    #     aux = str(seq_record.seq)
    #     seq = ''.join([aux[i:i + line_length] + "\n" for i in range(0, len(aux), line_length)])
    #
    #     output_handle.write(">%s %s\n%s\n" % (
    #         seq_record.id,
    #         seq_record.description,
    #         seq))
    #
    # output_handle.close()
    # input_handle.close()
