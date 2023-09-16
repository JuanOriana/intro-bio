import sys, os
from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline


def main():
    if len(sys.argv) != 2:
        print("Wrong amount of arguments, must be 1")
        exit(1)

    input_file = sys.argv[1]

    if input_file.split(".")[-1] != "fasta":
        print("Wrong format for input or output file")
        exit(1)

    file_name = os.path.splitext(os.path.basename(input_file))[0]
    output_file = f"./FASTA/{file_name}_MSA.fasta"

    mafft_cline = MafftCommandline(input=input_file)
    stdout, stderr = mafft_cline()
    with open(output_file, "w") as output_handle:
        output_handle.write(stdout)

    alignment = AlignIO.read(output_file, "fasta")
    print(alignment)
    # print("-----------------")

    # for seq in alignment:
    #     print(seq.seq)

    # summary = AlignInfo.SummaryInfo(alignment)
    # identity = summary.pos_specific_score_matrix()
    # print(identity)


if __name__ == '__main__':
    main()
