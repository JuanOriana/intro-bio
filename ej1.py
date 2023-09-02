from Bio import SeqIO 
import sys

if len(sys.argv) != 3:
    print("Wrong arguments amount")
    exit(1)

gbk_filename = sys.argv[1]
faa_filename = sys.argv[2]

line_length = 70


if (gbk_filename.split(".")[1] != "gbk" and gbk_filename.split(".")[1] != "gb") or faa_filename.split(".")[1] != "fasta":
    print("Wrong format for input or output file")
    exit(1)

input_handle = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank"):
    print("Dealing with GenBank record %s" % seq_record.id)
    aux = str(seq_record.seq)
    seq = ''.join([aux[i:i+line_length]+"\n" for i in range(0, len(aux), line_length)])

    output_handle.write(">%s %s\n%s\n" % (
            seq_record.id,
            seq_record.description,
            seq))

output_handle.close()
input_handle.close()