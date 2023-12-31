import sys
import os
from Bio.Blast import NCBIWWW
from Bio import SeqIO

if len(sys.argv) != 2:
    print("Wrong amount of arguments, must be 1")
    exit(1)

faa_filename = sys.argv[1]

if faa_filename.split(".")[-1] != "fasta":
    print("Wrong format for input or output file")
    exit(1)

blast_filename = "./BLAST/" + faa_filename.split(".")[0].split("/")[-1] + ".xml"
os.makedirs(os.path.dirname(blast_filename), exist_ok=True)

in_file = open(faa_filename)
record = SeqIO.parse(in_file, format="fasta")
save_file = open(blast_filename, "w+")

for rec in record:
    print(rec)
    result_handle = NCBIWWW.qblast("blastp", "nr", rec.format("fasta"))
    print("----------")
    print(result_handle)
    save_file.write(result_handle.read())
    result_handle.close()
else:
    save_file.close()



# if prot_type == 'protein':
#     blast = NcbiblastpCommandline(db="swissprot", evalue='0.0001', remote=remote)()
# elif prot_type == 'nucleotide':
#     blast = NcbiblastnCommandline(db="dbest", evalue='0.0001', remote=remote)()
# else:
#   print('Type must be protein or nucleotide')
#   exit(1)
# input_handle  = open(faa_filename, "r")
# output_handle = open(blast_filename, "w")

# # blast_result_record = NCBIXML.read(StringIO(blast))