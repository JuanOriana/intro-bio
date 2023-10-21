import sys
import os
from io import StringIO

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline

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
save_file.close()

# Specify the full file path to the query sequence
# query_sequence = os.path.abspath(faa_filename)
# blast_filename = os.path.join("BLAST", os.path.basename(faa_filename).split(".")[0] + ".xml")

# Update the BLAST command with the full file path
# blast = NcbiblastpCommandline(db="swissprot",  evalue='0.001', outfmt=5, out=blast_filename)()
# # input_handle  = open(faa_filename, "r")
# # output_handle = open(blast_filename, "w")
# # output_handle.write(str(blast))
# # output_handle.close()
# # blast_result_record = NCBIXML.read(StringIO(str(blast)))
# # print(blast_result_record)
