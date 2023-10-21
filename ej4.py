import sys
import os
# Para poder correr esto tuve que instalar emboss local ((linux) sudo apt install emboss

if len(sys.argv) != 3:
    print("ERROR: Parametros invalidos. Se requieren archivo_entrada.fasta archivo_salida.patmatmotifs")
    exit(1)

if sys.argv[1].split(".")[1] != "fasta" or sys.argv[2].split(".")[1] != "patmatmotifs":
    print("ERROR: Parametros invalidos. Se requieren archivo_entrada.fasta archivo_salida.patmatmotifs")
    exit(1)

os.system("transeq -sequence " + sys.argv[1] + " -outseq ej5.fas" )
# os.system("getorf -sequence " + sys.argv[1] + " -outseq ej5.fas" ) # encontrar ORFs en secuencias de nucleótidos
os.system("patmatmotifs -sequence ej5.fas -outfile " + sys.argv[2])