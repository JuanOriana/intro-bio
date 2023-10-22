import io
import os
from Bio import Entrez, SeqIO
import pandas as pd


def get_ncbi_description(accession: str):
    try:
        handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text")
        record = handle.read()
        handle.close()

        seq_record = SeqIO.read(io.StringIO(record), "genbank")
        description = seq_record.description

        return description
    except Exception as e:
        print(f"Error fetching information for {accession}: {e}")


def get_proteins_descriptions(df: pd.DataFrame):
    names = []
    accessions = []
    descriptions = []

    unique_values = set(
        df[df['Official Symbol Interactor A'] != 'CNGA3']['Official Symbol Interactor A'].unique()) | set(
        df['Official Symbol Interactor B'].unique())
    print("Distinct proteins:", len(unique_values))

    for value in unique_values:
        sub_df_a = df[(df['Official Symbol Interactor A'] == value) & (df['Official Symbol Interactor A'] != 'CNGA3')]
        sub_df_b = df[(df['Official Symbol Interactor B'] == value)]

        if not sub_df_a.empty:
            accession = sub_df_a.iloc[0]['REFSEQ Accessions Interactor A'].split('|')[0]
            description = get_ncbi_description(accession)
        elif not sub_df_b.empty:
            accession = sub_df_b.iloc[0]['REFSEQ Accessions Interactor B'].split('|')[0]
            description = get_ncbi_description(accession)
        else:
            accession = ""
            description = ""

        names.append(value)
        accessions.append(accession)
        descriptions.append(description)

    new_df = pd.DataFrame({'name': names, 'accession': accessions, 'description': descriptions})

    return new_df


def main():
    df = pd.read_csv("BioGRID/BIOGRID-GENE.csv", delimiter='\t', encoding='utf-8')
    df = df.drop_duplicates()
    print("Total number of interactions:", len(df))

    output_filename = './output/proteins_descriptions.csv'
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)

    get_proteins_descriptions(df).to_csv(output_filename, index=False)


if __name__ == '__main__':
    Entrez.email = "pdomingues@itba.edu.ar"
    main()
