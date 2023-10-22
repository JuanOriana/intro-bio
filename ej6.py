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

        comments = seq_record.annotations.get("comment", "")
        summary = ""

        for paragraph in comments.split('.\n'):
            if paragraph.startswith("Summary:"):
                summary = paragraph.strip("Summary: ").strip()
                break

        return description, summary
    except Exception as e:
        print(f"Error fetching information for {accession}: {e}")
        return None, None


def get_proteins_descriptions(df: pd.DataFrame):
    unique_values = set(
        df[df['Official Symbol Interactor A'] != gene_name]['Official Symbol Interactor A'].unique()) | set(
        df['Official Symbol Interactor B'].unique())
    print("Distinct proteins:", len(unique_values))

    proteins_info = []

    for value in unique_values:
        sub_df_a = df[(df['Official Symbol Interactor A'] == value) & (df['Official Symbol Interactor A'] != gene_name)]
        sub_df_b = df[(df['Official Symbol Interactor B'] == value)]

        if not sub_df_a.empty:
            accession = sub_df_a.iloc[0]['REFSEQ Accessions Interactor A'].split('|')
        elif not sub_df_b.empty:
            accession = sub_df_b.iloc[0]['REFSEQ Accessions Interactor B'].split('|')
        else:
            accession = ""

        for unique_accession in accession:
            description, summary = get_ncbi_description(unique_accession)
            if summary:
                break

        proteins_info.append({'name': value, 'description': description, 'summary': summary})

    new_df = pd.DataFrame(proteins_info)

    return new_df


def main():
    df = pd.read_csv("BioGRID/BIOGRID-GENE.csv", delimiter='\t', encoding='utf-8')
    df = df.drop_duplicates()
    print("Total number of interactions:", len(df))

    output_filename = './output/proteins_descriptions.csv'
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)

    get_proteins_descriptions(df).to_csv(output_filename, index=False)


if __name__ == '__main__':
    gene_name = 'CNGA3'
    Entrez.email = "pdomingues@itba.edu.ar"
    main()
