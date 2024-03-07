import re
import pandas as pd
import io
from Bio import SeqIO
from Bio.Seq import Seq




def find_first_of_multiple(seq, substrings):
    for sub in substrings:
        index = seq.find(sub)
        if index != -1:
            return index
    return -1

def translate_with_missing(sequence):

    first_three = str(Seq(sequence[:3]).translate())


    last_three = str(Seq(sequence[-3:]).translate())

    # Translate the middle in chunks of three, marking changes as 'missing'
    middle_chunks = [str(Seq(sequence[i:i+3]).translate()) if sequence[i:i+3] != 'NNN' else 'X' for i in range(3, len(sequence)-3, 3)]

    # Combine
    translated_sequence = first_three + ''.join(middle_chunks) + last_three

    return translated_sequence


def cdr3_extract(df, fasta_file):
    bed_columns = ['query', 'start', 'end', 'strand']

    bed_data = []

    for index, row in df.iterrows():
        v_length = row['v_length']
        query = row['Query']
        v_seq = row['v_seq']
        j_seq = row['j_seq']
        v_querystart = row['v_querystart']
        v_queryend = row['v_queryend']
        j_querystart = row['j_querystart']
        j_queryend = row['j_queryend']
        strand = row['strand']
        start = 0
        end = 0
        if any(pd.isna([v_length, query, v_seq, j_seq, v_querystart, v_queryend, j_querystart, j_queryend, strand])):
            continue
        if v_length > 309:
            if strand == 'Plus/Plus':
                start = v_querystart + 309
                end_index = find_first_of_multiple(j_seq, ['TGG', 'TTT', 'TTC'])
#                end_index = j_seq.find('TGG')
                end = j_querystart + end_index
            elif strand == 'Plus/Minus':
                index = find_first_of_multiple(j_seq, ['CCA', 'AAA', 'GAA'])
                start = j_querystart + index
                end = v_queryend - 309
        if v_length < 500:
            if strand == 'Plus/Plus':
                index_v_seq = find_first_of_multiple(v_seq[260:], ['TGC', 'TGT'])
 #               index_v_seq = v_seq[270:].find('TGC')
                start = v_querystart + 260 +index_v_seq
                index_j_seq = find_first_of_multiple(j_seq, ['TGG', 'TTT', 'TTC'])
#                index_j_seq = j_seq.find('TGG')
                end = j_querystart + index_j_seq + 2
            elif strand == 'Plus/Minus':
                index_j_seq = find_first_of_multiple(j_seq, ['CCA', 'AAA', 'GAA'])
#                index_j_seq = j_seq.find('CCA' if 'CCA' in j_seq else 'AAA') if isinstance(j_seq, str) else -1
                start = j_querystart + index_j_seq
                index_v_seq = find_first_of_multiple(v_seq[0:50], ['GCA', 'ACA'])
#                index_v_seq = v_seq[0:10].find('GCA')
                end = v_querystart + index_v_seq +2

        bed_data.append([query, start, end, strand])

    bed_df = pd.DataFrame(bed_data, columns=bed_columns)
    bed_df = bed_df.dropna(subset=['start', 'end', 'strand'])
#    bed_df.to_csv('cdr3test.bed', sep='\t', header=False, index=False)
    bed_df['start'] = bed_df['start'].astype(int)
    bed_df['end'] = bed_df['end'].astype(int)

    cdr_data = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        query = record.id
        if query in df['Query'].tolist():
 #           sequence = record.seq
            bed_row = bed_df[bed_df['query'] == query]
            if not bed_row.empty:
                start = int(bed_row['start'].iloc[0])
                end = int(bed_row['end'].iloc[0])
                strand = bed_row['strand'].iloc[0]
                sequence = record.seq
                cdr_sequence = str(sequence[start-1:end])
                if len(cdr_sequence) < 200:
                    if strand == 'Plus/Minus':
                        cdr_sequence = str(Seq(cdr_sequence).reverse_complement())
                    aa_sequence = translate_with_missing(cdr_sequence)
                    cdr_data.append([query, cdr_sequence, aa_sequence])

    cdr_df = pd.DataFrame(cdr_data, columns=['Query', 'cdr_sequence', 'aa_sequence'])
    return cdr_df
