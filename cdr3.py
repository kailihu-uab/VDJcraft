import re
import pandas as pd
import io
from Bio import SeqIO
from Bio.Seq import Seq



def find_first_of_multiple(seq, substrings):
    indices = []                        
    for sub in substrings:              
        i = seq.find(sub)               
        if i != -1:                     
            indices.append(i)           
    return min(indices) if indices else -1   

   
def find_last_of_multiple(seq, substrings):
    last_index = -1
    for sub in substrings:
        idx = seq.rfind(sub)
        if idx > last_index:
            last_index = idx
    return last_index
 

def safe_translate(seq):
    clean_seq = seq.replace('-', '').replace('N', 'A')
    n = len(clean_seq) - (len(clean_seq) % 3)
    if n <= 0:
        return ''
    return str(Seq(clean_seq[:n]).translate())



def translate_with_missing(sequence):

    first_three = str(Seq(sequence[:3]).translate())


    last_three = str(Seq(sequence[-3:]).translate())

    # Translate the middle in chunks of three, marking changes as 'missing'
    middle_chunks = [str(Seq(sequence[i:i+3]).translate()) if sequence[i:i+3] != 'NNN' else 'X' for i in range(3, len(sequence)-3, 3)]

    # Combine
    translated_sequence = first_three + ''.join(middle_chunks) + last_three

    return translated_sequence



def clean_cdr3_dataframe(df):
    def clean_row(row):
        query = row['Query']
        cdr3_seq = row['cdr_sequence']
        cdr3_aa = row['aa_sequence']

        if "*" not in cdr3_aa:
            return pd.Series([query, cdr3_seq, cdr3_aa, 'original'])

        seq = cdr3_seq.upper()
        start_codons = ["TGT", "TGC"]
        end_codons = ["TGG", "TAT", "TAC"]

        fixed_seq = None
        for start in start_codons:
            start_idx = seq.find(start)
            if start_idx != -1:
                for end in end_codons:
                    end_idx = seq.rfind(end)
                    if end_idx > start_idx and (end_idx - start_idx) % 3 == 0:
                        fixed_seq = seq[start_idx:end_idx+3]
                        break
            if fixed_seq:
                break
        if fixed_seq:
            fixed_aa = str(Seq(fixed_seq).translate(to_stop=False))
            if "*" not in fixed_aa:
                return pd.Series([query, fixed_seq, fixed_aa, 'motif_repaired'])

        for shift in [1, 2]:
            shifted_seq = seq[shift:]
            alt_aa = str(Seq(shifted_seq).translate(to_stop=False))
            if "*" not in alt_aa:
                return pd.Series([query, shifted_seq, alt_aa, 'frameshifted'])

        fixed_aa = cdr3_aa.replace("*", "?")
        return pd.Series([query, cdr3_seq, fixed_aa, 'star_replaced'])

    cleaned_df = df.apply(clean_row, axis=1)
    cleaned_df.columns = ['Query', 'cdr_sequence', 'aa_sequence', 'cdr3_fix_status']
    return cleaned_df



def cdr3_extract(df, fasta_file, output_bed, output_cdr3):
    df = pd.read_csv(df, sep='\t')
    df.columns = df.columns.str.strip()


    bed_columns = ['query', 'start', 'end', 'strand']

    bed_data = []

    for index, row in df.iterrows():
        v_length = row['v_length']
        query = row['Query']
        v_gene = row['v_gene']
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
        if v_length < 500:
            if strand == 'Plus/Plus':
                v_seq = v_seq.replace('-', '').replace('N', 'A')
                index_v_seq = find_last_of_multiple(v_seq[150:], ['TATTGC', 'TATTGT', 'TACTGC', 'TACTGT','TTCTGT'])
                if index_v_seq != -1:
                    start = v_querystart + 150 + index_v_seq + 3
                else:
                    index_v_seq = find_last_of_multiple(v_seq[150:], ['TGC', 'TGT'])
                    if index_v_seq != -1:
                       start = v_querystart + 150 + index_v_seq
                    else:
                        if v_gene.startswith(('IGKV', 'IGLV', 'TRAV', 'TRGV')):
                            start = v_querystart + 262
                        elif v_gene.startswith(('IGHV', 'TRBV', 'TRDV')):
                            start = v_querystart + 274
                        else:
                            start = None

                j_seq = j_seq.replace('-', '').replace('N', 'A')
                index_j_seq = find_first_of_multiple(j_seq, ['TGGGG', 'TTTGG', 'TTCGG'])
                if index_j_seq != -1:
                    end = j_querystart + index_j_seq + 2
                else:
                    index_j_seq = find_first_of_multiple(j_seq, ['TGG', 'TTT', 'TTC'])
                    if index_v_seq != -1:
                        end = j_querystart + index_j_seq + 2
                    else:
                        start = None
            elif strand == 'Plus/Minus':
                j_seq = j_seq.replace('-', '').replace('N', 'A')
                index_j_seq = find_last_of_multiple(j_seq, ['CCCCA', 'CCAAA', 'CCGAA'])
                if index_j_seq != -1:
                    start = j_querystart + index_j_seq + 2
                else:
                    index_j_seq = find_last_of_multiple(j_seq, ['CCA', 'AAA', 'GAA'])
                    start = j_querystart + index_j_seq
                v_seq = v_seq.replace('-', '').replace('N', 'A')
                index_v_seq = find_first_of_multiple(v_seq[0:60], ['GCAGTA', 'ACAGTA'])
                if index_v_seq != -1:
                    end = v_querystart + index_v_seq + 2
                else:
                    index_v_seq = find_first_of_multiple(v_seq[0:60], ['GCA', 'ACA'])
                    if index_v_seq != -1:
                        end = v_querystart + index_v_seq + 2
                    else:
                        if v_gene.startswith(('IGKV', 'IGLV', 'TRAV', 'TRGV')):
                            end = v_queryend - 261
                        elif v_gene.startswith(('IGHV', 'TRBV', 'TRDV')):
                            end = v_queryend - 273                    
                        else:
                            end = None

        bed_data.append([query, start, end, strand])

    bed_df = pd.DataFrame(bed_data, columns=bed_columns)
    bed_df = bed_df.dropna(subset=['start', 'end', 'strand'])
    bed_df['start'] = bed_df['start'].astype(int)
    bed_df['end'] = bed_df['end'].astype(int)
    bed_df.to_csv(output_bed, sep='\t', header=False, index=False)
    cdr_data = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        query = record.id
        if query in df['Query'].tolist():
            bed_row = bed_df[bed_df['query'] == query]
            if not bed_row.empty:
                start = int(bed_row['start'].iloc[0])
                end = int(bed_row['end'].iloc[0])
                strand = bed_row['strand'].iloc[0]
                sequence = record.seq
                cdr_sequence = str(sequence[start-1:end])
                if len(cdr_sequence) < 100:
                    if strand == 'Plus/Minus':
                        cdr_sequence = str(Seq(cdr_sequence).reverse_complement())
                    aa_sequence = translate_with_missing(cdr_sequence)
                    cdr_data.append([query, cdr_sequence, aa_sequence])

    cdr_df = pd.DataFrame(cdr_data, columns=['Query', 'cdr_sequence', 'aa_sequence'])
    cdr_df_revisit = clean_cdr3_dataframe(cdr_df)
    cdr_df_revisit = cdr_df_revisit[cdr_df_revisit['cdr_sequence'] != '']
    cdr_df_revisit.to_csv(output_cdr3, sep='\t', index=False, na_rep='out of frame')
    return cdr_df_revisit


