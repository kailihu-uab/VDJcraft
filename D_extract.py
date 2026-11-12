import math
import os

def extract_sequences(fasta_file, queries_df, dseq_file):
    queries = {}
    for _, row in queries_df.iterrows():
        queries[row['Query']] = (row['dstart'], row['dend'])
    sequences = {}
    with open(fasta_file, 'r') as infile, open(dseq_file, 'w') as outfile:
        sequence = ""
        current_query = ""

        for line in infile:
            if line.startswith('>'):
                if current_query in queries:
                    start, end = queries[current_query]
                    if not (math.isnan(start) or math.isnan(end)):
                         start = int(start)
                         end = int(end)
                         if start >= 1 and start < end:
                               outfile.write(f">{current_query}\n{sequence[start-1:end]}\n")
                current_query = line.strip()[1:]
                sequence = ""
            else:
                sequence += line.strip()

        # Process the last sequence after the loop ends
        if current_query in queries:
            start, end = queries[current_query]
            if not (math.isnan(start) or math.isnan(end)):
                start = int(start)
                end = int(end)
                if start >= 1 and start < end:
                       outfile.write(f">{current_query}\n{sequence[start-1:end]}\n")



def process_line(line, cutoff):
    columns = line.strip().split('\t')

    try:
        Query = columns[0]
        score = float(columns[3])
        gene_name = columns[1]
        length = columns[2]
        seq = columns[12]
        ident_percentage = columns[6].rstrip('%')
        identpct = float(ident_percentage)
    except (ValueError, IndexError):
        # header / blank / malformed line
        return None

    if identpct > cutoff:
        d_col = ['IGHD1', 'IGHD2', 'IGHD3', 'IGHD4','IGHD5','IGHD6','IGHD7','TRBD','TRDD']

        if any(condition in gene_name for condition in d_col):
            group_key = 'D'
        else:
            group_key = None

        return Query, gene_name,score,group_key, seq,line,identpct,length
    else:
        return None


def find_alt(max_line, candidates):
    """Genes (other than the selected one) with the same score, identpct and length as the top hit."""
    if max_line is None:
        return 'NA'
    key = (max_line[2], max_line[6], max_line[7])
    alts = []
    for c in candidates:
        if (c[2], c[6], c[7]) == key and c[1] != max_line[1] and c[1] not in alts:
            alts.append(c[1])
    return ','.join(alts) if alts else 'NA'


def top_extract(input_file,output_file, cutoff, output_alt=None):
    if output_alt is None:
        output_alt = os.path.splitext(output_file)[0] + '_alt.txt'

    current_query = None
    d_max_line = None
    d_lines = []
    current_key = None
    query = ''
    d_max = ''
    d_seq = ''

    with open(input_file,'r') as f,open(output_file,'w') as d, open(output_alt,'w') as dalt:
        d.write('Query'+'\t'+'d_gene'+'\t'+'d_seq'+'\n')
        dalt.write('Query'+'\t'+'alt_d'+'\n')
        for line in f:
            result = process_line(line, cutoff)
            if result:
                Query, gene_name,score, group_key, seq,processed_line,identpct,length = result
                score = float(score)
                if current_query is None or Query != current_query:

                    if current_query is not None:
                        query = current_query
                        d_max = d_max_line[1] if d_max_line else 'NA'
                        d_seq = d_max_line[4] if d_max_line else 'NA'
                        d.write(query+'\t'+d_max+'\t'+d_seq+'\n')
                        dalt.write(query+'\t'+find_alt(d_max_line, d_lines)+'\n')

                    current_query = Query
                    d_max_line = None
                    d_lines = []

                if group_key == 'D':
                    d_line = (Query,gene_name,score,group_key,seq,processed_line,identpct,length)
                    d_lines.append(d_line)
                    if d_max_line is None or score > d_max_line[2]:
                        d_max_line = d_line


        if current_query is not None: 
            query = current_query
            d_max = d_max_line[1] if d_max_line else 'NA'
            d_seq = d_max_line[4] if d_max_line else 'NA'
            d.write(query+'\t'+d_max+'\t'+d_seq+'\n')
            dalt.write(query+'\t'+find_alt(d_max_line, d_lines)+'\n')
