import math



def extract_sequences(fasta_file, queries_df, dseq_file):
#    dseq_file = 'dseq.fa'
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
                  #       print(f"Extracting for {current_query}, start: {start}, end: {end}")
                         if start >= 1 and start < end:
      #                         print(f"Extracting for {current_query}, start: {start}, end: {end}")
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

def process_line(line):
    columns = line.strip().split('\t')

    Query = columns[0]
    try:
        score = float(columns[3])
    except ValueError:
        return None
    gene_name = columns[1]
    seq = columns[12]
    ident_percentage = columns[6].rstrip('%')
    identpct = float(ident_percentage)

    if identpct > 80:
        d_col = ['IGHD1', 'IGHD2', 'IGHD3', 'IGHD4','IGHD5','IGHD6','IGHD7','TRBD','TRDD']

        if any(condition in gene_name for condition in d_col):
            group_key = 'D'
        else:
            group_key = None

        return Query, gene_name,score,group_key, seq,line
    else:
        return None
    
def top_extract(input_file,output_file):

    current_query = None
    d_max_line = None
    current_key = None
    query = ''
    d_max = ''
    d_seq = ''

    with open(input_file,'r') as f,open(output_file,'w') as d:
        d.write('Query'+'\t'+'d_gene'+'\t'+'d_seq'+'\n')
        for line in f:
            result = process_line(line)
            if result:
                Query, gene_name,score, group_key, seq,processed_line = result
                if current_query is None:
                        current_query = Query
                        current_key = group_key
                if Query == current_query:
                    query == current_query
                    if group_key == 'D' and (d_max_line is None or score > d_max_line[2]):
                        d_max_line = (Query,gene_name,score,group_key,seq,processed_line)
                else:
           #         if d_max_line:
                    query = d_max_line[0] if d_max_line else 'NA'
                    d_max = d_max_line[1] if d_max_line else 'NA'
                    d_seq = d_max_line[4] if d_max_line else 'NA'
                    d.write(query+'\t'+d_max+'\t'+d_seq+'\n')
     #           out = (query+'\t'+d_max+'\t'+d_seq+'\n').strip()
     #           print(out)

                    current_query = Query
                    d_max_line = (Query,gene_name,score,group_key,seq,processed_line) if group_key == 'D' else None


   #     if d_max_line:
        query = d_max_line[0] if d_max_line else 'NA'
        d_max = d_max_line[1] if d_max_line else 'NA'
        d_seq = d_max_line[4] if d_max_line else 'NA'
        d.write(query+'\t'+d_max+'\t'+d_seq+'\n')
   #     lastline = (query+'\t'+d_max+'\t'+d_seq+'\n').strip()
    #    print(d_max_line[1])
    f.close()
    d.close()
