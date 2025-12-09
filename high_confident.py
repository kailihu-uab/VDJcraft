#import pandas as pd
#import argparse

#parser = argparse.ArgumentParser()

#parser.add_argument('-p', '--VJIdent', type=int, default = 90, help="Enter preferred identifation percentage for VJ genes")
#parser.add_argument('-pd', '--DIdent', type=int, default = 95, help="Enter preferred identifation percentage for D genes")

#args = parser.parse_args()

def process_line(line,cutoff):
    columns = line.strip().split('\t')

    Query = columns[0]
    try:    
        score = float(columns[3])
    except ValueError:
        return None
    gene_name = columns[1]
    length = columns[2]
    strand = columns[9]
    querystart = columns[10]
    queryend = columns[11]
    seq = columns[12]
    ident_percentage = columns[6].rstrip('%')
    identpct = float(ident_percentage)

    if identpct > cutoff:
        v_col = ['IGHV', 'IGKV','IGLV', 'TRAV', 'TRBV', 'TRGV', 'TRDV']
        j_col = ['IGHJ', 'IGKJ','IGLJ', 'TRAJ', 'TRBJ']
        c_col = ['IGHA', 'IGHE', 'IGHG', 'IGHD chr14','IGHM','IGLC','IGKC','TRGC','TRBC','TRDC','TRAC']

        if any(condition in gene_name for condition in v_col):
            group_key = 'V'
        elif any(condition in gene_name for condition in j_col):
            group_key = 'J'
        elif any(condition in gene_name for condition in c_col):
            group_key = 'C'
        else:
            group_key = None

        return Query, gene_name,score, group_key,length,strand,querystart,queryend,seq,line
    else:
        return None



def top_extract(input_file,add_file,output_file, output_vj, out_vgene, cutoff):
    current_query = None
    v_max_line = None
    j_max_line = None
    c_max_line = None
    current_key = None
    c_max = ''
    query=''
    v_max=''
    j_max=''
    strand = ''
    v_le = ''
    v_querystart=v_queryend=j_querystart=j_queryend=''
    v_seq=''
    j_seq= ''
    cmax_line = ''
    with open(input_file, 'r') as file, open(add_file, 'w') as add, open(output_file, 'w') as out, open(output_vj, 'w') as outvj, open(out_vgene, 'w') as outv:
        out.write('Query'+'\t'+'v_gene'+'\t'+'j_gene'+'\t'+'c_gene'+'\t'+'v_seq'+'\t'+'j_seq'+'\n')
        outvj.write('Query'+'\t'+'v_gene'+'\t'+'j_gene'+'\t'+'strand'+'\t'+'v_length'+'\t'+'v_querystart'+'\t'+'v_queryend'+'\t'+'j_querystart'+'\t'+'j_queryend'+'\t'+'v_seq'+'\t'+'j_seq'+'\n')
        outv.write('Query'+'\t'+'v_gene'+'\t'+'strand'+'\t'+'v_seq'+'\n')
        for line in file:
            result = process_line(line,cutoff)
            if result:
                Query, gene_name,score, group_key,length,strand,querystart,queryend,seq,line = result
                if current_query is None or Query != current_query:
                    # Output the data for the previous query
                    if current_query is not None:
                        # Write outputs for the previous query
                        v_max = v_max_line[1] if v_max_line else 'NA'
                        v_seq = v_max_line[8] if v_max_line else 'NA'
                        strand = v_max_line[5] if v_max_line else 'NA'
                        v_le = v_max_line[4] if v_max_line else 'NA'
                        v_querystart = v_max_line[6] if v_max_line else 'NA'
                        v_queryend = v_max_line[7] if v_max_line else 'NA'

                        j_max = j_max_line[1] if j_max_line else 'NA'
                        j_seq = j_max_line[8] if j_max_line else 'NA'
                        j_querystart = j_max_line[6] if j_max_line else 'NA'
                        j_queryend = j_max_line[7] if j_max_line else 'NA'

                        c_max = c_max_line[1] if c_max_line else 'NA'

                        add.write(f"{v_max_line[9] if v_max_line else ''}\t{j_max_line[9] if j_max_line else ''}\t{c_max_line[8] if c_max_line else ''}\n")
                        out.write(f"{current_query}\t{v_max}\t{j_max}\t{c_max}\t{v_seq}\t{j_seq}\n")
                        outvj.write(f"{current_query}\t{v_max}\t{j_max}\t{strand}\t{v_le}\t{v_querystart}\t{v_queryend}\t{j_querystart}\t{j_queryend}\t{v_seq}\t{j_seq}\n")
                        outv.write(f"{current_query}\t{v_max}\t{strand}\t{v_seq}\n")

                    # Reset for the new query
                    current_query = Query
                    v_max_line = None
                    j_max_line = None
                    c_max_line = None

                # Update max lines for the current query
                if group_key == 'V' and (v_max_line is None or score > v_max_line[2]):
                    v_max_line = (Query, gene_name, score, group_key, length, strand, querystart, queryend, seq, line)
                elif group_key == 'J' and (j_max_line is None or score > j_max_line[2]):
                    j_max_line = (Query, gene_name, score, group_key, length, strand, querystart, queryend, seq, line)
                elif group_key == 'C' and (c_max_line is None or score > c_max_line[2]):
                    c_max_line = (Query, gene_name, score, group_key, length, strand, querystart, queryend, line) 

        # Output the data for the last query
        if current_query is not None:
                    v_max = v_max_line[1] if v_max_line else 'NA'
                    v_seq = v_max_line[8] if v_max_line else 'NA'
                    strand = v_max_line[5] if v_max_line else 'NA'
                    v_le = v_max_line[4] if v_max_line else 'NA'
                    v_querystart = v_max_line[6] if v_max_line else 'NA'
                    v_queryend = v_max_line[7] if v_max_line else 'NA'

                    j_max = j_max_line[1] if j_max_line else 'NA'
                    j_seq = j_max_line[8] if j_max_line else 'NA'
                    j_querystart = j_max_line[6] if j_max_line else 'NA'
                    j_queryend = j_max_line[7] if j_max_line else 'NA'

                    c_max = c_max_line[1] if c_max_line else 'NA'



                    add.write(f"{v_max_line}\t{j_max_line}\t{c_max_line}\n" if v_max_line and j_max_line and c_max_line else "NA\n")
                    out.write(f"{current_query}\t{v_max}\t{j_max}\t{c_max}\t{v_seq}\t{j_seq}\n")
                    outvj.write(f"{current_query}\t{v_max}\t{j_max}\t{strand}\t{v_le}\t{v_querystart}\t{v_queryend}\t{j_querystart}\t{j_queryend}\t{v_seq}\t{j_seq}\n")
                    outv.write(f"{current_query}\t{v_max}\t{strand}\t{v_seq}\n")



