#import pandas as pd
#import argparse
import os

#parser = argparse.ArgumentParser()

#parser.add_argument('-p', '--VJIdent', type=int, default = 90, help="Enter preferred identifation percentage for VJ genes")
#parser.add_argument('-pd', '--DIdent', type=int, default = 95, help="Enter preferred identifation percentage for D genes")

#args = parser.parse_args()

def process_line(line,cutoff):
    columns = line.strip().split('\t')

    try:
        Query = columns[0]
        score = float(columns[3])
        gene_name = columns[1]
        length = columns[2]
        strand = columns[9]
        querystart = columns[10]
        queryend = columns[11]
        seq = columns[12]
        ident_percentage = columns[6].rstrip('%')
        identpct = float(ident_percentage)
    except (ValueError, IndexError):
        # header / blank / malformed line
        return None

    if identpct > cutoff:
        v_col = ['IGHV', 'IGKV','IGLV', 'TRAV', 'TRBV', 'TRGV', 'TRDV']
        j_col = ['IGHJ', 'IGKJ','IGLJ', 'TRAJ', 'TRBJ', 'TRGJ', 'TRDJ']
        c_col = ['IGHA', 'IGHE', 'IGHG', 'IGHD chr14','IGHM','IGLC','IGKC','TRGC','TRBC','TRDC','TRAC']

        if any(condition in gene_name for condition in v_col):
            group_key = 'V'
        elif any(condition in gene_name for condition in j_col):
            group_key = 'J'
        elif any(condition in gene_name for condition in c_col):
            group_key = 'C'
        else:
            group_key = None

        return Query, gene_name,score, group_key,length,strand,querystart,queryend,seq,identpct,line
    else:
        return None


def update_line(prev_line, score, identpct):
    if prev_line is None:
        return True

    prev_score = prev_line[2]
    prev_ident = prev_line[9]

    if score > prev_score:
        return True

    if 500 < score <= prev_score and identpct > prev_ident:
        return True

    return False


def find_alt(max_line, candidates):
    """Genes (other than the selected one) with the same score, identpct and length as the top hit."""
    if max_line is None:
        return 'NA'
    key = (max_line[2], max_line[9], max_line[4])
    alts = []
    for c in candidates:
        if (c[2], c[9], c[4]) == key and c[1] != max_line[1] and c[1] not in alts:
            alts.append(c[1])
    return ','.join(alts) if alts else 'NA'


def top_extract(input_file,add_file,output_file, output_vj, out_vgene, cutoff, score_cutoff=300, output_alt=None):
    if output_alt is None:
        output_alt = os.path.splitext(output_file)[0] + '_alt.txt'

    current_query = None
    v_max_line = None
    j_max_line = None
    c_max_line = None
    v_lines, j_lines = [], []

    with open(input_file, 'r') as file, open(add_file, 'w') as add, open(output_file, 'w') as out, open(output_vj, 'w') as outvj, open(out_vgene, 'w') as outv, open(output_alt, 'w') as outalt:
        out.write('Query'+'\t'+'v_gene'+'\t'+'j_gene'+'\t'+'c_gene'+'\t'+'v_seq'+'\t'+'j_seq'+'\n')
        outvj.write('Query'+'\t'+'v_gene'+'\t'+'j_gene'+'\t'+'strand'+'\t'+'v_length'+'\t'+'v_querystart'+'\t'+'v_queryend'+'\t'+'j_querystart'+'\t'+'j_queryend'+'\t'+'v_seq'+'\t'+'j_seq'+'\n')
        outv.write('Query'+'\t'+'v_gene'+'\t'+'strand'+'\t'+'v_seq'+'\n')
        outalt.write('Query'+'\t'+'alt_v'+'\t'+'alt_j'+'\n')

        def write_query():
            # Write outputs for current_query (shared by the loop and the last query)
            if v_max_line is not None and v_max_line[2] <= score_cutoff:
                return  # skip this query

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

            add.write(f"{v_max_line[10].rstrip() if v_max_line else ''}\t{j_max_line[10].rstrip() if j_max_line else ''}\t{c_max_line[10].rstrip() if c_max_line else ''}\n")
            out.write(f"{current_query}\t{v_max}\t{j_max}\t{c_max}\t{v_seq}\t{j_seq}\n")
            outvj.write(f"{current_query}\t{v_max}\t{j_max}\t{strand}\t{v_le}\t{v_querystart}\t{v_queryend}\t{j_querystart}\t{j_queryend}\t{v_seq}\t{j_seq}\n")
            outv.write(f"{current_query}\t{v_max}\t{strand}\t{v_seq}\n")
            outalt.write(f"{current_query}\t{find_alt(v_max_line, v_lines)}\t{find_alt(j_max_line, j_lines)}\n")

        for line in file:
            result = process_line(line,cutoff)
            if result:
                Query, gene_name,score,group_key,length,strand,querystart,queryend,seq,identpct,line = result
                if current_query is None or Query != current_query:
                    # Output the data for the previous query
                    if current_query is not None:
                        write_query()

                    # Reset for the new query
                    current_query = Query
                    v_max_line = None
                    j_max_line = None
                    c_max_line = None
                    v_lines, j_lines = [], []

                # Update max lines for the current query
      #          if group_key == 'V' and (v_max_line is None or score > v_max_line[2] or identpct > v_max_line[9]):
                if group_key == 'V':
                    v_lines.append(result)
                    if update_line(v_max_line, score, identpct):
                        v_max_line = result
                elif group_key == 'J':
                    j_lines.append(result)
                    if j_max_line is None or score > j_max_line[2]:
                        j_max_line = result
                elif group_key == 'C' and (c_max_line is None or score > c_max_line[2]):
                    c_max_line = result

        # Output the data for the last query
        if current_query is not None:
            write_query()
