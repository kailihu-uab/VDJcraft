def process_line(line):
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

    if identpct > 80:
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
    
def top_extract(input_file,add_file,output_file, output_vj, out_vgene):
    current_query = None
    v_max_line = None
    j_max_line = None
    c_max_line = None
    current_key = None
    c_max = ''
    query=''
    v_max=''
    strand = ''
    v_le = ''
    v_querystart=v_queryend=j_querystart=j_queryend=''
    v_seq=''
    j_seq= ''
    with open(input_file, 'r') as file, open(add_file, 'w') as add, open(output_file, 'w') as out, open(output_vj, 'w') as outvj, open(out_vgene, 'w') as outv:
        out.write('Query'+'\t'+'v_gene'+'\t'+'j_gene'+'\t'+'c_gene'+'\t'+'v_seq'+'\t'+'j_seq'+'\n')
        outvj.write('Query'+'\t'+'v_gene'+'\t'+'j_gene'+'\t'+'strand'+'\t'+'v_length'+'\t'+'v_querystart'+'\t'+'v_queryend'+'\t'+'j_querystart'+'\t'+'j_queryend'+'\t'+'v_seq'+'\t'+'j_seq'+'\n')
        outv.write('Query'+'\t'+'v_gene'+'\t'+'strand'+'\t'+'v_seq'+'\n')
        for line in file:
            result = process_line(line)
            if result:
                Query, gene_name,score, group_key,length,strand,querystart,queryend,seq,line = result
                if current_query is None:
                        current_query = Query
                        current_key = group_key
                if Query == current_query:
                    query = current_query
                    if group_key == 'V' and (v_max_line is None or score > v_max_line[2]):
                        v_max_line = (Query,gene_name,score,group_key,length,strand,querystart,queryend,seq,line)
                    elif group_key == 'J' and (j_max_line is None or score > j_max_line[2]):
                        j_max_line = (Query,gene_name,score,group_key,length,strand,querystart,queryend,seq,line)
                    elif group_key == 'C' and (c_max_line is None or score > c_max_line[2]):
                        c_max_line = (Query,gene_name,score,group_key,length,strand,querystart,queryend,line)
                else:
                    query = current_query
             #       if v_max_line:
                    if v_max_line:
                        vmax_line = v_max_line[9]
                    v_max = v_max_line[1] if v_max_line else 'NA'
                    strand = v_max_line[5] if v_max_line else 'NA'
                    v_le = v_max_line[4] if v_max_line else 'NA'
                    v_querystart=v_max_line[6] if v_max_line else 'NA'
                    v_queryend=v_max_line[7] if v_max_line else 'NA'
                    v_seq = v_max_line[8] if v_max_line else 'NA'
             #       if j_max_line:
                    if j_max_line:
                        jmax_line = j_max_line[9]
            #        jmax_line = j_max_line[9] if j_max_line
                    j_max = j_max_line[1] if j_max_line else 'NA'
                    j_querystart=j_max_line[6] if j_max_line else 'NA'
                    j_queryend=j_max_line[7] if j_max_line else 'NA'
                    j_seq = j_max_line[8] if j_max_line else 'NA'
               #     if c_max_line:
                    if c_max_line:
                        cmax_line = c_max_line[8]
                    c_max = c_max_line[1] if c_max_line else 'NA'

                    add.write(vmax_line+jmax_line+cmax_line)
                    out.write(query+'\t'+v_max+'\t'+j_max+'\t'+c_max+'\t'+v_seq+'\t'+j_seq+'\n')
                    outvj.write(query+'\t'+v_max+'\t'+j_max+'\t'+strand+'\t'+v_le+'\t'+v_querystart+'\t'+v_queryend+'\t'+j_querystart+'\t'+j_queryend+'\t'+v_seq+'\t'+j_seq+'\n')
                    outv.write(query+'\t'+v_max+'\t'+strand+'\t'+v_seq+'\n')

                    current_query = Query
                    current_key = group_key
                    v_max_line = (Query,gene_name,score,group_key,length,strand,querystart,queryend,seq,line) if group_key == 'V' else None
                    j_max_line = (Query,gene_name,score,group_key,length,strand,querystart,queryend,seq,line) if group_key == 'J' else None
                    c_max_line = (Query,gene_name,score,group_key,length,strand,querystart,queryend,line) if group_key == 'C' else None


    #    if v_max_line:
   #     max_line = line if v_max_line or j_max_line or c_max_line else 'NA'
        query = v_max_line[0] if v_max_line else 'NA'
        v_max = v_max_line[1] if v_max_line else 'NA'
        v_le = v_max_line[4] if v_max_line else 'NA'
        strand = v_max_line[5] if v_max_line else 'NA'
        v_querystart=v_max_line[6] if v_max_line else 'NA'
        v_queryend=v_max_line[7] if v_max_line else 'NA'
        v_seq = v_max_line[8] if v_max_line else 'NA'
   #     if j_max_line:
        j_max = j_max_line[1] if j_max_line else 'NA'
        j_querystart=j_max_line[6] if j_max_line else 'NA'
        j_queryend=j_max_line[7] if j_max_line else 'NA'
        j_seq = j_max_line[8] if j_max_line else 'NA'
 #       if c_max_line:
        c_max = c_max_line[1] if c_max_line else 'NA'
   #     add.write(max_line+'\n')
        out.write(query+'\t'+v_max+'\t'+j_max+'\t'+c_max+'\t'+v_seq+'\t'+j_seq)
        outvj.write(query+'\t'+v_max+'\t'+j_max+'\t'+strand+'\t'+v_le+'\t'+v_querystart+'\t'+v_queryend+'\t'+j_querystart+'\t'+j_queryend+'\t'+v_seq+'\t'+j_seq)
        outv.write(query+'\t'+v_max+'\t'+strand+'\t'+v_seq)

#    print(lastline)
        file.close()
        out.close()
        outvj.close()                 