
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

    if identpct > 0:
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

        return Query, gene_name,score, group_key,length,strand,querystart,queryend,seq,ident_percentage,line
    else:
        return None


def top_extract(input_file, output_vjc):
    current_query = None
    v_max_line = None
    j_max_line = None
    c_max_line = None
    vmax_line = ''
    jmax_line = ''
    cmax_line = ''
    current_key = None
    c_max = ''
    query=''
    v_max=''
    strand = ''
    v_le = ''
    v_querystart=v_queryend=j_querystart=j_queryend=''
    v_seq=''
    j_seq= ''
    v_pct = ''
    j_pct = ''
    c_pct = ''
    cmax_line = ''
    with open(input_file, 'r') as file, open(output_vjc, 'w') as outvjc:
        outvjc.write('Query'+'\t'+'Gene'+'\t'+'Length'+'\t'+'Score'+'\t'+'Expect'+'\t'+'Identities'+'\t'+'Iden_percentage'+'\t'+'Gaps'+'\t'+'Gap_percentage'+'\t'+'Strand'+'\t'+'querystart'+'\t'+'queryend'+'\t'+'queryseq'+'\t'+'subjectstart'+'\t'+'subjectend'+'\t'+'subjectseq'+'\n')
        for line in file:
            result = process_line(line)
            if result:
                Query, gene_name,score, group_key,length,strand,querystart,queryend,seq,ident_percentage,line = result
                if current_query is None:
                        current_query = Query
                        current_key = group_key
                if Query == current_query:
                    query = current_query
                    if group_key == 'V' and (v_max_line is None or score > v_max_line[2]):
                        v_max_line = (Query,gene_name,score,group_key,length,strand,querystart,queryend,seq,ident_percentage,line)
                    elif group_key == 'J' and (j_max_line is None or score > j_max_line[2]):
                        j_max_line = (Query,gene_name,score,group_key,length,strand,querystart,queryend,seq,ident_percentage,line)
                    elif group_key == 'C' and (c_max_line is None or score > c_max_line[2]):
                        c_max_line = (Query,gene_name,score,group_key,length,strand,querystart,queryend,ident_percentage,line)
                else:
                    query = current_query
             #       if v_max_line: 
                    if v_max_line:
                       vmax_line = v_max_line[10]
                       outvjc.write(vmax_line)
                    v_max = v_max_line[1] if v_max_line else 'NA'
                    strand = v_max_line[5] if v_max_line else 'NA'
                    v_le = v_max_line[4] if v_max_line else 'NA'
                    v_querystart=v_max_line[6] if v_max_line else 'NA'
                    v_queryend=v_max_line[7] if v_max_line else 'NA'
                    v_seq = v_max_line[8] if v_max_line else 'NA'
                    v_pct = v_max_line[9] if v_max_line else 'NA'
             #       if j_max_line:
                    if j_max_line:
                        jmax_line = j_max_line[10]
                        outvjc.write(jmax_line)
            #        jmax_line = j_max_line[9] if j_max_line
                    j_max = j_max_line[1] if j_max_line else 'NA'
                    j_querystart=j_max_line[6] if j_max_line else 'NA'
                    j_queryend=j_max_line[7] if j_max_line else 'NA'
                    j_seq = j_max_line[8] if j_max_line else 'NA'
                    j_pct = j_max_line[9] if j_max_line else 'NA'
               #     if c_max_line:
                    if c_max_line:
                        cmax_line = c_max_line[9]
                        outvjc.write(cmax_line)
                    c_max = c_max_line[1] if c_max_line else 'NA'
                    c_pct = c_max_line[8] if c_max_line else 'NA'
                    

                    current_query = Query
                    current_key = group_key
                    v_max_line = (Query,gene_name,score,group_key,length,strand,querystart,queryend,seq,ident_percentage,line) if group_key == 'V' else None
                    j_max_line = (Query,gene_name,score,group_key,length,strand,querystart,queryend,seq,ident_percentage,line) if group_key == 'J' else None
                    c_max_line = (Query,gene_name,score,group_key,length,strand,querystart,queryend,ident_percentage,line) if group_key == 'C' else None
          

    #    if v_max_line:
   #     max_line = line if v_max_line or j_max_line or c_max_line else 'NA'
        query = v_max_line[0] if v_max_line else 'NA'
        v_max = v_max_line[1] if v_max_line else 'NA'
        v_le = v_max_line[4] if v_max_line else 'NA'
        strand = v_max_line[5] if v_max_line else 'NA'
        v_querystart=v_max_line[6] if v_max_line else 'NA'
        v_queryend=v_max_line[7] if v_max_line else 'NA'
        v_seq = v_max_line[8] if v_max_line else 'NA'
        vmax_line = v_max_line[10] if v_max_line else 'NA'
        outvjc.write(vmax_line)
   #     if j_max_line:
        j_max = j_max_line[1] if j_max_line else 'NA'
        j_querystart=j_max_line[6] if j_max_line else 'NA'
        j_queryend=j_max_line[7] if j_max_line else 'NA'
        j_seq = j_max_line[8] if j_max_line else 'NA'
        jmax_line = j_max_line[10] if j_max_line else 'NA'
        outvjc.write(jmax_line)
 #       if c_max_line:
        c_max = c_max_line[1] if c_max_line else 'NA'
        cmax_line = c_max_line[9] if c_max_line else 'NA'
        outvjc.write(cmax_line)

