import re
import pandas as pd
import io
from Bio import SeqIO
from Bio.Seq import Seq


def extract_cdr1_sequence(df):
    result = []

    for index, row in df.iterrows():
        query = row['Query']
        v_gene = row['v_gene']
        v_seq = row['v_seq']
        strand = row['strand']

        start_index = -1
        end_index = -1
        cdr1_sequence = ""
        if strand == 'Plus/Plus':
            start_range = range(60, 91)
            end_range = range(100, 121)

            for i in start_range:
                if isinstance(v_seq, str) and v_seq[i:i + 3] in ['TGC', 'TGT']:
                    start_index = i+12
                    break

            for i in end_range:
                if isinstance(v_seq, str) and v_seq[i:i + 3] == 'TGG':
                    end_index = i-6
                    break

        elif strand == 'Plus/Minus':
            start_range = range(170, 190)
            end_range = range(220, 240)


            for i in start_range:
                if isinstance(v_seq, str) and v_seq[i:i + 3] == 'CCA':
                    start_index = i+9  # Adjust to the actual index in v_seq
                    break

            for i in end_range:
                if isinstance(v_seq, str) and v_seq[i:i + 3] in ['GCA', 'ACA']:
                    end_index = i-7  # Adjust to the actual index in v_seq
                    break

        if start_index != -1 and end_index != -1:
            cdr1_sequence = v_seq[start_index:end_index]
   #     else:
   #         if strand == 'Plus/Plus':
   #             cdr1_sequence = v_seq[75:99]
   #         elif strand == 'Plus/Minus':
   #             cdr1_sequence = v_seq[len(v_seq) - 194:len(v_seq) - 218]

        result.append([query, strand, v_gene, cdr1_sequence])

    result_df = pd.DataFrame(result, columns=['query', 'strand', 'v_gene', 'cdr1_sequence'])
    return result_df


def cdr2_extract(df, column1, column2, column3,column4,column5,column6):
    cdr2=''
    cdr2rows = []
    for index, row in df.iterrows():

        val1 = str(row[column1])
        val2 = str(row[column2])
        val3 = str(row[column3])
        val4 = str(row[column4])
        val5 = str(row[column5])
        val6 = str(row[column6])

        if val2 == 'Plus/Plus':
       #     val4_head1 = val4[0:6]
            start_index1 = val3[150:180].find(val4)
            a=len(val4)
            if start_index1 != -1:
                cdr2 = (val3[start_index1+150:start_index1+150+a])
            else:
                start_range = range(100, 121)
                for i in start_range:
                    if isinstance(val3, str) and val3[i:i + 3] == 'TGG':
                        start_index = i+45
                        a=len(val4)
                        cdr2 = (val3[start_index:start_index+a])
                        break
                else:
                    cdr2 = 'NA'
            if not hasattr(val1, '__iter__'):
                val1 = [val1]
            if not hasattr(cdr2, '__iter__'):
                cdr2 = [cdr2]
        if val2 == 'Plus/Minus':
        #    val5_head1 = val5[0:6]
            start_index1 = val3[120:180].find(val5)
            a=len(val5)
            if start_index1 != -1:
                cdr2 = (val3[start_index1+120:start_index1+120+a])
            else:
                end_range = range(len(val3) - 120, len(val3) - 100)
                for i in end_range:
                    if isinstance(val3, str) and val3[i:i + 3] == 'CCA':
                        end_index = i-42
                        a=len(val5)
                        cdr2 = (val3[end_index-a:end_index])
                        break
                else:
                    cdr2 = 'NA'
            if not hasattr(val1, '__iter__'):
                val1 = [val1]
            if not hasattr(cdr2, '__iter__'):
                cdr2 = [cdr2]
            cdr2 = {'Query': val6, 'v_gene': val1, 'cdr2': cdr2}
            cdr2rows.append(cdr2)
            dfcdr2 = pd.DataFrame(cdr2rows)
    return dfcdr2


