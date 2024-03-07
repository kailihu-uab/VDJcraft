import pandas as pd


#file1_path = 'HG0268/vdjc_original.csv'

#df1 = pd.read_csv(file1_path, sep='\t', header=0) #names=['query', 'chr','position', 'v_gene', 'd_gene', 'j_gene','c_gene'])


def error_correction(df):
    df = df.sort_values(by='position')

    
    clusters = []

    # Iterate through cluster by position with a difference range of 100
    current_cluster = []
    for index, row in df.iterrows():
        if not current_cluster or abs(row['position'] - current_cluster[-1]['position']) <= 100:
            current_cluster.append(row)
        else:
            clusters.append(pd.DataFrame(current_cluster))
            current_cluster = [row]

    if current_cluster:
        clusters.append(pd.DataFrame(current_cluster))

    # Iterate through clusters
    for cluster in clusters:
        common_values = cluster[['v_gene', 'd_gene', 'j_gene', 'c_gene']].mode().loc[0].to_dict()

        differing_rows = cluster[(cluster[['v_gene', 'd_gene', 'j_gene', 'c_gene']] != common_values).sum(axis=1) == 1]

        for index in differing_rows.index:
            cluster.loc[index, ['v_gene', 'd_gene', 'j_gene', 'c_gene']] = common_values

        # Update
        df.update(cluster)
        df['position'] = df['position'].astype(int)
    return df

#with open('HG0268/vdjc_original.csv','r') as vdjc, open('HG0268/table.csv','r') as T:
#        df_vdjc = pd.read_csv(vdjc, sep='\t', header=0)

#        df_pos = pd.read_csv(T, sep='\t', header=0)

#        pos_vdjc = pd.merge(df_pos, df_vdjc, on = 'Query', how='inner')

#        pos_vdjc.to_csv('0268pos.txt', sep='\t', index=False)

#        corrected_df = error_correction(pos_vdjc)

    #    corrected_df['position'] = corrected_df['position'].astype(int)
#        corrected_df.to_csv('0268correct.txt', sep='\t', index=False)
# Apply the function
#result_df = error_correction(df1)
#result_df['position'] = result_df['position'].astype(int)
#result_df.to_csv('2106correct.txt', sep='\t', index=False)
# Print the result
#print(result_df)

