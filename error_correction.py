import pandas as pd




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


