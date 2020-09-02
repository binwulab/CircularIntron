"""Convert Excel tables to usable .bed files."""
import pandas as pd
import numpy as np


def df2bed(fname, df, chrom, start, end, name, strand):
    """Turn pandas DataFrame into bed file.

    chrom, start, end, name, strand are all column names.
    """
    with open(fname, 'w') as f:
        df = df.astype(str)
        for i, row in df.iterrows():
            line = [row[chrom],
                    row[start],
                    row[end],
                    row[name],
                    '1',
                    row[strand]]
            f.write('\t'.join(line))
            f.write('\n')


hela_sheet = 'HeLa(hg19)'
hela = pd.read_excel('pnas.1808816115.sd01.xlsx', sheet_name=hela_sheet)
hela['uid'] = hela.index + 2
hela['uid'] = 'hela_' + hela['uid'].astype(str)

sel = hela['BP distance from 3\'ss'] == 'unknown'
hela.loc[sel, 'BP distance from 3\'ss'] = np.nan
hela['BP distance from 3\'ss'] = hela['BP distance from 3\'ss'].astype(float)

sel = ((hela.strand == '+') &
       (~hela['BP distance from 3\'ss'].isna()))
hela.loc[sel, 'end'] = -1*hela.loc[sel, 'BP distance from 3\'ss']\
                       + hela.loc[sel, 'end']

sel = ((hela.strand == '-') &
       (~hela['BP distance from 3\'ss'].isna()))
hela.loc[sel, 'start'] = hela.loc[sel, 'BP distance from 3\'ss']\
                         + hela.loc[sel, 'start']

hela[['start', 'end']] = hela[['start', 'end']].astype(int)

df2bed('talhouarne_hela.bed', hela,
       'chromosome', 'start', 'end', 'uid', 'strand')
