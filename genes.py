import pandas as pd
import numpy as np

sig = pd.read_csv('sig_genes.csv')
#print(up_reg.columns)

sig = sig.drop(['Gene Name'], axis=1)
#print(up_reg.columns)

#genes = up_reg[[' Gene Name']].copy()
#genes.to_csv('genes.txt', index=False, header=False)

gene_names = pd.read_csv('genenames.txt')

new = pd.concat([gene_names, sig], axis=1, join='inner')

#print(new.columns)
#print(new)
#print(new.head(100))

up_reg = pd.read_csv('up.csv')
up_reg = up_reg['Gene Name'].tolist()

down_reg = pd.read_csv('down.csv')
down_reg = down_reg['Gene Name'].tolist()

#TOP 50 DIFFERENTIALLY EXPRESSED GENES !!!
#print(new['Probe Set Name'])
up_regulate = new.loc[new['Probe Set Name'].isin(up_reg)]

#remove duplicates
up_regulate = up_regulate.drop_duplicates(subset='Gene name')
up_regulated = up_regulate.sort_values(by='Fold Change (log2)', ascending=False)
print(up_regulated.head(50), '\n')
print(up_regulated.shape)
#up_regulate.to_csv('up.txt')

down_regulate = new.loc[new['Probe Set Name'].isin(down_reg)]

#remove duplicates
down_regulate = down_regulate.drop_duplicates(subset="Gene name")
down_regulated = down_regulate.sort_values(by='Fold Change (log2)', ascending=True)
print(down_regulated.head(50),'\n')
print(down_regulated.shape)
#down_regulate.to_csv('down.txt')
