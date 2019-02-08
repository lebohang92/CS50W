import pandas as pd
import numpy as np

#Creates new sig dataset with gene name added
class read_csv:
    def __init__(self, sig, names):
        self.sig = sig
        self.gene_names = names

    def out_csv(self):
        print(self.sig.columns)
        self.sig = self.sig.drop([' Gene Name'], axis=1)
        sig_genes = pd.concat([self.gene_names, self.sig], axis=1, join='inner')
        sig_genes = sig_genes.drop_duplicates(subset='Gene Name')

        return sig_genes

#upregulated genes
class up_regulated:
    def __init__(self, up_reg):
        self.up_reg = up_reg
        self.sig = sig

    def up_regulated(self):

        up_reg_genes = pd.read_csv('up.csv')

        read_in = read_csv()
        new = read_in.out_csv()

        up_regulate = self.sig.loc[self.sig['Probe Set Name'].
                            isin(self.up_reg['Gene Name'])]

        up_reg_genes = up_regulate.sort_values(by=' Fold Change (log2)',
                                                ascending=False)

        up_regulate.to_csv('up_reg.csv ', index=False)

        up_list = up_reg_genes['Gene Name']
        up_list = up_list.to_csv('up_lst.csv', index=False, header=True)

        return up_reg_genes

#down regulated genes
class down_regulated:
    def __init__(self):
        self.down_reg = pd.read_csv('down.csv')

    def down_regulated(self):
        self.down_reg = self.down_reg['Gene Name'].tolist()

        read_in = read_csv()
        new = read_in.out_csv()

        down_regulate = new.loc[new['Probe Set Name'].isin(self.down_reg)]
        down_regulate = down_regulate.drop_duplicates(subset="Gene Name")
        down_regulated = down_regulate.sort_values(by=' Fold Change (log2)',
                                                    ascending=True)
        down_regulate.to_csv('down_reg.csv', index=False)

        down_list = down_regulated['Gene Name']
        down_list = down_list.to_csv('down_lst.csv', index=False,
                                        header='Gene Name')

        return down_regulated

#searching through resistant genes
class resistant_genes:
    def __init__(self):
        self.resis = pd.read_csv('resistant_genes.tsv', sep='\t')
        self.resis2 = pd.read_excel('resis_genes_2.xlsx')

    def search_genes(self):
        res = self.resis['search_term']
        res.sort_values(ascending=True)
        resis_gene = res.to_csv('resistant_genes.csv', index=False,
                                header=['Gene Name'])

        res_genes = pd.read_csv('resistant_genes.csv')
        res_genes = res_genes.drop_duplicates()

        res_genes.to_csv('resistant_genes.csv', index=False)
        new = pd.read_csv('up_reg.csv')

        resis_genes = new.loc[new['Gene Name'].isin(res_genes['Gene Name'])]
        resis_genes2 = new.loc[new['Gene Name'].isin(self.resis2['Gene'])]

        resis_genes_final = pd.concat([resis_genes, resis_genes2],
                                    ignore_index=True)

        resis_genes_final = resis_genes_final.drop_duplicates()

        list = resis_genes_final['Gene Name']
        list.to_csv('resistant_genes.csv', index=False, header=['Gene Name'])

        up_regulate = resis_genes_final.sort_values(by=' Fold Change (log2)',
                                                    ascending=False)

        resis_genes_final.to_csv('resistant_gene_list.csv', index=False)

        return resis_genes_final

#searching up regulated genes
class gene_search:
    def __init__(self):
        self.genes = pd.read_csv('up_lst.csv')
        self.search = list(map(str,
                        input('Enter Genes to search: ').upper().split()))

    def search_into(self):

        for gene in self.search:
            a = gene in self.genes.values
            if a is True:
                print(gene, 'present')
            else:
                print(gene, 'absent')


if __name__ == '__main__':
    a = read_csv()
    a.out_csv()

    b = up_regulated()
    b.up_regulated()

    c = down_regulated()
    c.down_regulated()

    d = resistant_genes()
    d.search_genes()

    e = gene_search()
    e.search_into()
