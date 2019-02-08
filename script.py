import pandas as pd

class read_csv:
    def __init__(self, sig, names):
        self.sig = sig
        self.names = names

    def rewrite(self):
        self.sig = self.sig.drop(['Gene Name'], axis=1)
        new = pd.concat([self.names, self.sig], axis=1, join='inner')
        new = new.drop_duplicates(subset='Gene Name')

        return new

class up_regulated:
    def __init__(self, up_reg, sig_genes):
        self.up_reg = up_reg
        self.sig_genes = sig_genes

    def up_reg(self):
        self.up_reg.to_csv('up.csv', header=['Gene Name'], index=False)
        self.up_reg = pd.read_csv('up.csv')
        self.up_reg = self.up_reg['Gene Name'].tolist()

        up_regulate = self.sig_genes.loc[self.sig_genes['Probe Set Name']
                                            .isin(self.up_reg)]

        up_regulated = up_regulate.sort_values(by='Fold Change (log2)',
                                                ascending=False)

        up_list = up_regulated['Gene Name']
        up_list = up_list.to_csv('up_reg_list.csv', index=False, header=True)

        return up_regulated.to_csv('up_reg_full.csv', index=False,
                                    header=True)

class down_regulated:
    def __init__(self, down_reg, sig_genes):
        self.down_reg = down_reg
        self.sig_genes = sig_genes

    def down_reg(self):
        self.down_reg = self.down_reg['Gene Name'].tolist()

        down_regulate = self.sig_genes.loc[self.sig_genes['Probe Set Name']
                                            .isin(self.down_reg)]

        down_regulate = down_regulate.drop_duplicates(subset="Gene Name")
        down_regulated = down_regulate.sort_values(by='Fold Change (log2)',
                                                    ascending=True)

        down_list = down_regulated['Gene Name']
        down_list = down_list.to_csv('down_reg_list.csv', index=False,
                                        header='Gene Name')

        return down_regulated.to_csv('down_reg_full.csv', index=False,
                                        header=True)

class resistant_genes:
    def __init__(self, resis, resis_2, up_reg):
        self.resis = resis
        self.resis2 = resis_2
        self.resis_genes = 'resistant_gene_list.csv'
        self.up_reg = up_reg

    def search_genes(self):
        self.resis = pd.DataFrame(self.resis)
        res = self.resis['search_term']
        res.sort_values(ascending=True)

        resis_gene = res.to_csv(self.resis_genes, index=False,
                                    header=['Gene Name'])

        res_genes = pd.read_csv(self.resis_genes)
        res_genes = res_genes.drop_duplicates()
        res_genes.to_csv(self.resis_genes, index=False)

        resis_genes = self.up_reg.loc[self.up_reg['Gene Name'].
                                        isin(self.resis['Gene Name'])]

        resis_genes2 = self.up_reg.loc[self.up_reg['Gene Name'].
                                        isin(self.resis2['Gene'])]

        resis_genes_final = pd.concat([resis_genes, resis_genes2],
                                        ignore_index=True)

        resis_genes_final = resis_genes_final.drop_duplicates()

        list = resis_genes_final['Gene Name']
        list.to_csv('resistant_genes_list.csv', index=False)

        up_regulate = resis_genes_final.sort_values(by='Fold Change (log2)',
                                                        ascending=False)

        return up_regulate.to_csv('resistant_genes.csv',
                                    index=False)

class gene_search:
    def __init__(self, genes, search):
        self.genes = genes
        self.search = search

    def search_into(self):

        for gene in self.search:
            a = gene in self.genes.values
            if a is True:
                print(gene, 'present')
            else:
                print(gene, 'absent')

if __name__ == '__main__':

#MDA_MB_231 DATASET
    sig_mda = pd.read_csv('/Users/lm/datasets/mda/Sig_genes.csv')
    names_mda = pd.read_csv('/Users/lm/datasets/mda/gene_names.txt')

    read = read_csv(sig_mda, names_mda)
    new_file = read_csv.rewrite(read)

    up_reg_mda = pd.read_csv('/Users/lm/datasets/mda/up.csv')
    up_genes = up_regulated(up_reg_mda, new_file)
    up_regulated.up_reg(up_genes)

    down_reg_mda = pd.read_csv('/Users/lm/datasets/mda/down.csv')
    down_genes = down_regulated(down_reg_mda, new_file)
    down_regulated.down_reg(down_genes)

    resis_mda = pd.read_csv('/Users/lm/datasets/mda/resistant_genes.tsv',
                            sep='\t')

    resis_2_mda = pd.read_excel('/Users/lm/datasets/mda/resis_genes_2.xlsx')
    up_reg = pd.read_csv('/Users/lm/datasets/mda/up_reg.csv')

    res_genes = resistant_genes(resis_mda, resis_2_mda, up_reg)
    resistant_genes.search_genes(res_genes)

    gene = pd.read_csv('/Users/lm/datasets/mda/up_lst.csv')
    search = list(map(str, input('Enter Genes to search MDA dataset: ')
                    .upper().split()))

    search_gene = gene_search(gene, search)
    gene_search.search_into(search_gene)


'''
#TAM Datasets
    sig_tam = pd.read_csv('/Users/lm/datasets/tam/Sig_genes.csv')
    names_tam = pd.read_csv('/Users/lm/datasets/tam/gene_names.txt')

    read = read_csv(sig_tam, names_tam)
    new_file = read_csv.rewrite(read)

    up_reg_tam = pd.read_csv('/Users/lm/datasets/tam/up.csv')
    up_genes = up_regulated(up_reg_tam, new_file)
    up_regulated.up_reg(up_genes)

    down_reg_tam = pd.read_csv('/Users/lm/datasets/tam/down.csv')
    down_genes = down_regulated(down_reg_tam, new_file)
    down_regulated.down_reg(down_genes)

    resis_tam = pd.read_csv('/Users/lm/datasets/tam/up_Resistant.tsv',
                            sep='\t')

    up_reg_tam = pd.read_csv('/Users/lm/datasets/tam/up_reg.csv')

    res_genes = resistant_genes(resis_tam, resis_2_mda, up_reg_tam)
    resistant_genes.search_genes(res_genes)

    gene = pd.read_csv('/Users/lm/datasets/tam/up_lst.csv')
    search = list(map(str, input('Enter Genes to search TAM dataset: ').
                        upper().split()))

    search_gene = gene_search(gene, search)
    gene_search.search_into(search_gene)
'''
