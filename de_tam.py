import pandas as pd
import numpy as np
import os
import GEOparse
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import quantile_transform
from sklearn.preprocessing import normalize
import rpy2.robjects as robjects

#Read .SOFT files from NCBI GEO
class read_geo:
    def __init__(self, samples, mcf10):
        self.samples = samples
        self.mcf10 = mcf10
        self.con_mcf10 = []
        self.case = []
        self.control = []
        self.control_final = pd.DataFrame()
        self.mcf10_final = pd.DataFrame()
        self.case_final = pd.DataFrame()
        self.matrix = pd.DataFrame()

    def read_into(self):
        records = GEOparse.get_GEO(self.samples)
        pheno = records.phenotype_data[["title", "source_name_ch1"]]
        print(pheno)
        tam_samples = records.pivot_samples('VALUE')

        mcf10_records = GEOparse.get_GEO(self.mcf10)
        pheno = mcf10_records.phenotype_data[["title", "source_name_ch1"]]
        print(pheno)
        mcf10_samples = mcf10_records.pivot_samples('VALUE')

        n_cases = int(input('Number of case samples: '))
        [self.case.append(input('GSE sample: ')) for i in range(n_cases)]

        n_controls = int(input('Number of control samples: '))
        [self.control.append(input('GSE sample: ')) for i in range(n_controls)]

        n_mcf10 = int(input('Number of mcf10 samples: '))
        [self.con_mcf10.append(input('GSE sample: ')) for i in range(n_mcf10)]

        self.control_final = tam_samples[self.control]
        self.case_final = tam_samples[self.case]
        self.mcf10_final = mcf10_samples[self.con_mcf10]

        final_dataset = pd.concat([self.case_final, self.control_final,
                                self.mcf10_final], axis=1, sort=False)

        df_or = pd.DataFrame(final_dataset)
        #df_norm = quantile_transform(df_or, n_quantiles=1000,
                                    #random_state=None)

        df_norm = normalize(df_or, norm='l2', axis=1, copy=True,
                                return_norm=False)

        df_log = np.log(df_or)
        df = pd.DataFrame(df_log, columns = df_or.columns,
                            index=tam_samples.index)
        self.matrix = df

        df = df.to_csv('combined.csv')

#Histogram distribution of data after normalisation
    def qcontrol(self):
        self.matrix.hist()
        graph = sns.despine(offset=10, trim=True)
        plt.show()

#Differential expression
class DE:
    def __init__(self):
        self.df = pd.DataFrame()

    def limma(self):

        robjects.r('data_ <- read.csv(file="combined.csv")')
        robjects.r('shape <- dim(data_)')
        robjects.r('data_ <- subset(data_, select = -ID_REF)')
        robjects.r('data_ <- data.matrix(data_)')

        robjects.r('library(limma)')
        robjects.r('design <- model.matrix(~ 0+factor(c(1,1,1,2,2,2,3,3,3)))')
        robjects.r('colnames(design) <- c("case", "control", "normal")')
        robjects.r('fit <- lmFit(data_, design)')

        robjects.r('contrast.matrix <- makeContrasts(control-case,'
                    'normal-control, normal-case, levels=design)')

        robjects.r('fit2 <- contrasts.fit(fit, contrast.matrix)')
        robjects.r('fit2 <- eBayes(fit2)')
        robjects.r('table <- topTable(fit2, coef=1, adjust="BH",'
                    'number=shape[1])')
        robjects.r('table <- table[order(as.numeric(rownames(table)))'
                    ',,drop=FALSE]')

        robjects.r('write.csv(table, file = "sig_new.csv", row.names=FALSE)')

        self.df = pd.read_csv('sig_new.csv')
        #print(self.df.head(10))

#Creates new sig dataset with gene name added
class read_csv:
    def __init__(self):
        self.sig = pd.read_csv('Sig_genes.csv')
        self.gene_names = pd.read_csv('gene_names.txt')

    def out_csv(self):
        print(self.sig.columns)
        self.sig = self.sig.drop(['Gene Name'], axis=1)
        sig_genes = pd.concat([self.gene_names, self.sig], axis=1,
                                join='inner')
        sig_genes = sig_genes.drop_duplicates(subset='Gene Name')
        sig_genes.to_csv('sig_genes_final.csv', index=False)

        return sig_genes

#upregulated genes
class up_regulated:
    def __init__(self):
        self.up_reg = pd.read_csv('up.csv')

    def up_regulated(self):
        self.up_reg.to_csv('up.csv', header=['Gene Name'], index=False)
        up_reg_genes = pd.read_csv('up.csv')

        read_in = read_csv()
        new = read_in.out_csv()

        up_regulate = new.loc[new['Probe Set Name'].
                            isin(up_reg_genes['Gene Name'])]

        up_reg_genes = up_regulate.sort_values(by='Fold Change (log2)',
                                                ascending=False)
        up_regulate.to_csv('up_reg.csv', index=False)

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
        down_regulated = down_regulate.sort_values(by='Fold Change (log2)',
                                                    ascending=True)
        down_regulate.to_csv('down_reg.csv', index=False)

        down_list = down_regulated['Gene Name']
        down_list = down_list.to_csv('down_lst.csv', index=False,
                                        header='Gene Name')

        return down_regulated

#searching through resistant genes
class resistant_genes:
    def __init__(self):
        self.resis = pd.read_csv('/Users/lm/documents/datasets/mda/'
                                    'resistant_genes.tsv', sep='\t')
        self.resis2 = pd.read_excel('/Users/lm/documents/datasets/mda'
                                        '/resis_genes_2.xlsx')
        self.resis3 = pd.read_excel('/Users/lm/documents/datasets/mda'
                                        '/resis_genes_3.xlsx')
        self.resis4 = pd.read_excel('/Users/lm/documents/datasets/mda'
                                        '/resis_genes_4.xlsx')
        self.resis5 = pd.read_excel('/Users/lm/documents/datasets/mda'
                                        '/resis_genes_5.xlsx')

    def search_genes(self):
        print(self.resis5.columns)
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
        resis_genes3 = new.loc[new['Gene Name'].
                                isin(self.resis3['Drug Target'])]
        resis_genes4 = new.loc[new['Gene Name'].
                                isin(self.resis4['Gene Name '])]
        resis_genes5 = new.loc[new['Gene Name'].
                                isin(self.resis5['Gene Name'])]

        resis_genes_final = pd.concat([resis_genes, resis_genes2, resis_genes3,
                                        resis_genes4, resis_genes5],
                                    ignore_index=True)

        resis_genes_final = resis_genes_final.drop_duplicates()
        resis_genes_final = resis_genes_final.sort_values(by=
                                                        'Fold Change (log2)',
                                                        ascending=False)

        list = resis_genes_final['Gene Name']
        list.to_csv('resistant_genes.csv', index=False, header=['Gene Name'])

        up_regulate = resis_genes_final.sort_values(by='Fold Change (log2)',
                                                    ascending=False)

        resis_genes_final.to_csv('resistant_gene_list.csv', index=False)

        return resis_genes_final

#searching up regulated genes
class gene_search:
    def __init__(self):
        self.genes = pd.read_csv('up_lst.csv')
        self.search = list(map(str, input('Enter Genes to search: ').
                            upper().split()))

    def search_into(self):

        for gene in self.search:
            a = gene in self.genes.values
            if a is True:
                print(gene, 'present')
            else:
                print(gene, 'absent')

#searching the 'RNU6' genes
class search_rnu6:
    def __init__(self):
        self.resis_genes = pd.read_csv('resistant_gene_list.csv')
        self.up = pd.read_csv('up_reg.csv')
        self.down = pd.read_csv('down_reg.csv')

    def search_genes(self):
        query_resis = self.resis_genes['Gene Name'].str.contains('RNU6')
        query_up = self.up['Gene Name'].str.contains('RNU6')
        query_down = self.down['Gene Name'].str.contains('RNU6')

        resis = self.resis_genes[query_resis]['Gene Name']
        up = self.up[query_up]['Gene Name']
        final_list = pd.concat([resis, up])
        final_list = final_list.drop_duplicates()

        return final_list.to_csv('rnu6_genes.csv', index=False,
                                    header=['Gene Name'])

if __name__ == '__main__':

#TAM DATASETS
    a = read_csv()
    a.out_csv()

    b = up_regulated()
    b.up_regulated()

    c = down_regulated()
    c.down_regulated()

    d = resistant_genes()
    d.search_genes()

    e = search_rnu6()
    e.search_genes()

    samples = 'GSE67916'
    mcf10 = 'GSE15065'
    files = read_geo(samples, mcf10)
    read_geo.read_into(files)
    read_geo.qcontrol(files)

    f = DE()
    f.limma()

    #print('Because, only I can')
