import pandas as pd
import numpy as np
import GEOparse
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import quantile_transform
from sklearn.preprocessing import normalize
import rpy2.robjects as robjects
from pathlib import Path

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
        my_file = Path("/Users/lm/documents/datasets/tam/sig_new.csv")
        try:
            my_abs_path = my_file.resolve(strict=True)

        except FileNotFoundError:
            records = GEOparse.get_GEO(self.samples)
            pheno = records.phenotype_data[["title", "source_name_ch1"]]
            print(pheno)
            tam_samples = records.pivot_samples('VALUE')

            mcf10_records = GEOparse.get_GEO(self.mcf10)
            pheno = mcf10_records.phenotype_data[["title", "source_name_ch1"]]
            print(pheno)
            mcf10_samples = mcf10_records.pivot_samples('VALUE')

            n_cases = int(input('Number of case samples: '))
            [self.case.append(input('GSE sample: '))
            for i in range(n_cases)]

            n_controls = int(input('Number of control samples: '))
            [self.control.append(input('GSE sample: '))
            for i in range(n_controls)]

            n_mcf10 = int(input('Number of mcf10 samples: '))
            [self.con_mcf10.append(input('GSE sample: '))
            for i in range(n_mcf10)]

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
        #plt.show()

#Differential expression
class DE:
    def __init__(self):
        self.df = pd.DataFrame()

    def sig_genes(self):
        my_file = Path("/Users/lm/documents/datasets/tam/sig_new.csv")
        try:
            my_abs_path = my_file.resolve(strict=True)

        except FileNotFoundError:
            robjects.r('data_ <- read.csv(file="combined.csv")')
            robjects.r('shape <- dim(data_)')
            robjects.r('data_ <- subset(data_, select = -ID_REF)')
            robjects.r('data_ <- data.matrix(data_)')

            robjects.r('library(limma)')
            robjects.r('design <- model.matrix(~ 0+factor'
                        '(c(1,1,1,2,2,2,3,3,3)))')
            robjects.r('colnames(design) <- c("case", "control", "normal")')
            robjects.r('fit <- lmFit(data_, design)')

            robjects.r('contrast.matrix <- makeContrasts(control-case,'
                                'normal-control, normal-case, levels=design)')

            robjects.r('fit2 <- contrasts.fit(fit, contrast.matrix)')
            robjects.r('fit2 <- eBayes(fit2)')
            robjects.r('table <- topTable(fit2, coef=1, adjust="BH",'
                                'number=Inf)')

            robjects.r('table <- table[order(as.numeric(rownames(table)))'
                                ',,drop=FALSE]')

            robjects.r('write.csv(table, file = "sig_new.csv",'
                        'row.names=FALSE)')

            self.df = pd.read_csv('sig_new.csv')
            #print(self.df.head(10)

#Creates new sig dataset with gene name added
class geneID:
    def __init__(self):
        self.sig = pd.read_csv('sig_new.csv')
        self.probes = pd.read_csv('combined.csv')

    def gene_symb(self):
        df = pd.concat([self.probes[['ID_REF']], self.sig], axis=1, sort=False)
        df.to_csv('new.csv')

        #robjects.r('library(BiomartRT)')
        #robjects.r('affyids = c("202763_at","209310_s_at","207500_at")')
        #robjects.r('getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol"),'
                        #'filters = "affy_hg_u133_plus_2",'
                        #'values = affyids,'
                        #'mart = ensembl)')

if __name__ == '__main__':

#Read GEO files
    samples = 'GSE67916'
    mcf10 = 'GSE15065'
    files = read_geo(samples, mcf10)
    read_geo.read_into(files)
    #read_geo.qcontrol(files)

#Differential expression
    diff = DE()
    diff.sig_genes()

#Reading Sig gene list
    difflst = geneID()
    difflst.gene_symb()



#COMPARING MULTIPLE GEO MATRICES 
> targets <- matrix(c("GSM65523", "noHS", "HS",
>                    "GSM65567", "noHS", "HS"), ncol=3, byrow=TRUE,
>                  dimnames=list(NULL, c("SlideNumber", "Cy3",
"Cy5")))
> design <- modelMatrix(targets, ref="noHS")
> lmFit(exprs, design)
