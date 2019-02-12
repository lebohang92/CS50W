import pandas as pd
import numpy as np
import GEOparse
from sklearn.preprocessing import quantile_transform
import rpy2.robjects as robjects
robjects.r('library(limma)')
robjects.r('library(hgu133plus2.db)') #MCF_7_TAM
robjects.r('library(pd.hugene.1.0.st.v1)') #MDA
robjects.r('library(illuminaHumanv4.db)') #TNBC
robjects.r('library(illuminaHumanWGDASLv3.db)') #MCF_7_Stem_cells, could be not sure...
from pathlib import Path

#Reading .SOFT files from NCBI GEO
class read_geo:
    def __init__(self, exp, con):
        self.exp = exp
        self.con = con

    def read_files(self):
        my_file = Path("/Users/lm/documents/datasets/project/diff/dataset.csv")

        try:
            my_abs_path = my_file.resolve(strict=True)

        except:
            exper = GEOparse.get_GEO(self.exp)
            pheno = exper.phenotype_data[["title", "source_name_ch1"]]
            print(pheno)
            exper = exper.pivot_samples('VALUE')

            norm = GEOparse.get_GEO(self.con)
            pheno = norm.phenotype_data[["title", "source_name_ch1"]]
            print(pheno)
            norm = norm.pivot_samples('VALUE')

            exp = [str(i) for i in input('Experiment (GSE): ').split()]
            con = [str(i) for i in input('Control (GSE): ').split()]
            nor = [str(i) for i in input('Normal (GSE): ').split()]

            exp_sample = exper[exp]
            con_sample = exper[con]
            norm_sample = norm[nor]

            dataset = pd.concat([exp_sample, con_sample, norm_sample],axis=1, sort=False)

            proc = quantile_transform(dataset, n_quantiles=1000,random_state=None)
            fold = np.log(proc)
            out_data = pd.DataFrame(fold, columns=dataset.columns,index=exper.index)

            return out_data.to_csv('dataset.csv')

#Differential Expression using Limma
class diff:
    def diff(self):
        robjects.r('data <- subset(read.csv(file="dataset.csv"),select = -ID_REF)')
        robjects.r('design <- model.matrix(~ 0+factor(c(1,1,1,2,2,2,3,3,3)))')
        robjects.r('colnames(design) <- c("Case", "Control", "Normal")')
        robjects.r('fit <- lmFit(data, design)')
        robjects.r('contrast.matrix <- makeContrasts(Control-Case,Normal-Control, Normal-Case, levels=design)')
        robjects.r('fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))')
        robjects.r('table <- topTable(fit2, coef=1, adjust="BH",number=Inf)')
        robjects.r('table <- table[order(as.numeric(rownames(table))),,drop=FALSE]')
        robjects.r('write.csv(table, file="sig_genes.csv", row.names=FALSE)')

        genes = pd.read_csv('dataset.csv')
        genes = pd.DataFrame(genes.iloc[:,0], columns=['ID_REF'])

        df = pd.read_csv('sig_genes.csv')
        df = pd.concat([genes, df], axis=1)
        df.to_csv('sig_genes.csv', index=False)

#Up-regulated genes: log2 >1 & p>0.05; Down-regulated genes: log2<1 & p>0.05
class fold:
    def fold_change(self):
        sig = pd.read_csv('sig_genes.csv')
        up = sig[(sig['logFC']>0.1) & (sig['P.Value']>0.05)]
        down = sig[(sig['logFC']<0.1) & (sig['P.Value']>0.05)]

        up.to_csv('up.csv', index=False)
        down.to_csv('down.csv', index=False)

#Annotating Probe ID list into Gene Symbol
class annotation:
    def anno(self):
        my_file = Path("/Users/lm/documents/datasets/project/diff/file.csv")

        try:
            my_abs_path = my_file.resolve(strict=True)

        except FileNotFoundError:
            robjects.r('data <- subset(read.csv(file="dataset.csv"),select = ID_REF, row.names=NULL)')
            robjects.r('probes <- as.vector(as.matrix(data))')
            robjects.r('names <- select(hgu133plus2.db, probes, c("SYMBOL"))')
            robjects.r('write.csv(names, file="file.csv", row.names=FALSE)')

        genes = pd.read_csv('file.csv')
        df = pd.read_csv('sig_genes.csv')
        df = pd.concat([genes, df], axis=1)

#Create a single dataset of collected reistant genes from multiple datasets
class resistant_genes:
    def __init__(self, up, do, genes):
        self.d1 = pd.read_csv('/Users/lm/documents/datasets/project/resistant/res1.csv')
        self.d2 = pd.read_csv('/Users/lm/documents/datasets/project/resistant/res2.csv')
        self.d3 = pd.read_csv('/Users/lm/documents/datasets/project/resistant/res3.csv')
        self.d4 = pd.read_csv('/Users/lm/documents/datasets/project/resistant/res4.csv')
        self.d5 = pd.read_csv('/Users/lm/documents/datasets/project/resistant/res5.csv')
        self.up = up
        self.do = do
        self.genes = genes

#Select the correct column from resistant gene databases, create a single dataset of resistant genes
    def res_genes(self):
        final = pd.concat([self.d1, self.d2, self.d3, self.d4, self.d5])
        final = final.drop_duplicates(subset='Gene Name')
        final.to_csv('final_dataset.csv', index=False)

#Map up and down regulated genes
    def map(self):
        up_reg = self.genes.loc[self.genes['PROBEID'].isin(self.up['ID_REF'])]
        down_reg = self.genes.loc[self.genes['PROBEID'].isin(self.do['ID_REF'])]
        res_genes = pd.read_csv('final_dataset.csv')
        up_res = up_reg.loc[up_reg['SYMBOL'].isin(res_genes['Gene Name'])]
        up_res = up_res.drop_duplicates(subset='SYMBOL')

        #print(up_res.shape)

#Pathway analysis
#Network analysis
#Build training dataset
#Model training
#Model testing
#Model evaluation

if __name__ == '__main__':

#MCF_7_TAM datasets
    read_geo.read_files(read_geo('GSE67916', 'GSE15065'))
    diff = diff()
    diff.diff()

    fold = fold()
    fold.fold_change()

    ann = annotation()
    ann.anno()

    genes = pd.read_csv('file.csv')
    up = pd.read_csv('up.csv')
    do = pd.read_csv('down.csv')

    resistant_genes.res_genes(resistant_genes(up, do, genes))
    resistant_genes.map(resistant_genes(up, do, genes))

    #print('God is Great')
