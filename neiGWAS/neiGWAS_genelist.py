# install python module using reticulate (if not installed)
# reticulate::py_install("pandas")

import pandas as pd
import numpy as np
import re
gff = pd.read_csv("../data/TAIR10_GFF3_genes.gff.gz", sep="\t", header=None)
gff.columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff_gene = gff[(gff.type=="gene")|(gff.type=="transposable_element_gene")|(gff.type=="pseudogene")]
des = pd.read_csv("../AthDescription/Araport11_genes.201606.transcript.rep_ERCC_Virus7457_GFP_GUS.txt.gz", sep="\t", compression="gzip")
    
# change "P_nei" into "P_self" when compiling self-genotype effects
def SNP2GENES_lmm(Chr, pos, window):
    start = int(pos) - int(window)
    end = int(pos) + int(window)
    seqID = "Chr" + str(Chr)
    gff_sub = gff_gene[gff_gene.seqid==str(seqID)]
    target = gff_sub[(gff_sub.start>start)&(gff_sub.end<end)]
    a = "AT"
    b = str(Chr)
    c = "G"
    pattern = a + b + c
    text = target.loc[:,"attributes"]
    res_list = pd.DataFrame(index=[], columns=[])
    for i in text:
        matchOB = re.search(pattern, i)
        atg = i[matchOB.start():(matchOB.end()+5)]
        desRes = des[des.GeneID==atg]
        nrows = len(desRes.index)
        res = pd.DataFrame(index=[], columns=[])
        res["Chr"] = np.repeat(Chr, nrows)
        res["Pos"] = np.repeat(pos, nrows)
        res["MAF"] = np.repeat(float(gwas_out[(gwas_out.Chr==int(Chr))&(gwas_out.Position==int(pos))].MAF), nrows)
        res["beta"] = np.repeat(float(gwas_out[(gwas_out.Chr==int(Chr))&(gwas_out.Position==int(pos))].beta_nei), nrows)
        res["P"] = np.repeat(float(gwas_out[(gwas_out.Chr==int(Chr))&(gwas_out.Position==int(pos))].P_nei), nrows)
        res["Locus"] = np.array(desRes.GeneID)
        res["Type"] = np.array(desRes.GeneModelType)
        res["Symbol"] = np.array(desRes.GeneSymbol)
        res["ShortDescription"] = np.array(desRes.ShortDescription)
        res["CuratorSummary"] = np.array(desRes.CuratorSummary)
        res_list = res_list.append(res)
    return res_list

f_path = "../output/"
f_name = "JPNneiGWAS_richnessS1"
f_input = f_path + f_name + ".csv"
gwas_out = pd.read_csv(f_input)
padj = np.percentile(gwas_out["P_nei"].dropna(),0.1)
gwas_out_p = gwas_out[gwas_out["P_nei"]<padj]

print(padj)

gene_list = pd.DataFrame(index=[], columns=[])
for j in range(0, gwas_out_p.shape[0]):
    res = SNP2GENES_lmm(gwas_out_p.iat[j,1], gwas_out_p.iat[j,2], 10000) # window = 10 kb
    gene_list = gene_list.append(res)

f_out = "../geneList/" + f_name + "_0001_P_nei_10kb.txt" # export the results outside the GitHub repository
gene_list.to_csv(f_out,sep="\t")


