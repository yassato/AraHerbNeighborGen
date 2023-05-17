##########################
# GO enrichment analysis #
##########################

# load library and functions
library(GO.db)

# install this utility package, "GOfisher", beforehand
# devtools::install_github("https://github.com/yassato/GOfisher") 
library(GOfisher) 

library(rrvgo)
library(org.At.tair.db)

# load a list of AGI and GO SLIM. 
load("../AthGO/ulgSLIM.TAIR_211201")

rrvgo_fig = function(gl) {
  fisher.res = ng.mft(ulg, gl)
  res = ng.GOfisher2REVIGO(fisher.res[fisher.res[,"xtt"]>1,],gl,ulg)
  res = na.omit(res)
  res = res[res$qvalue<0.05,]
  simMatrix = calculateSimMatrix(res$ID, orgdb="org.At.tair.db",ont="BP", method="Rel")
  scores = setNames(-log10(res$qvalue), res$ID)
  reducedTerms = reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.At.tair.db")
  treemapPlot(reducedTerms)
  # scatterPlot(simMatrix, reducedTerms, labelSize=2)
  heatmapPlot(simMatrix,
              reducedTerms,
              annotateParent=TRUE,
              annotationLabel="parentTerm",
              fontsize=4)
  return(reducedTerms)
}

genes = read.csv("../output/HolesCHZ_glmnetLassoMAF5_mean_10kbAGI.csv")
gl_pos = genes[genes$nei_beta>0,"Locus"]
gl_pos = unique(gl_pos)
gl_neg = genes[genes$nei_beta<0,"Locus"]
gl_neg = unique(gl_neg)

res1 = rrvgo_fig(gl_pos)
res2 = rrvgo_fig(gl_neg)

# write REVIGO tables
write.csv(res1,file="../output/REVIGO_LASSO_pos2.csv")
write.csv(res2,file="../output/REVIGO_LASSO_neg2.csv")

# search particular GOs from a gene list
#GO:0006952, defense response; GO:0009607, response to biotic stimulus; 
#GO:0009611, response to wounding; GO:0002213, defense response to insect; GO:0009625, response to insect
#GO:0019761, glucosinolate biosynthetic process; GO:0019748, secondary metabolic process
#GO:0009695, jasmonic acid biosynthetic process, "GO:0009753", response to jasmonic acids
#GO:0010026, trichome differentiation; GO:0048629, trichome patterning
target = "GO:0009695"
withgo = ulg[ulg[,"GOid"]==target, "locus"]
des[intersect(gl_pos, withgo),]
# ComputationalDescription   GeneSymbol NormalizationGroup
# AT1G67560 PLAT/LH2 domain-containing lipoxygenase family protein;(source:Araport11) ATLOX6, LOX6               data
# AT3G45140                                         lipoxygenase 2;(source:Araport11) ATLOX2, LOX2               data
