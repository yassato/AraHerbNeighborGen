#########################################
# selection scan by EHH for A. thaliana #
#########################################

# set target chr. no.
tar = 5

# load annotation
# Pos: position at bp
# Alyr: A. lyrata allele
# Atha: A. thaliana allele
# Score: alignment score from the original .maf
chr = read.csv("../data/chr1_ann.csv",header=TRUE)
chr = chr[,-1]
Chr = rep(tar,nrow(chr))
chr = cbind(Chr,chr)
chr = cbind(c(1:nrow(chr)),chr)
colnames(chr) = c("ID", "chr", "pos", "Alyr", "Atha", "Score")

#load genotypes
geno = readRDS("../data/sub_snpMAF5.rds")
position = readRDS("../data/positionsMAF5.rds")
position = as.data.frame(position)

geno = geno[position$chr==tar,]
position = position[position$chr==tar,]

TF = (paste0(chr$chr,"-",chr$pos) %in%  paste0(position$chr,"-",position$position))
id = (1:nrow(chr))[TF]
map = chr[TF,]

start = proc.time()

hap_i = function(i) {
  search_row = which(position$position==chr[i,]$pos)
  if(length(search_row)!=0) { print(i)
    geno_i = as.numeric(geno[search_row,])
    
    if(as.character(chr[i,]$Alyr)==as.character(chr[i,]$Atha)) {
      geno_i[geno_i==1] = 2
      geno_i[geno_i==0] = 1
    } else {
      geno_i[geno_i==0] = 2
      geno_i[geno_i==1] = 1
    }
  }
  return(geno_i)
}

hap = mapply(hap_i,id)
hap = cbind(c(1:ncol(geno)), hap)

end = proc.time() - start

# save intermediate input files
saveRDS(map, file=paste0("../output/chr",tar,"_AthaMapInfo.rds"), compress=TRUE, version=2)
saveRDS(hap, file=paste0("../output/chr",tar,"_AthaHaploInfo.rds"), compress=TRUE, version=2)

write.table(map, file=paste0("../output/map",tar,".inp"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep=" ")
write.table(hap, file=paste0("../output/chr",tar,".hap"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep=" ")

# load library and calc. EHH and then iHS
# note: scan_hh() is not enough. iHS is required to compare two alleles.
library(rehh)
hap_data=data2haplohh(hap_file=paste0("../output/chr",tar,".hap"), 
                      map_file=paste0("../output/map",tar,".inp"),
                      recode.allele=FALSE, chr.name=tar)

res_hh = scan_hh(hap_data, threads=1)
res_ihs = ihh2ihs(res_hh)
res_ihs = res_ihs$ihs

write.csv(res_hh, file=paste0("../output/scan_hh200_chr",tar,".csv"))
write.csv(res_ihs, file=paste0("../output/scan_ihs200_chr",tar,".csv"))


