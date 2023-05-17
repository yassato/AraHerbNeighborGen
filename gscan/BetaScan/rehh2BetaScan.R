#################################################
# convert rehh inputs as an input for Beta scan #
#################################################

# set target chr. no.
tar = 5

hap = readRDS(file=paste0("../output/chr",tar,"_AthaHaploInfo.rds"))
map = readRDS(file=paste0("../output/chr",tar,"_AthaMapInfo.rds"))

hap2 = hap[,-1]
hap2 = hap2 - 1

num_derived = as.numeric(apply(hap2,2,sum))
num_hap = rep(nrow(hap2),ncol(hap2))

out = cbind(map$pos,num_derived,num_hap)
write.table(out, file=paste0("../output/BetaScanChr",tar,"in.txt"),
            quote=FALSE, col.names=FALSE, row.names=FALSE, sep=" ")

