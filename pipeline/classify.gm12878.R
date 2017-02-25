library(grocaptss)
source("common.R")

peaks.plus = read.table("hg19.gm12878.new_hmm2b.post2.pair_plus.bed")
peaks.minus = read.table("hg19.gm12878.new_hmm2b.post2.pair_minus.bed")

# remove chrY
peaks.plus = peaks.plus[peaks.plus[,1] != "chrY", ]
peaks.minus = peaks.minus[peaks.minus[,1] != "chrY", ]

bwSet = mainBwSet.grocap(basePath = commonDatasetPath, cellLine="gm12878")

res = stable.unstable.classify(peaks.plus, peaks.minus, bwSet)

stable.unstable.write("hg19.gm12878.new_hmm2b.post2", res)

#
# also create singleton set using thresholds from stable/unstable classification
#

single.gm = read.table("hg19.gm12878.new_hmm2b.post2.single.bed")
res.gm = single.filter(single.gm, bwSet, res$thresh.plus, res$thresh.minus)

write.table(res.gm$all, file="hg19.gm12878.new_hmm2b.post2.singletons.bed", sep='\t', quote=F, col.names=F, row.names=F)
write.table(res.gm$stable, file="hg19.gm12878.new_hmm2b.post2.singletons.S.bed", sep='\t', quote=F, col.names=F, row.names=F)
write.table(res.gm$unstable, file="hg19.gm12878.new_hmm2b.post2.singletons.U.bed", sep='\t', quote=F, col.names=F, row.names=F)
