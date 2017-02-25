library(grocaptss)
source("common.R")

cat("GM12878 TSS HMM parsing ...\n")
bwSet.gm = mainBwSet.grocap(basePath = commonDatasetPath, cellLine="gm12878")
bwBckSet.gm = backBwSet.grocap(basePath = commonDatasetPath, cellLine="gm12878")

# norm = compute.normalization(bwSet.gm)
# save(norm, file="../data/gm12878.totals.norm.Rdata")
# norm.gm.back = compute.normalization(bwBckSet.gm)
# save(norm.gm.back, file="../data/gm12878.totals.norm.back.Rdata")

load("../data/gm12878.totals.norm.Rdata")
load("../data/gm12878.totals.norm.back.Rdata")
scale.factor = (norm$GROcap.plus - norm$GROcap.minus) / (norm.gm.back$GROcap.plus - norm.gm.back$GROcap.minus)

res = process.genome.qhmm(bwSet.gm, bwBckSet.gm, bwSet.gm$GROcap.plus$chroms,scale.factor)

cat("Saving results ...\n")
bed.gm = res$preds
write.track.bed(bed.gm, filename="hg19.gm12878.new_hmm2b.bed", "hmm2b.preds.gm12878")

#
# create special tracks for M2 'peak' areas
peak.col.plus = do.call("paste", c(as.list(col2rgb(gcol$cage)), list(sep=',')))
peak.col.minus = do.call("paste", c(as.list(col2rgb(gcol$groseq)), list(sep=',')))

peaks.gm = res$peaks
write.track.bed(peaks.gm, filename="hg19.gm12878.new_hmm2bp.bed", "hmm2bp.preds.gm12878.peaks", color.plus = peak.col.plus, color.minus = peak.col.minus)
