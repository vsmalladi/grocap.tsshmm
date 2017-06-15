library(grocaptss)
source("common.R")

#cat("GM12878 TSS HMM parsing ...\n")
#bwSet.gm = mainBwSet.grocap(basePath = commonDatasetPath, cellLine="gm12878")
#bwBckSet.gm = backBwSet.grocap(basePath = commonDatasetPath, cellLine="gm12878")


cat("MCF-7 TSS HMM parsing ...\n")
bwSet.mcf7 = mainBwSet.grocap(basePath = commonDatasetPath, cellLine="VEH")
bwBckSet.mcf7 = backBwSet.grocap(basePath = commonDatasetPath, cellLine="VEH")

norm = compute.normalization(bwSet.mcf7)
save(norm, file="../data/mcf7.totals.norm.Rdata")
norm.mcf7.back = compute.normalization(bwBckSet.mcf7)
save(norm.mcf7.back, file="../data/mcf7.totals.norm.back.Rdata")

load("../data/mcf7.totals.norm.Rdata")
load("../data/mcf7.totals.norm.back.Rdata")
scale.factor = abs((norm$GROcap.plus - norm$GROcap.minus) / (norm.mcf7.back$GROcap.plus - norm.mcf7.back$GROcap.minus))

res = process.genome.qhmm(bwSet.mcf7, bwBckSet.mcf7, bwSet.mcf7$GROcap.plus$chroms, scale.factor)

cat("Saving results ...\n")
bed.mcf7 = res$preds
write.track.bed(bed.mcf7, filename="../peaks/hg19.mcf7.new_hmm2b.bed", "hmm2b.preds.mcf7")

#
# create special tracks for M2 'peak' areas
peak.col.plus = do.call("paste", c(as.list(col2rgb(gcol$cage)), list(sep=',')))
peak.col.minus = do.call("paste", c(as.list(col2rgb(gcol$groseq)), list(sep=',')))

peaks.mcf7 = res$peaks
write.track.bed(peaks.mcf7, filename="../peaks/hg19.mcf7.new_hmm2bp.bed", "hmm2bp.preds.mcf7.peaks", color.plus = peak.col.plus, color.minus = peak.col.minus)
