library(grocaptss)
source("common.R")

#cat("GM12878 TSS HMM parsing ...\n")
#bwSet.gm = mainBwSet.grocap(basePath = commonDatasetPath, cellLine="gm12878")
#bwBckSet.gm = backBwSet.grocap(basePath = commonDatasetPath, cellLine="gm12878")


cat("MCF-7 ES TSS HMM parsing ...\n")
bwSet.mcf7_e2 = mainBwSet.grocap(basePath = commonDatasetPath, cellLine="E2")
bwBckSet.mcf7_e2 = backBwSet.grocap(basePath = commonDatasetPath, cellLine="E2")

norm = compute.normalization(bwSet.mcf7_e2)
save(norm, file="../data/mcf7_e2.totals.norm.Rdata")
norm.mcf7_e2.back = compute.normalization(bwBckSet.mcf7_e2)
save(norm.mcf7_e2.back, file="../data/mcf7_e2.totals.norm.back.Rdata")

load("../data/mcf7_e2.totals.norm.Rdata")
load("../data/mcf7_e2.totals.norm.back.Rdata")
scale.factor = abs((norm$GROcap.plus - norm$GROcap.minus) / (norm.mcf7_e2.back$GROcap.plus - norm.mcf7_e2.back$GROcap.minus))

res = process.genome.qhmm(bwSet.mcf7_e2, bwBckSet.mcf7_e2, bwSet.mcf7_e2$GROcap.plus$chroms, scale.factor)

cat("Saving results ...\n")
bed.mcf7_e2 = res$preds
write.track.bed(bed.mcf7_e2, filename="../peaks/hg19.mcf7_e2.new_hmm2b.bed", "hmm2b.preds.mcf7_e2")

#
# create special tracks for M2 'peak' areas
peak.col.plus = do.call("paste", c(as.list(col2rgb(gcol$cage)), list(sep=',')))
peak.col.minus = do.call("paste", c(as.list(col2rgb(gcol$groseq)), list(sep=',')))

peaks.mcf7_e2 = res$peaks
write.track.bed(peaks.mcf7_e2, filename="../peaks/hg19.mcf7_e2.new_hmm2bp.bed", "hmm2bp.preds.mcf7_e2.peaks", color.plus = peak.col.plus, color.minus = peak.col.minus)
