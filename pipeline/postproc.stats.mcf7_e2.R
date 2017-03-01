#
library(grocaptss)
source("common.R")

# load data
bwSet.mcf7_e2 = mainBwSet.grocap(basePath = commonDatasetPath, cellLine="E2")
bwBckSet.mcf7_e2 = backBwSet.grocap(basePath = commonDatasetPath, cellLine="E2")

#
load("../data/mcf7_e2.totals.norm.Rdata")
load("../data/mcf7_e2.totals.norm.back.Rdata")
scale.factor = abs((norm$GROcap.plus - norm$GROcap.minus) / (norm.gm.back$GROcap.plus - norm.gm.back$GROcap.minus))


# load predictions
preds.hmm = read.table("hg19.mcf7_e2.new_hmm2b.bed", skip = 1)
preds.post2 = read.table("hg19.mcf7_e2.new_hmm2b.post2.bed")


# step type coverage
#
n.cores = 1
chroms = levels(preds.hmm[,1])
stepcov.preds.hmm = step.coverage(preds.hmm, chroms, bwSet.mcf7_e2, bwBckSet.mcf7_e2, scale.factor, n.cores = n.cores)
stepcov.preds.post2 = step.coverage(preds.post2, chroms, bwSet.mcf7_e2, bwBckSet.mcf7_e2, scale.factor, n.cores = n.cores)

# plot step coverage information
pdf("sup.s1.panelC.new.pdf", width=14, height = 3)
opar = par(no.readonly=TRUE)
par(mar=c(3, 4, 2, 1) + 0.1, font = 2, font.lab = 2)
par(mfrow=c(1,2))

step.data = rbind(
  stepcov.preds.hmm$covered / stepcov.preds.hmm$totals,
  stepcov.preds.post2$covered / stepcov.preds.post2$totals)

barplot(step.data, beside = T, ylim=c(0,1), names.arg = c("TAP+ > 0", "TAP+ > TAP-\n(TAP- > 0)", "TAP+ > TAP-\n(TAP- = 0)", "TAP+ < TAP-"),
        col=c(gcol$fwd, gcol$rev),
        legend.text = c("Viterbi", "Post-proc. "),
        ylab = "Fraction of steps covered",
        border = NA,
        space = c(0.1, 0.1, 0.5, 0.1, 0.5, 0.1, 0.5, 0.1),
        font = 2,
        lwd = 2,
        args.legend = list(x = "topright", box.lwd=2, border=NA, title = "PRO-cap HMM"))

# genome coverate
#
genome.size = 2897316137
chrom.hmm.gm = chromHMM.regions(basePath = commonDatasetPath, cellLine="gm12878")

gcov.chrom.hmm = bitmask.total(rbind(chrom.hmm.gm$proms, chrom.hmm.gm$enh))
gcov.preds.hmm = bitmask.total(preds.hmm)
gcov.preds.post2 = bitmask.total(preds.post2)

gcov.frac = c(gcov.chrom.hmm, gcov.preds.hmm, gcov.preds.post2)/genome.size

# plot

barplot(gcov.frac, beside = T, names.arg = c("ChromHMM\n(Prom.,Enh.)", "GRO-cap HMM\nViterbi", "GRO-cap HMM\nPost-proc."), col=c(gcol$cage, gcol$fwd, gcol$rev),
        ylab = "Fraction of genome covered",
        border = NA,
        lwd = 2,
        space = 0.08,
        font = 2)

par(mfrow=c(1,1))
par(opar)

dev.off()



# library coverage
#
n.cores = 1
chroms = levels(preds.hmm[,1])
readcov.preds.hmm = read.coverage(preds.hmm, chroms, bwSet.mcf7_e2, scale.factor, n.cores = n.cores)
readcov.preds.post2 = read.coverage(preds.post2, chroms, bwSet.mcf7_e2, scale.factor, n.cores = n.cores)

read.cov = c(
  readcov.preds.hmm$covered / readcov.preds.hmm$totals,
  readcov.preds.post2$covered / readcov.preds.post2$totals)

#
# load pairs
pairs.plus = read.table("hg19.mcf7_e2.new_hmm2b.post2.pair_plus.bed")
pairs.minus = read.table("hg19.mcf7_e2.new_hmm2b.post2.pair_minus.bed")


readcov.pairs = read.coverage(rbind(pairs.plus, pairs.minus), chroms, bwSet.gm, scale.factor, n.cores = n.cores)


#
# chromHMM prom/enh fraction breakdown
proms.preds.hmm = contained.fraction(container = chrom.hmm.gm$proms, regions = preds.hmm)
proms.preds.post2 = contained.fraction(container = chrom.hmm.gm$proms, regions = preds.post2)
proms.pairs = pairs.in.regions(pairs.plus, pairs.minus, container = chrom.hmm.gm$proms)

enh.preds.hmm = contained.fraction(container = chrom.hmm.gm$enh, regions = preds.hmm)
enh.preds.post2 = contained.fraction(container = chrom.hmm.gm$enh, regions = preds.post2)
enh.pairs = pairs.in.regions(pairs.plus, pairs.minus, container = chrom.hmm.gm$enh)



# chromHMM inclusion
#
active.chrom.hmm.gm = chromHMM.regions(basePath = commonDatasetPath, cellLine="gm12878", prom.states = "1", enh.states = "45")
act.proms.hmm = overlaped.mask(active.chrom.hmm.gm$proms, preds.hmm)
act.proms.post2 = overlaped.mask(active.chrom.hmm.gm$proms, preds.post2)

# size distribution
#
sizes.hmm = preds.hmm[,3] - preds.hmm[,2]
sizes.post2 = preds.post2[,3] - preds.post2[,2]

pdf("../output/sup.s1.panelB.new.pdf", width=14, height = 3)
opar = par(no.readonly=TRUE)
par(mfrow=c(2,1), mar = c(3, 1, 2, 1), font = 2, font.lab = 2, font.axis = 2, lwd = 2)
boxplot(sizes.hmm, ylim=c(0, 1600), main="Viterbi", horizontal = TRUE)
boxplot(sizes.post2, ylim=c(0, 1600), main = "Post-processed", horizontal = TRUE)
mtext("TSS region length (bp)", side = 1, line = 1.5)
par(mfrow=c(1,1))
par(opar)
dev.off()
