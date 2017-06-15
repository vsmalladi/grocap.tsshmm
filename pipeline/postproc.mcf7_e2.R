library(grocaptss)
source("common.R")

cat(" * Loading data ...\n")

bwSet.mcf7_e2 = mainBwSet.grocap(basePath = commonDatasetPath, cellLine="E2")
bwBckSet.mcf7_e2 = backBwSet.grocap(basePath = commonDatasetPath, cellLine="E2")

preds.mcf7_e2 = read.table("hg19.mcf7_e2.new_hmm2b.bed", skip=1, colClasses = c("factor", "integer", "integer", "factor", "integer", "factor", "integer", "integer", "factor"))
peaks.mcf7_e2 = read.table("hg19.mcf7_e2.new_hmm2bp.bed", skip=1)

# scale.factor = 4.335949
load("../data/mcf7_e2.totals.norm.Rdata")
load("../data/mcf7_e2.totals.norm.back.Rdata")
scale.factor = abs((norm$GROcap.plus - norm$GROcap.minus) / (norm.mcf7_e2.back$GROcap.plus - norm.mcf7_e2.back$GROcap.minus))

cat(" * Applying filter ...\n")
preds.mcf7_e2.filtered = filter.split(preds.mcf7_e2, peaks.mcf7_e2, scale.factor, 100)

cat(" * Post filter expansion ...\n")
preds.mcf7_e2.filtered.post1 = edge.expand(preds.mcf7_e2.filtered, bwSet.mcf7_e2, bwBckSet.mcf7_e2, scale.factor)

cat(" * Post expansion trim ...\n")
preds.mcf7_e2.filtered.post2b = edge.trim2(preds.mcf7_e2.filtered.post1, bwSet.mcf7_e2, bwBckSet.mcf7_e2, scale.factor)

cat(" * Post filter trim (skip expansion) ...\n")
preds.mcf7_e2.filtered.post2a = edge.trim2(preds.mcf7_e2.filtered, bwSet.mcf7_e2, bwBckSet.mcf7_e2, scale.factor)

cat(" * Saving post-processed ...\n")
write.bed(preds.mcf7_e2.filtered, "../peaks/hg19.mcf7_e2.new_hmm2b.filtered.bed")
write.bed(preds.mcf7_e2.filtered.post1, "../peaks/hg19.mcf7_e2.new_hmm2b.post1.bed")
write.bed(preds.mcf7_e2.filtered.post2b, "../peaks/hg19.mcf7_e2.new_hmm2b.post2b.bed")
write.bed(preds.mcf7_e2.filtered.post2a, "../peaks/hg19.mcf7_e2.new_hmm2b.post2a.bed")
# actually used
write.bed(preds.mcf7_e2.filtered.post2a, "../peaks/hg19.mcf7_e2.new_hmm2b.post2.bed")


cat(" * Creating pairs ...\n")
pairs.150.gm.b = create.pairs(preds.mcf7_e2.filtered.post2b, 150)
pairs.150.gm.a = create.pairs(preds.mcf7_e2.filtered.post2a, 150)

cat(" * Saving pairs ...\n")
write.bed(pairs.150.gm.b$minus, "../peaks/hg19.mcf7_e2.new_hmm2b.post2b.pair_minus.bed")
write.bed(pairs.150.gm.b$plus, "../peaks/hg19.mcf7_e2.new_hmm2b.post2b.pair_plus.bed")

write.bed(pairs.150.gm.a$minus, "../peaks/hg19.mcf7_e2.new_hmm2b.post2a.pair_minus.bed")
write.bed(pairs.150.gm.a$plus, "../peaks/hg19.mcf7_e2.new_hmm2b.post2a.pair_plus.bed")
# actually used
write.bed(pairs.150.gm.a$minus, "../peaks/hg19.mcf7_e2.new_hmm2b.post2.pair_minus.bed")
write.bed(pairs.150.gm.a$plus, "../peaks/hg19.mcf7_e2.new_hmm2b.post2.pair_plus.bed")


cat(" * Creating single ...\n")
single.gm.b = create.single(preds.gm.filtered.post2b, 250)
single.gm.a = create.single(preds.gm.filtered.post2a, 250)

cat(" * Writing single ...\n")
write.bed(single.gm.b, "../peaks/hg19.gm12878.new_hmm2b.post2b.single.bed")
write.bed(single.gm.a, "../peaks/hg19.gm12878.new_hmm2b.post2a.single.bed")
# actually used
write.bed(single.gm.a, "../peaks/hg19.gm12878.new_hmm2b.post2.single.bed")
