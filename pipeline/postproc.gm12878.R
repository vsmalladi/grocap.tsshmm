library(grocaptss)
source("common.R")

cat(" * Loading data ...\n")

bwSet.gm = mainBwSet.grocap(basePath = commonDatasetPath, cellLine="gm12878")
bwBckSet.gm = backBwSet.grocap(basePath = commonDatasetPath, cellLine="gm12878")

preds.gm = read.table("hg19.gm12878.new_hmm2b.bed", skip=1, colClasses = c("factor", "integer", "integer", "factor", "integer", "factor", "integer", "integer", "factor"))
peaks.gm = read.table("hg19.gm12878.new_hmm2bp.bed", skip=1)

# scale.factor = 4.335949
load("../data/gm12878.totals.norm.Rdata")
load("../data/gm12878.totals.norm.back.Rdata")
scale.factor = (norm$GROcap.plus - norm$GROcap.minus) / (norm.gm.back$GROcap.plus - norm.gm.back$GROcap.minus)

cat(" * Applying filter ...\n")
preds.gm.filtered = filter.split(preds.gm, peaks.gm, scale.factor, 100)

cat(" * Post filter expansion ...\n")
preds.gm.filtered.post1 = edge.expand(preds.gm.filtered, bwSet.gm, bwBckSet.gm, scale.factor)

cat(" * Post expansion trim ...\n")
preds.gm.filtered.post2b = edge.trim2(preds.gm.filtered.post1, bwSet.gm, bwBckSet.gm, scale.factor)

cat(" * Post filter trim (skip expansion) ...\n")
preds.gm.filtered.post2a = edge.trim2(preds.gm.filtered, bwSet.gm, bwBckSet.gm, scale.factor)

cat(" * Saving post-processed ...\n")
write.bed(preds.gm.filtered, "hg19.gm12878.new_hmm2b.filtered.bed")
write.bed(preds.gm.filtered.post1, "hg19.gm12878.new_hmm2b.post1.bed")
write.bed(preds.gm.filtered.post2b, "hg19.gm12878.new_hmm2b.post2b.bed")
write.bed(preds.gm.filtered.post2a, "hg19.gm12878.new_hmm2b.post2a.bed")
# actually used
write.bed(preds.gm.filtered.post2a, "hg19.gm12878.new_hmm2b.post2.bed")


cat(" * Creating pairs ...\n")
pairs.150.gm.b = create.pairs(preds.gm.filtered.post2b, 150)
pairs.150.gm.a = create.pairs(preds.gm.filtered.post2a, 150)

cat(" * Saving pairs ...\n")
write.bed(pairs.150.gm.b$minus, "hg19.gm12878.new_hmm2b.post2b.pair_minus.bed")
write.bed(pairs.150.gm.b$plus, "hg19.gm12878.new_hmm2b.post2b.pair_plus.bed")

write.bed(pairs.150.gm.a$minus, "hg19.gm12878.new_hmm2b.post2a.pair_minus.bed")
write.bed(pairs.150.gm.a$plus, "hg19.gm12878.new_hmm2b.post2a.pair_plus.bed")
# actually used
write.bed(pairs.150.gm.a$minus, "hg19.gm12878.new_hmm2b.post2.pair_minus.bed")
write.bed(pairs.150.gm.a$plus, "hg19.gm12878.new_hmm2b.post2.pair_plus.bed")


cat(" * Creating single ...\n")
single.gm.b = create.single(preds.gm.filtered.post2b, 250)
single.gm.a = create.single(preds.gm.filtered.post2a, 250)

cat(" * Writing single ...\n")
write.bed(single.gm.b, "hg19.gm12878.new_hmm2b.post2b.single.bed")
write.bed(single.gm.a, "hg19.gm12878.new_hmm2b.post2a.single.bed")
# actually used
write.bed(single.gm.a, "hg19.gm12878.new_hmm2b.post2.single.bed")
