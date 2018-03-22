library(tidyverse)
library(phangorn)
library(DECIPHER)
setwd('..')
seqtab.nochim0 <- read.csv('data1/seqtab.nochimMar2018.csv')
seqtab.nochim <- seqtab.nochim0[,-1]
rownames(seqtab.nochim) <- seqtab.nochim0[,1]

seqs <- colnames(seqtab.nochim)
names(seqs) <- seqs

print('alligning sequences')
ptm <- proc.time()
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)
proc.time() - ptm

print('initial tree builing')
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)

print('optimizing tree')
ptm <- proc.time() # slow step
fitGTRopt <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

proc.time() - ptm #~7000 seconds

print('saving tree')
library(phyloseq)
pt <- phy_tree(fitGTRopt$tree)

write.tree(pt, file="data1/phylogeny096Mar2018tre.tre")
