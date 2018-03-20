## # dada2 processing of 096 pilot study
## following tutorial at
## https://benjjneb.github.io/dada2/tutorial.html

## Two changes had to happen from last time
## https://benjjneb.github.io/dada2/faq.html#can-i-use-dada2-with-my-454-or-ion-torrent-data

## filterAndTrim(..., maxLen=XXX) # XXX depends on the chemistry # lets set to 500, well bring trunclen down to 350.
## https://github.com/benjjneb/dada2/issues/275
## dada(..., HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)

library('dada2'); packageVersion("dada2")
library('ggplot2'); packageVersion('ggplot2')

setwd('..')

path <- "for_dada2"
fnFs00 <- list.files(path)
fnFs0 <- fnFs00[grep("fastq", list.files(path))]
sample.names <- sapply(fnFs0, function(file) gsub('\\.fasta',"", file))
fnFs <- file.path(path, fnFs0)

filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))

print('set seed to 33')
set.seed(33)

print('filter and trim sequences')
out <- filterAndTrim(fnFs, filtFs, trimLeft = 22, truncLen=300, maxLen=500, maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                 compress=TRUE, multithread=TRUE)

print('learing errors')
ptm = proc.time()
errF <- learnErrors(filtFs, multithread=TRUE)
proc.time() - ptm

save(errF, file = 'proc096/errF-Mar2018.Rdata')

pdf('figures/errF-Mar2018.pdf') # nreads shoudn't matter, there are only 306339 of them
plotErrors(errF, nominalQ=TRUE)
dev.off()

print('dereplicating sequences')
derepFs <- derepFastq(filtFs, verbose=TRUE)
names(derepFs) <- sample.names

print('run dada function')
ptm <- proc.time()
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32, pool = TRUE)
proc.time() - ptm

seqtab <- makeSequenceTable(dadaFs)

print('remove chimeras')
ptm <- proc.time()
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
proc.time() - ptm # 3 seconds

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "tabled", "nonchim")
rownames(track) <- sample.names
#track[10:20,]
track

write.csv(track, 'proc096/trackMar2018.csv')

write.csv(seqtab.nochim, "data1/seqtab.nochimMar2018.csv")

