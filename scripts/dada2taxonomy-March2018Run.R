library('dada2'); packageVersion("dada2")

setwd('..')

seqtab.nochim0 <- read.csv("data1/seqtab.nochimMar2018.csv", row.names = 1)
seqtab.nochim1 <- as.matrix(seqtab.nochim0)

## Part 2: Assign Taxonomy
set.seed(33)
seqs <- getSequences(seqtab.nochim1)
seqsRC <- dada2:::rc(seqs)
seqsRChead <- seqsRC[1:10]

ptm <- proc.time()
taxa <- assignTaxonomy(seqsRC, "Training/rdp_train_set_16.fa.gz", multithread=TRUE, tryRC=TRUE, minBoot=80)
proc.time() - ptm

## The following also works
# taxa <- assignTaxonomy(seqtab.nochim1, "Training/rdp_train_set_16.fa.gz", multithread=TRUE, tryRC=TRUE, minBoot=80)


# reverse complement because thats how this study went

ptm <- proc.time()
genus.species.rc <- assignSpecies(seqsRC, "Training/rdp_species_assignment_16.fa.gz", allowMultiple = TRUE)
proc.time() - ptm


 
taxaDf <- data.frame(taxa)
taxaDf$seq <- rownames(taxaDf)

gsDf <- data.frame(genus.species.rc)
gsDf$seq <- dada2:::rc(rownames(gsDf))
gsDf$seq <- rownames(gsDf)

library(dplyr); packageVersion('dplyr')
taxaGS <- left_join(taxaDf, gsDf, by = c('seq'))

# make sure genera all have same names
gencheck <- na.omit(taxaGS[,c("Genus.x", "Genus.y", 'Species')])

gencheck[which(as.character(gencheck$Genus.x) != as.character(gencheck$Genus.y)),]

# some of the Clostridum genera have fancy genus names that don't show up in the assign sepcies thing


colnames(taxaGS)

taxaGS2 <- taxaGS[c("seq", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus.x", "Genus.y", "Species")] %>% rename(Genus = Genus.x) %>% rename('Sequence' = 'seq')

colnames(taxaGS2)

write.csv(taxaGS2, 'data1/TaxaMar2018.csv', row.names = FALSE)
