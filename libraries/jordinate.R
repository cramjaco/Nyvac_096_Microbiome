library(phyloseq)
library(vegan)
library(dplyr)
## Load in test data, comment out if this becomes a library.
#load("pspg4.Rdata")

## ## Testing commented out
## ## Building an ordination function for capscale objects.
## data("GlobalPatterns")

## ## Process and downsample as in https://joey711.github.io/phyloseq/plot_ordination-examples.html
## GP = GlobalPatterns
## wh0 = genefilter_sample(GP, filterfun_sample(function(x) x > 5), A=0.5*nsamples(GP))
## GP1 = prune_taxa(wh0, GP)
## GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))
## phylum.sum = tapply(taxa_sums(GP1), tax_table(GP1)[, "Phylum"], sum, na.rm=TRUE)
## top5phyla = names(sort(phylum.sum, TRUE))[1:5]
## GP1 = prune_taxa((tax_table(GP1)[, "Phylum"] %in% top5phyla), GP1)

## GP1Log <- transform_sample_counts(GP1, function(x) log(1+x))
## GP1Cap <- ordinate(GP1Log, method = "CAP", formula = ~SampleType)
## GP1Scores <- vegan::scores(GP1Cap, display = c("sp", "wa", "bp"), choices = c(1,2))
## ## Process input data

## 

##
## Process Taxa
proc_taxa <- function(physeq, ordscores, CAPNames = c("CAP1", "CAP2")){
    specdf <- merge(ordscores$species, tax_table(physeq), all.x = TRUE, by = 0)
    specdf <- mutate(specdf, DO = sqrt(get(CAPNames[1])^2 + get(CAPNames[2])^2))
    ## specdfTop <- top_n(specdf, 5, DO)
    specdfReps <- unique(rbind(
        top_n(specdf, 2, get(CAPNames[1])),
        top_n(specdf, 2, get(CAPNames[2])),
        top_n(specdf, -2, get(CAPNames[1])),
        top_n(specdf, -2, get(CAPNames[2]))
        ))
    return(list(specdf = specdf, specdfReps = specdfReps))
}

## test <- proc_taxa(GP1Log, GP1Scores, CAPNames = c("CAP1", "CAP2"))

## Process Sites
proc_site <- function(physeq, ordscores, CAPNames = c("CAP1", "CAP2")){
    sitedf <- merge(ordscores$sites, sample_data(physeq), all.x = TRUE, by = 0)
    sitedf <- mutate(sitedf, DO = sqrt(get(CAPNames[1])^2 + get(CAPNames[2])^2))
    return(list(sitedf = sitedf))
}

## test <- proc_site(GP1Log, GP1Scores)

## Process biplot vectors
proc_biplot <- function(physeq, ordscores, CAPNames = c("CAP1", "CAP2")){
    biplotdf <- data.frame(ordscores$biplot)
    biplotdf <- mutate(biplotdf, DO = sqrt(get(CAPNames[1])^2 + get(CAPNames[2])^2), Row.names = rownames(biplotdf))
    head(biplotdf)
    biplotdfReps <- unique(rbind(
        top_n(biplotdf, 2, get(CAPNames[1])),
        top_n(biplotdf, 2, get(CAPNames[2])),
        top_n(biplotdf, -2, get(CAPNames[1])),
        top_n(biplotdf, -2, get(CAPNames[2]))
    ))
    list(biplotdf = biplotdf, biplotdfReps = biplotdfReps)
}

proc_all <- function(physeq, ordscores, CAPNames = c("CAP1", "CAP2")){
    
    pt <- proc_taxa(physeq,ordscores, CAPNames)
    ps <- proc_site(physeq, ordscores, CAPNames)
    bp <- proc_biplot(physeq, ordscores, CAPNames)
    c(pt, ps, bp)
    
}

## test <- proc_all(GP1Log, GP1Scores)

## Setup
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",  "#0072B2", "#D55E00", "#CC79A7", "#F0E442")
cbPalette1 <- c("#E69F00", "#56B4E9", "#009E73",  "#0072B2", "#D55E00", "#CC79A7", "#F0E442")

### Plot
triplot_cap <- function(physeq, ord, choices = c(1,2), scaling = 3, sitefill = NULL, ...){
    ## Necessary libraries
    require(vegan)
    require(ggplot2)
    require(viridis)
    require(dplyr)
    ## Check input functions
    if(!is(physeq, "phyloseq")){
        stop("physeq must be of class 'phyloseq'")
    }
    if(!is(ord, 'capscale')){
        stop("ord must be of class 'capscale'")
    }

    ## Names of capscale choices
    CAPaName <- paste("CAP", choices[1], sep = "")
    CAPbName <- paste("CAP", choices[2], sep = "")
    
    ## Some math
    ord.evals <- ord$CCA$eig
    sco <- vegan::scores(ord, display = c("sp", "wa", "bp"), choices = choices,
                  scaling = scaling)
    proc <- proc_all(physeq, sco, CAPNames = c(CAPaName, CAPbName))

 

    ## Axis labels
    xlab <- paste(CAPaName, sprintf("[%1.1f%%]", 100 * summary(ord)$cont$importance[2,choices[1]]))
    xlab
    ylab <- paste(CAPbName, sprintf("[%1.1f%%]", 100 * summary(ord)$cont$importance[2,choices[2]]))
    ylab
    
    ## Plot body
    p <- ggplot() +
                                        # Sites (Participants)
        geom_point(data = proc$sitedf, aes_string(
                                           x = CAPaName,
                                           y = CAPbName,
                                           fill = sitefill
                                       ), pch = 21, size = 3) +
                                        # Taxa
        geom_point(data = proc$specdf, aes_string(
                                           x = CAPaName,
                                           y = CAPbName,
                                           color = 'Phylum'
                                       ), pch = 2)+
        geom_text(data = proc$specdfReps, aes(
                                              x = get(CAPaName) * 1.05,
                                              y = get(CAPbName) * 1.05,
                                              color = Phylum,
                                              label = Genus
                                          ), pch = 2, cex = 3) +
        scale_colour_manual(values=cbPalette1) +
        scale_fill_viridis(direction = 1) +
        
                                        # Biplot Arrows; Antigens
        geom_segment(data = proc$biplotdf, aes(
                                               x = 0, y = 0,
                                               xend = get(CAPaName), 
                                               yend = get(CAPbName)
                                           ),
                     size = 0.5,
                     arrow = arrow(length = unit(0.02, "npc"))
                     ) +
        geom_text(data = proc$biplotdfReps, aes(
                                           x = get(CAPaName) * 1.05,
                                           y = get(CAPbName) * 1.05,
                                           label = Row.names
                                           ), cex = 4)+
                                        # Adjust Axes
        coord_fixed(sqrt(ord.evals[choices[2]] / ord.evals[choices[1]]))+
        labs(x = xlab, y = ylab)
        
    p
}

##
#triplot_cap(GP1Log, GP1Cap)
## triplot_cap(GP1Log, GP1Cap, sitefill = "DO", scaling = 1, choices = c(2,3))
#triplot_cap(1, GP1Cap)
    

    

log_ratio_x <- function(x){
    loc_ncol = dim(x)[2]
    loc_nrow = dim(x)[1]
    nout = (loc_ncol * (loc_ncol - 1))/2
    outmat = matrix(, nrow = loc_nrow, ncol = nout)
    to_colnames = matrix(,nrow = 0, ncol = nout)
    x_colnames = colnames(x)
    k = 0
    for(i in 1:(loc_ncol-1)){
        i
        for(j in (i+1):loc_ncol){
            k = k + 1
            loc_num = x[,i]
            loc_denom = x[,j]
            loc_lr = log(loc_num/loc_denom)
            outmat[,k] = loc_lr
            to_colnames[k] = paste(x_colnames[i], x_colnames[j], sep = "-")
        }
        }
    colnames(outmat) = to_colnames
    outmat
    }

log_ratio_x2 <- function(x){
    loc_ncol = dim(x)[2]
    loc_nrow = dim(x)[1]
    nout = (loc_ncol * (loc_ncol - 1))
    outmat = matrix(, nrow = loc_nrow, ncol = nout)
    to_colnames = matrix(,nrow = 0, ncol = nout)
    x_colnames = colnames(x)
    k = 0
    for(i in 1:loc_ncol){
        for(j in 1:loc_ncol){
            if(i !=j){
                k = k + 1
                loc_num = x[,i]
                loc_denom = x[,j]
                loc_lr = log(loc_num/loc_denom)
                outmat[,k] = loc_lr
                to_colnames[k] = paste(x_colnames[i], x_colnames[j], sep = "-")
            }
        }
    }
    colnames(outmat) = to_colnames
    outmat
}
