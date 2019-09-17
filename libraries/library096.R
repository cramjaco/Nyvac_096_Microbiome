## Color-blind friendly color balettes
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",  "#0072B2", "#D55E00", "#CC79A7", "#F0E442")

### Custom phyloseq operators

## The following three functions are for joining pieces of phyloseq objects

## Joins two sets of sample data and retuns them as type sample data
sample_join <- function(x, y, ...) {
    require(tibble)
    require(phyloseq)
    loc <- left_join(rownames_to_column(data.frame(x)), y, ...)
    loc2 <- column_to_rownames(loc, var = "rowname")
    loc3 <- sample_data(loc2)
    return(loc3)
}

## Joins two sets of taxa data, and  returns them as type taxa data.
## One must include "by" in the ...
taxa_join <- function(x, y, ...) {
    require(tibble)
    require(phyloseq)
    loc <- left_join(rownames_to_column(data.frame(x)), y, ...)
    loc2 <- column_to_rownames(loc, var = "rowname")
    loc3 <- tax_table(loc2)
    return(loc3)
}

## Joins either sample data or taxa data into a phyloseq object
phylo_join <- function(ps, y, type = 'sample', ...) {
    require(phyloseq)
    switch(type, 
           sample = {
        pssd <- sample_data(ps)
        lj <- sample_join(pssd, y, ...)
        mp <- merge_phyloseq(ps, sample_data(lj)) 
    # I don't know why lj can't just go in by itself, it should already be sample data
    }, 
           taxa = {
        pssd <- tax_table(ps)
        lj <- taxa_join(pssd, y, ...)
        mp <- merge_phyloseq(ps, tax_table(lj)) 
    })
    
    return(mp)
}

## Applies a transformation, usually clr, to a whole otu table.
transform_otu_table <- function(ps, f, ...){
    otu1 <- f(otu_table(ps))
    otu1
    ps1 <- ps
    ps1@otu_table <- otu1
    ps1
}

## Sometimes my phyloseq objects have a "tag" taxonomic field that I am re-defining.
## This removes that "tag" field.
remove_tag_phyloseq <- function(ps){
    # tx <- ps %>% tax_table %>% as.data.frame %>% rownames_to_column %>%
    # dplyr::select(-tag) %>% column_to_rownames(var = 'rowname') %>% as.matrix %>% tax_table
    tx0 <- as.data.frame(ps@tax_table@.Data)
    tx1 <- rownames_to_column(tx0)
    tx <- tx1 %>%
    dplyr::select(-tag) %>% column_to_rownames(var = 'rowname') %>% as.matrix %>% tax_table
    
    ps@tax_table <- tx
    
    ps
}

## The tip_glom_saveid function agglomerates phyloseq objects as in phyloseq::tip_glom, with two key additions
## (1) It saves the taxonomic information of the agglomerated taxa into an oldGroups field in the taxonomic table.
## (2) It allows the user to specify "k" (the number of agglomerated groups the user wants) rather than "h" (phylogenetic tree height for agglomeration).

tip_glom_saveid <- function(physeq, hcfun = agnes, h = NULL, k = NULL, ...){

    ## Necessary libraries
    require(phyloseq)
    require(ape)
    require(cluster)
    require(dplyr)
    require(tibble)

    ## Original phyloseq object
    physeq0 <- physeq

    ## From phyloseq::tip_glom
    dd = as.dist(cophenetic.phylo(phy_tree(physeq)))
  #
    if(!is.null(h) & !is.null(k)){
        error("you must specify h XOR k, not both")
    }else if(!is.null(h)){
        psclust = cutree(as.hclust(hcfun(dd, ...)), h = h) 
    }else if(!is.null(k)){
        psclust = cutree(as.hclust(hcfun(dd, ...)), k = k) 
    } else {
        error("you must specify h or k, not neither")
    }
      
    ## Key showing which OTUs are in which groups
    taxGlomKey <- data.frame(key = psclust) %>% rownames_to_column

    ## Key comparing new otus to old otus
    nGrp = max(taxGlomKey$key) 

    grps = vector("list", nGrp)
    for (loc_key in 1:nGrp){
        grps[[loc_key]] = taxGlomKey %>% filter(key == loc_key) %>% pull(rowname)
    }
    

    ## Turn list of lists into list of concatenated strings.
    oldGroups <- sapply(grps,
                        function(x){
                            unlist(
                                paste(x, collapse = ', ')
                            )
                        }
                        )
    
    ## From phyloseq::tip_glom - do the actual merging
    cliques = levels(factor(psclust))[tapply(psclust, factor(psclust), 
                                             function(x) {
                                                 length(x) > 1
                                             })]
    for (i in cliques) {
        physeq = merge_taxa(physeq, eqtaxa = names(psclust)[psclust == 
                                                            i])
    }

    ## the following crashes in r 6.1, rewrite as per https://github.com/joey711/phyloseq/issues/983 21 Aug 2019
    #glomReps <- physeq %>% tax_table %>% as.data.frame %>% rownames
    glomReps <- rownames(as.data.frame(physeq@tax_table@.Data))
    
    
    ## order that the old groups should go in in the new object
    sapply(glomReps, function(y){
        (which(sapply(grps, function(x){y %in% x})))             
    }) -> grpOrder
    
    ## rewrite 21 Aug 2019
    #physeq %>% tax_table %>% data.frame(oldGroups = oldGroups[grpOrder]) %>% as.matrix %>% tax_table -> gtt
    gtt0 <- as.data.frame(physeq@tax_table@.Data)
    gtt0$oldGroups = oldGroups[grpOrder]
    gtt <- tax_table(as.matrix(gtt0))
    
    physeq@tax_table <- gtt
    
    physeq
    
}

# Fix indenting
     
## Discard taxa not appearing in at least some fraction of the samples.
prevalence_filter_taxa <- function(ps, thresh = 0.1){
    Lprevdf = apply(X = otu_table(ps),
                 MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
    # Add taxonomy and total read counts to this data.frame
    Lprevdf = data.frame(Prevalence = Lprevdf,
                      TotalAbundance = taxa_sums(ps),
                      tax_table(ps))
    LprevalenceThreshold = thresh * nsamples(ps)
    keepTaxa = rownames(Lprevdf)[(Lprevdf$Prevalence >= LprevalenceThreshold)]
    keepTaxa
    prune_taxa(keepTaxa, ps)
}

## A function that just passes its input through. Useful with pipes.
pass <- function(x){x}
     
## recode a vector of continuous values into a 0 1 vector where 1s are above or equal to median 
## and zeros are below
medcode <- function(vec){
    sapply(vec, function(x){
      if(is.finite(x)){
    if(x < median(na.omit(vec))){0}else{1}}else{NaN}
    }
           )
           }
     
## recode a vector of continuous values into a 0 1 vector where "High" are above or
##equal to median and "Low" are below
medcode_hl <- function(vec){
    sapply(vec, function(x){
    if(x < median(na.omit(vec))){"Low"}else{"High"}
    }
           )
           }

# do nothing if a vector of zeros and ones
medcode2 <- function(vec)
    if(all(vec %in% c(0,1))){vec}else{medcode(vec)}
 
## Box-cox transform a variable of interest.
jac_box_cox <- function(vec){
    require(car)
    pt <- car::powerTransform(vec)
    car::bcPower(vec, pt$roundlam)
    }
    

## take a phyloseq object, and replace the old rownames of the tax_table
##and columnames of the otu table and the tip labels
## of the phylo tree with a new set of names
swap.phyloseq.taxnames <- function(ps, oldname = 'Sequence', newname = 'tag'){
                                        # make sure this is a phyloseq object with otu and taxonomy table
    if(class(ps) != 'phyloseq'){ stop('ps must be a phyloseq object')}
    if(is.null(tax_table(ps))){stop ('ps must have a tax_table')}
    if(is.null(tax_table(ps))){stop ('ps must have a tax_table')}
    
    ps2 <- ps
    
                                        # Get old and new names
    #tt <- tax_table(ps)
    tt <- as.data.frame(ps@tax_table@.Data)
    oldNames <- rownames(tt)
    newNames <- as.vector(tt[,newname])
    
                                        # Rename taxonomy table
    tt <- data.frame(tt, OldRN = oldNames)
    colnames(tt)[which(colnames(tt) == 'OldRN')] <- oldname
    tt <- tax_table(as.matrix(tt))
    ps2@tax_table <- tt
    
    taxa_names(ps2) <- newNames
    
    ps2
}

                                        # Bootstrap proportionality value phi
boot_phi <- function(data, indices){
                                        # Rel is a sample X Taxon matrix of relative abundances. 
                                        # No zeros. Rather, one shold pad with detection threshold, eg 0.001.
                                        # returns sol. This has phi values for every pair of taxa. Phi is the vector we are trying to bootstrap.
    
    require(dplyr)    
    rel.boot <- data[indices,]
    sol <- rel.boot %>% make_proportionality_matrix %>% 
        as.data.frame %>%
        rownames_to_column("TaxonX") %>% gather(TaxonY, phi, -TaxonX) %>%
        filter(TaxonX != TaxonY)
    sol$phi
}

## Convert relative abundance data into a proportionality matrix
make_proportionality_matrix <- function(LRel){
    LRel.clr <- unwarn_null(compositions::clr(LRel))

    LRel.vlr <- unwarn_null(variation(acomp(LRel)))

    LRel.clr.var <- apply(LRel.clr, 2, compositions::var)

    LRel.phi <- sweep(LRel.vlr, 2, LRel.clr.var, FUN="/")
    LRel.phi
}

## Suppress warnings
## for reasons I dont' understand, this actually represses all error messages,
## not just specified warnstr
unwarn <- function(f, warnstr = ""){
    withCallingHandlers(f,
                        warning = function(w){
                            if(any(grepl(warnstr, w))){
                                invokeRestart("muffleWarning")
                            }
                        }
                        )
    
}

                                        # library compositions is "NULL cannot have attributes warning central
                                        # http://romainfrancois.blog.free.fr/index.php?post/2009/05/20/Disable-specific-warnings

unwarn_null_core <- function(w){
    if(any(grepl("as NULL cannot have attributes", w))){
        invokeRestart("muffleWarning")
    }
}

unwarn_null <- function(f){
    withCallingHandlers(f,
                        warning = unwarn_null_core)
    
}

                                        #Convert "male" and "female" type strings into ones for males and 0 for females.
testIsMale <- function(sex){
    if(sex %in% c('Male', 'male', 'm', 'M')){
        1
    }else if(sex %in% c('Female', 'female', 'f', 'F')){
        0
    }else{
        NaN
    }
}

                                        # Apply testIsMale over a vector, rather than just to one thing.
testIsMaleVec <- function(sexVec){sapply(sexVec, testIsMale)}

### Taxon tagging function
## It has more than a few pieces

tag_taxon <- function(tvdf){
                                        # take a vector of taxonomic identities phylum at left, species at right
                                        # return the right-most name
    species <- tvdf %>% dplyr::select(Species) %>% unlist %>% as.character
    
    genus2 <- tvdf %>% dplyr::select(Genus, Genus.y) %>% unlist
    genus <- genus2 %>% na.omit %>% .[1] %>% as.character
    
    lowTax <- tvdf %>% dplyr::select(Kingdom:Family) %>% unlist
    lowTax1 <- lowTax %>% na.omit
    lowTaxTag <- lowTax1 %>% .[length(.)] %>% as.character

    if(!is.na(species)){
        out <- paste(genus, species)
    }else if(!is.na(genus)){
        out <- genus
    }else{
        out <- lowTaxTag
    }
    out
}

fix_EColi <- function(x){
    if(x == 'Escherichia/Shigella albertii/boydii/coli/coli,/dysenteriae/enterica/fergusonii/flexneri/sonnei/vulneris'){
        out <- 'Escrichia coli et al.'
    }else{
        out <- x
    }
    out
}

fix_LongName <- function(x){
    if(nchar(x) > 35){
        out <- paste(strtrim(x, width = 30), 'et al.')
    }else{
        out <- x
    }
    out
}

tag_taxon2 <- function(tvdf){
    pre <- tax_taxon(tvdf)
    out <- fix_EColi(pre)
    out
}

fix_EColi_in_vec <- function(vec){
    sapply(vec, fix_EColi)
}

fix_LongName_in_vec <- function(vec){
    sapply(vec, fix_LongName)
}

### Coming up with good names for agglomerated phyloseq groups (or even species level OTUs)

tag_phyloseq <- function(ps){
    ## name and deduplicate the names for each taxonomic group in a phyloseq object
    ## requires a phyloseq object with a tax table with Genus, Genus.y, and Species in it
    
    ## this always yells "Setting row names on a tibble is deprecated." I'd disable that warning if I could.
    require(phyloseq)
    require(dplyr)
    require(purrrlyr)
    
    sequences <- rownames(tax_table(ps))
    
    #ps %>% tax_table %>% as.data.frame %>%
tx0 <- as.data.frame(ps@tax_table@.Data)
tx0 %>%
        by_row(tag_taxon, .to = "tag") %>%
        mutate(tag = fix_LongName_in_vec(fix_EColi_in_vec(tag))) %>%
        mutate(tag = make.unique(as.character(tag)), sequence = sequences) %>%
        column_to_rownames(var = 'sequence') %>% #
        as.matrix %>%
        tax_table -> TTMod
    
    psOut <- ps
    psOut@tax_table <- TTMod
    psOut
    
}
## apply dplyr::mutate to the sample data of a phyloseq object.
mutate_phyloseq_sample <- function(ps, ...){
    require(phyloseq)
    require(dplyr)
    ps %>% sample_data %>% as('data.frame') %>% rownames_to_column %>%
        mutate(...) %>%
        column_to_rownames('rowname') %>%
        sample_data %>%
        pass-> sd
    
    ps1 <- ps
    ps1@sam_data <- sd
    
    ps1
    
}

# https://support.bioconductor.org/p/74637/
# I am doing the robust method because sometimes p-values don't extend to one (in the local tests). Rather thay only extend to .8 or so. This may be an artifact of OTU compositionality.
robust_qvalue <- function(p){
    tryCatch({qvalue(p)
             },error = function(e){
        qvalue(p, lambda = seq(0.05, max(p), 0.05))
    })
}

p2q <- function(p){
    require(qvalue)
    tryCatch({
    qres <- robust_qvalue(p)
    qres$q
        },error = function(e){
        rep(NaN, length(p))
    })
}
 
format_round <- function(num, digits = 3){
    format(round(num, digits), nsmall = digits)
}