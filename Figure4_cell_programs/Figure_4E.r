# Code for Figure 4E

library(enrichR)
library(data.table)
library(cowplot)
library(ggplot2)
source('colors.r')

# ----------------
# Useful functions
# ----------------

fast_enrich = function(genes, regex='GO_.*2017$|KEGG.*2016|Reactome.*2016|Panther_2016', collapse=FALSE){
    
    # This function uses EnrichR to calculate gene set enrichments for a set of genes
    # genes = input list of genes
    # regex = regex matching the gene sets to use
    # collapse = combine results using rbind
    
    # For full list of databases: listEnrichrDbs()
    # Computes enrichment with Fisher exact test
    # Also uses random gene sets to correct the Fisher test
    
    # Select databases to use
    dbs = grep(regex, listEnrichrDbs()[,3], value=T)
    
    # Run enrichment test
    res = enrichr(genes, dbs)
    
    # Fix each list
    res = sapply(res, function(a) {
        
        # Sort by adjusted p-value
        a = as.data.table(a)[,.(Term, Overlap, P.value, Adjusted.P.value, Genes)][order(Adjusted.P.value)]
	
        # Get overlap statistics
        a[, NM := as.integer(gsub('/.*', '', Overlap))]
        a[, NQ := length(genes)]
        a[, NT := as.integer(gsub('.*/', '', Overlap))]
	
    }, simplify=F)
    
    # Return results
    if(collapse == TRUE){res = do.call(rbind, res)[order(P.value)]}
    res
}


load_kegg = function(names.use=NULL, names.rmv=NULL, do.flatten=FALSE, no_spaces=FALSE, regex=NULL){

    # This function conveniently loads different KEGG terms
    # names.use = KEGG terms to use
    # names.rmv = KEGG terms to remove
    # do.flatten = return flat list
    # no_species = remove spaces from KEGG term names
    # regex = regex specifying KEGG terms to keep
    
    kegg = readRDS('kegg.rds')
    
    kegg = sapply(names(kegg), function(A)
               sapply(names(kegg[[A]]), function(B)
                   sapply(names(kegg[[A]][[B]]), function(C) {
                       if(!is.null(regex)){if((! grepl(regex, A) & ! grepl(regex, B)) & ! grepl(regex, C)){return(NULL)}}
                       if(!is.null(names.use)){if(!(A %in% names.use | B %in% names.use | C %in% names.use)){return(NULL)}}
                       if(!is.null(names.rmv)){if(A %in% names.rmv | B %in% names.rmv | C %in% names.rmv){return(NULL)}}
                       as.character(na.omit(kegg[[A]][[B]][[C]]))
                   }, simplify=F),
               simplify=F),
           simplify=F)
    
    if(no_spaces == TRUE){
        for(A in names(kegg)){
            for(B in names(kegg[[A]])){
                names(kegg[[A]][[B]]) = gsub(' ', '_', names(kegg[[A]][[B]]))
            }
            names(kegg[[A]]) = gsub(' ', '_', names(kegg[[A]]))
        }
        names(kegg) = gsub(' ', '_', names(kegg))
    }
    
    if(do.flatten == TRUE){
        kegg = unlist(unlist(kegg, recursive=FALSE), recursive=FALSE)
        kegg = kegg[!sapply(kegg, is.null)]
    }
    kegg
}


# ------------
# 1. Load data
# ------------
# Load the DE statistics for double-positive (DP) vs. double-negative (DN) cells
# These were calculated using MAST (see Methods), but are provided here for convenience

# First, define the datasets that we will use. Full references available in Supplementary Tables
# ----------------------------------------------------------------------------------------------
# gut = Smillie et al, Cell, 2020
# lung = Regev/Rajagopal
# lca = Teichmann
# nose = Ordovas-Montanes et al, Nature, 2018
datasets.use = c('gut', 'lung', 'lca', 'nose')

# Load gene modules for DP cells for gut, lung, and nose
# ------------------------------------------------------
mod = sapply(datasets.use, function(a) readRDS(paste0(a, '.ACE2_TMPRSS2.rds')), simplify=F)
mod = unlist(mod, recursive=F)
mod = mod[lengths(mod) > 0]
mod = sapply(mod, function(a) a$mast, simplify=F)


# -----------------------------------------------
# 2. Use enrichR to calculate enriched KEGG terms
# -----------------------------------------------

# This function filters the DE results as follows:
# - require that the adjusted p-value (padjH) be less than 0.05
# - require that the fraction of expressing cells (alpha) be greater than 0.10
# - require that the estimated fold-change (mastfc) be greater than 0
# - sort results by decreasing estimated fold-change
filter_de = function(x){x[padjH < .05][alpha > .1][mastfc > 0][order(-mastfc)]}

# This function uses enrichR to calculate the enriched KEGG terms for each DE list
# Input arguments:
# - regex = Regex that matches the datasets we want to use
# - k = number of genes to use from each DE list for KEGG enrichment
get_kegg = function(regex='.*', k=25){
    # Select top k genes from each module matching the provided regex
    genes.use = grep('^RP', sort(unique(sapply(mod[grep(regex, names(mod))], function(a) filter_de(a)[1:k]$gene))), value=T, invert=T)
    # Use enrichR to calculate the top enriched KEGG terms
    fast_enrich(genes.use)
}

# Run KEGG enrichments for all tissues, gut, lung, and nose
k.all = get_kegg('.*', k=25)$KEGG
k.gut = get_kegg('gut', k=50)$KEGG
k.lung = get_kegg('lung|lca', k=50)$KEGG
k.nose = get_kegg('nose', k=50)$KEGG


# --------------------------------------------------------
# 3. Make barplot of top 5 enriched KEGG terms (Figure 4E) 
# --------------------------------------------------------

# Remove "Human Diseases" KEGG terms
disease_terms = unname(unlist(sapply(load_kegg()$'Human Diseases', names)))

make_barplot = function(kegg, disease=FALSE, nterms=10, fill, lab){
    # This function uses ggplot to make a barplot of the most significant KEGG terms

    # Input arguments:
    # kegg = KEGG enrichment table
    # disease = include disease terms?
    # nterms = number of terms to include
    # fill = fill color
    # lab = labels
    
    # filter kegg results
    x = kegg[Adjusted.P.value < .05][, Term := gsub(' Homo.*', '', Term)][,.(Term, -log10(Adjusted.P.value))]
    
    # remove disease terms
    i = gsub(' Homo.*', '', x$Term) %in% disease_terms    
    if(disease == FALSE){
        x = x[! i]
    } else {
        x = x[i]
    }
    
    # select number of terms
    x = x[1:min(nterms, nrow(x))]
    x$Term = factor(x$Term, levels=rev(x$Term))
    
    # make ggplot barplot
    ggplot(x) + geom_bar(aes(x=Term, y=V2), stat='identity', fill=fill) + coord_flip() + xlab(lab) + ylab('-log10(Adjusted p-value)')
    
}

# Make barplots of top 5 enriched KEGG terms for all tissues (disease + non-disease), nose, and lung
p1 = make_barplot(k.all, disease=TRUE, nterms=5, fill=set.colors[[1]], lab='All (disease)') + ylab('')
p2 = make_barplot(k.all, disease=FALSE, nterms=5, fill=set.colors[[2]], lab='All (non-disease)') + ylab('')
p3 = make_barplot(k.nose, disease=FALSE, nterms=5, fill=set.colors[[3]], lab='Nose') + ylab('')
p4 = make_barplot(k.lung, disease=FALSE, nterms=5, fill=set.colors[[4]], lab='Lung')
ps = plot_grid(p1,p2,p3,p4,nrow=4,align='v')
save_plot(ps, file='Figure_4E.pdf', nrow=2, ncol=2)

