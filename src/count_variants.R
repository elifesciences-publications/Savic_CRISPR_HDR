library("readxl")
library("writexl")
library("CrispRVariants")
library("GenomicAlignments") # Needed for variant classification
library("rtracklayer") # For reading in guide sequences
library("Rsamtools") # For reading in reference sequences
library("dplyr") # Batch convert columns to factor
library("tibble") # For convenience function rownames_to_column

# --------------------------------------------------------------------------
# Read the metadata table, guide locations and reference sequences

#metadata <- read_excel("metadata.xlsx")
metadata <- read_excel("metadata/metadata.xlsx")

metadata <- metadata %>% mutate_at(c("Experiment", "Amplicon",
                                     "Replicate", "Condition"),
                                   funs(factor(.)))

bam_fnames <- sprintf("bam/%s/%s.bam", metadata$Experiment, metadata$Sample_ID)
cat(sprintf("All files exist: %s\n", all(file.exists(bam_fnames))))
metadata$bam_fname <- bam_fnames

amplicons <- FaFile("amplicons/amplicons.fa")
guide_locs <- import("guides.bed")
# guide_locs <- guide_locs + 45
names(guide_locs) <- guide_locs$name
guide_seqs <- getSeq(amplicons, guide_locs)

metadata <- split(metadata, metadata$Experiment, drop = TRUE)

# --------------------------------------------------------------------------
# Create CrisprSet objects

corrected <- c("hCXCR4" = "SNV:-1T,1C", "hEMX1" = "SNV:-1C,1T,2C",
               "hHBB" = "SNV:-1C,1C", "mPcsk9" = "SNV:-1G,1C",
               "mR26_mESC" = "SNV:-1G,1C", "mutRFP_Reporter" = "SNV:2T,3A")
               
csets <- lapply(seq_along(metadata), function(i){
  print(i)
  nm <- names(metadata)[i]
  md <- metadata[[i]]

  ref_nm <- unique(md$Amplicon)
  ref <- guide_seqs[ref_nm][[1]] # Subsetting to make a DNAString not DNAStringSet
  gd <- guide_locs[ref_nm]
  cset <- readsToTarget(md$bam_fname, target = gd, reference = ref, verbose = FALSE,
                        names = md$Sample_name, target.loc = width(gd) - 6)
  
  amplicon <- unique(md$Amplicon)
  if (amplicon %in% names(corrected)){
    cset_alns <- alns(cset)
    allele_labels <- mcols(unlist(cset_alns))$allele
    corrected_lab <- corrected[[amplicon]]
    allele_labels[allele_labels == corrected_lab] <- "Corrected"
    cset$setCigarLabels(labels = relist(allele_labels, cset_alns))
  }
  
  cset
})

names(csets) <- names(metadata)

# --------------------------------------------------------------------------
# Make figures

# lower cutoff for indels 

snv_indel <- list(c(0.05,0.01), c(0.5,0.01), c(0.5,0.5))
top.n <- 50



order_alleles <- function(alleles, amplicon){
    allele_ord <- c("no variant")
    if ("Corrected" %in% names(corrected)){
      allele_ord <- c(allele_ord, "Corrected")
    }
    alleles <- alleles[! alleles %in% allele_ord]
    is_snv <- grepl("SNV", alleles)
    allele_ord <- c(allele_ord, alleles[!is_snv], alleles[is_snv])
    
    # Other last
    if ("Other" %in% allele_ord){
      allele_ord <- c(allele_ord[! allele_ord == "Other"], "Other")
    }
    allele_ord
}

  
# Note that variantCounts function works by rowSums, not rowMeans
dummy <- lapply(names(csets), function(nm){
  print(nm)
  md <- metadata[[nm]]
  cset <- csets[[nm]]
  vc <- variantCounts(cset, result = "proportions")
  rmaxs <- apply(vc, 1, max)
  top_alleles <- rownames(vc)[order(rmaxs, decreasing = TRUE)[1:top.n]]
  indel <- grepl("I|D", rownames(vc))
  
  prev <- NULL
  for (cutoffs in snv_indel){
    snv_cutoff <- cutoffs[1]
    indel_cutoff <- cutoffs[2]
    
    outfname <- sprintf("figures/%s_snv_%s_indel%s_sortby_batch.pdf",
                        nm, gsub("\\.", "_", snv_cutoff),
                        gsub("\\.", "_", indel_cutoff))

    keep <- rmaxs >= snv_cutoff | (rmaxs >= indel_cutoff & indel) & 
            rownames(vc) %in% top_alleles
    alleles <- rownames(vc)[keep]
    alleles <- order_alleles(alleles, unique(md$Amplicon))   
    
    # some figure parameters chosen by trial and error
    xangle <- ifelse(length(alleles) > 10, 45, 25)
    xangle <- ifelse(nm %in% c("HEK293T_EMX1"), 90, xangle)
    
    ptsz <- ifelse(nm %in% c("Reporter_mutRFP") & indel_cutoff == 0.01, 1.5, 2)
    atsz <- ifelse(nm %in% c("Reporter_mutRFP") & indel_cutoff == 0.01, 6, 8)
    
    pa_args <- list(style = "mismatches", legend.cols = 4,
                    ins.size = 3, legend.text.size = 6,
                    alleles = alleles, top.n = nrow(vc),
                    min.insertion.freq = 15, axis.text.size = atsz)
    hm_args <- list(type = "proportions", x.angle = xangle, top.n = nrow(vc),
                    group = md$Replicate,
                    plot.text.size = ptsz, header = "counts", alleles = alleles)
  
    fig_height <- max(min(0.35 * length(alleles), 11), 3.5)
    
    # Figure in main paper
    if (nm %in% c("Reporter_mutRFP","Reporter_mutRFP-On_target") & indel_cutoff == 0.5){ 
      fig_height <- 5
      pa_args$plot.text.size <- 3
      hm_args$plot.text.size <- 3
      
      pdf(outfname, width = 10, height = fig_height, useDingbats = FALSE)

      plotVariants(cset, plotAlignments.args = pa_args,
                   plotFreqHeatmap.args = hm_args,
                   col.wdth.ratio = c(1.1,1),
                   left.plot.margin =  grid::unit(c(0.2,0,8,0.5), "lines"))
      dev.off()

    } else {
      p <- do.call(plotAlignments, c(list(obj = cset), pa_args))
      p <- p + ggtitle(gsub("_", " ", nm))
      q <- do.call(plotFreqHeatmap, c(list(obj = cset), hm_args))
      gene_p <- grid::grid.rect(gp=grid::gpar(col="white"), draw = FALSE)

      pdf(outfname, width = 10, height = fig_height, useDingbats = FALSE)

      CrispRVariants:::arrangePlots(gene_p, p,q,
                                  col.wdth.ratio = c(1,1),
                                  row.ht.ratio = c(0.5, fig_height - 0.5),
                                  left.plot.margin =  grid::unit(c(0.2,0,8,0.5), "lines"))
      dev.off()
    }
  }
})


# --------------------------------------------------------------------------
# Tables of all variants

# Write tables of all variants into multisheet excel file
vcs <- lapply(csets, function(cset){
          data.frame(variantCounts(cset)) %>% rownames_to_column("Sample")
})


write_xlsx(vcs, path = "counts_tables/all_variants.xlsx")


# --------------------------------------------------------------------------
# Tables of variants by category

# Classify variants according to whether the correction was successful
classify_params <- list("hHBB" = list(allele.loc = IRanges(17, 18),
                                      ref.allele = "GT", target.allele = "CC"),

                        "hEMX1" = list(allele.loc = IRanges(17, 19),
                                       ref.allele = "AGG", target.allele = "CTC"),

                        "hCXCR4" = list(allele.loc = IRanges(17, 18), 
                                        ref.allele = "GA", target.allele = "TC"),

                        "mESC_Pcsk9" = list(allele.loc = IRanges(16, 17),
                                            ref.allele = "CA", target.allele = "GC"),

                        "mutRFP_Reporter" = list(allele.loc = IRanges(18,19),
                                                 ref.allele = "TA", target.allele = "CT"),

                        "mR26_mESC" = list(allele.loc = IRanges(16,17),
                                           ref.allele = "GA", target.allele = "CT"),
                        
                        "mPcsk9" = list(allele.loc = IRanges(17,18),
                                        ref.allele = "CA", target.allele = "GC")
                        )


classify_variants <- function(cset, allele.loc, ref.allele, target.allele){
    alns <- alns(cset)
    all_alns <- unlist(unname(alns))
    sqs <- CrispRVariants:::seqsToAln(cigar(all_alns), mcols(all_alns)$seq,
                                      cset$target, aln_start = start(all_alns))

    configs <- substr(sqs, start(allele.loc), end(allele.loc))
    is_ref <- configs == ref.allele
    is_target <- configs == target.allele

    dels <- cigarRangesAlongReferenceSpace(cigar(all_alns), ops = "D")
    ins <- cigarRangesAlongQuerySpace(cigar(all_alns), ops = "I")
    delta_size <- sum(width(ins)) - sum(width(dels))
    is_indel <- lengths(dels) > 0 | lengths(ins) > 0
    is_frameshift <- is_indel & abs(delta_size) %% 3 != 0
    
    allele_labs <- rep("", length(all_alns))
    allele_labs[is_target & is_frameshift] <- "Corrected plus frameshift indel"
    allele_labs[is_target &! is_frameshift] <- "Corrected plus inframe indel"
    allele_labs[is_target &! is_indel] <- "Perfect correction"
    allele_labs[!is_target & is_frameshift] <- "Frameshift indel"
    allele_labs[!is_target & ! is_frameshift] <- "Inframe indel"
    allele_labs[!is_target & ! is_indel & ! is_ref] <- "Mismatch"
    allele_labs[is_ref] <- "Reference"

    allele_levs <- c("Reference", "Perfect correction",
                     "Corrected plus inframe indel",
                     "Corrected plus frameshift indel",
                     "Inframe indel", "Frameshift indel", "Mismatch")
    allele_labs <- relist(allele_labs, alns)
    
    temp <- lapply(allele_labs, function(x){
                    x <- factor(x, levels = allele_levs)
                    table(x)
                  })
    result <- data.frame(do.call(rbind, temp))
    result$Total <- rowSums(result)
    result
}



by_cat <- lapply(names(csets), function(nm){
  cset <- csets[[nm]]
  gd <- cset$target$name
  if (! gd %in% names(classify_params)) return()

  pms <- classify_params[[gd]]
  pms$cset <- cset
  res <- do.call(classify_variants, pms)
  res
})

names(by_cat) <- names(csets)
by_cat <- by_cat[! sapply(by_cat, is.null)]
by_cat <- lapply(by_cat, function(bc) bc %>% rownames_to_column("Sample"))

write_xlsx(by_cat,  path = "counts_tables/by_category.xlsx")

