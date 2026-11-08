#' Read and filter a VEP-annotated VCF file
#'
#' @description
#' Reads a VEP-annotated VCF file for a given sample and applies a series of
#' filters: genotype quality (GQ), read depth (DP), allelic bias (VAF >= 0.25),
#' variant consequence (HIGH, MODERATE, splicing), synonymous variants, and a
#' blacklist of recurrent artefacts. Returns a \code{data.frame} ready to be
#' passed to \code{priorBestVariant()}.
#' 
#' #' One row is returned per variant and gene. Multi-allelic records are split
#' with \code{bcftools norm -m -any}, which does not split the CSQ field, so
#' only the CSQ entries of the row's own allele are used. These entries are
#' grouped by gene (\code{SYMBOL}, or \code{Gene} when empty) and, in each
#' gene, the entry kept is chosen by: protein-coding biotype, then the most
#' severe consequence (VEP consequence order), then MANE Select, canonical
#' transcript, Ensembl gene ID (\code{ENSG}) and the original order. The CSQ
#' field of each row keeps only that entry, so \code{priorBestVariant()} uses
#' the same consequence. Entries without gene (e.g. intergenic) give one row.
#'
#' @param sampleName Character. Sample identifier as written in the VCF file.
#' @param filter Character. FILTER field value to select variants (e.g.
#'   \code{"PASS"}). Use \code{""} to skip this filter. Default \code{""}.
#' @param geneQuality Numeric. Minimum genotype quality (GQ) threshold.
#'   Default \code{20}.
#' @param readDepth Numeric. Minimum read depth (DP) threshold. Default \code{10}.
#' @param vcfFile Character. Path to the bgzipped, VEP-annotated VCF file.
#' @param assembly Character. Genome assembly version. Either \code{"assembly38"}
#'   (default, hg38) or \code{"assembly37"} (hg19).
#' @param distSplicThreshold Integer or character. Maximum distance (bp) from the
#'   exon-intron boundary for intronic variants to be retained. Use \code{"full"}
#'   or \code{Inf} to include all intronic variants regardless of distance.
#'   Default \code{"full"}.
#' @param synonymous Logical. If \code{TRUE} (default), synonymous variants are
#'   included. If \code{FALSE}, synonymous variants without splicing consequence
#'   are discarded.
#'
#' @return A \code{data.frame} with the filtered variants, with columns
#'   \code{CHROM}, \code{POS}, \code{REF}, \code{ALT}, \code{FILTER},
#'   \code{INFO}, \code{GT}, and optionally \code{GQ}, \code{DP}, \code{AD},
#'   \code{GT_full} depending on the FORMAT fields present in the VCF.
#'   The attribute \code{vcfFile} stores the path to the original VCF file,
#'   required by \code{priorBestVariant()}.
#'
#' @export
#' @importFrom stringi stri_replace_first_regex stri_split_fixed stri_list2matrix stri_locate_first_regex stri_sub
#' @importFrom data.table fread
#' @importFrom utils read.csv
#'
#' @examples
#' \dontrun{
#' variants <- readVCF(sampleName = "sample1",
#'                     vcfFile    = "sample.vep.vcf.gz",
#'                     assembly   = "assembly38",
#'                     geneQuality = 20,
#'                     readDepth   = 10,
#'                     distSplicThreshold = "full",
#'                     synonymous  = TRUE)
#' }
readVCF <- function(sampleName = "", filter="", geneQuality=20, readDepth=10,
                    vcfFile = "", assembly="assembly38", distSplicThreshold="full", synonymous=TRUE) {
 
  # 1. Load Blacklist
  dest <- system.file(paste("extdata", assembly, sep="/"), package="ClinPrior")
  blacklistFile <- list.files(path=dest, full.names=TRUE)[grep("blacklist", list.files(dest))]
  blacklist <- read.csv(blacklistFile, sep="\t", header=FALSE)
  ensemble  <- do.call(paste, c(blacklist, sep="-"))
  
  blacklist_parts <- strsplit(ensemble, "-")
  bl_chrompos <- sapply(blacklist_parts, function(x) paste(x[1], x[2], sep="-"))
  bl_ref      <- sapply(blacklist_parts, function(x) x[3])
  bl_alt      <- sapply(blacklist_parts, function(x) x[4])
  
  # 1. Check which FORMAT fields exist in the header
  header_fmt <- system(paste("bcftools view -h", vcfFile, "| grep '^##FORMAT'"), intern=TRUE)
  has_GQ <- any(grepl("ID=GQ", header_fmt))
  has_DP <- any(grepl("ID=DP[,>\"]", header_fmt))
  has_AD <- any(grepl("ID=AD", header_fmt))
  
  # 2. Build bcftools filter -i logic
  format_filters <- c()
  if(has_GQ) format_filters <- c(format_filters, paste0("FORMAT/GQ[0] >= ", geneQuality))
  if(has_DP) format_filters <- c(format_filters, paste0("FORMAT/DP[0] >= ", readDepth))
  if(filter != "") format_filters <- c(format_filters, paste0("FILTER=\"", filter, "\""))
  
  full_filter <- if(length(format_filters) > 0) paste(format_filters, collapse=" && ") else ""
  
  # 3. Build FORMAT fields for query and column names
  fmt_fields  <- "[%GT"
  col_names   <- c("CHROM","POS","REF","ALT","FILTER","INFO","GT")
  if(has_GQ) { fmt_fields <- paste0(fmt_fields, "\t%GQ"); col_names <- c(col_names, "GQ") }
  if(has_DP) { fmt_fields <- paste0(fmt_fields, "\t%DP"); col_names <- c(col_names, "DP") }
  if(has_AD) { fmt_fields <- paste0(fmt_fields, "\t%AD"); col_names <- c(col_names, "AD") }
  # Afegeix al final de fmt_fields, DINS dels brackets, el camp complet
  # bcftools suporta multiples expressions dins de [ ]
  # Check which optional FORMAT fields exist before including in GT_full
  has_PGT <- any(grepl("ID=PGT", header_fmt))
  has_PID <- any(grepl("ID=PID", header_fmt))
  has_PL  <- any(grepl("ID=PL",  header_fmt))
  has_PS  <- any(grepl("ID=PS",  header_fmt))
  
  gt_full_parts <- c("%GT")
  if(has_GQ)  gt_full_parts <- c(gt_full_parts, "%GQ")
  if(has_DP)  gt_full_parts <- c(gt_full_parts, "%DP")
  if(has_AD)  gt_full_parts <- c(gt_full_parts, "%AD")
  if(has_PGT) gt_full_parts <- c(gt_full_parts, "%PGT")
  if(has_PID) gt_full_parts <- c(gt_full_parts, "%PID")
  if(has_PL)  gt_full_parts <- c(gt_full_parts, "%PL")
  if(has_PS)  gt_full_parts <- c(gt_full_parts, "%PS")
  
  fmt_fields <- paste0(fmt_fields, "\t", paste(gt_full_parts, collapse=":"), "]")
  col_names  <- c(col_names, "GT_full")
  
  # 4. Create temporary file and run system command
  tmpTxt <- tempfile(fileext=".txt")
  
  # Build command using bcftools filter -i for all filters
  cmd <- paste0("bcftools view -s ", sampleName, " ", vcfFile,
                if(full_filter != "") paste0(" | bcftools filter -i '", full_filter, "'") else "",
                " | bcftools norm -m -any",
                " | bcftools view -e 'ALT=\"*\"'",
                " | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO\t", fmt_fields, "\n'",
                " > ", tmpTxt)
  
  # Run the pipe with pipefail, so a failure of any bcftools step stops here
  errTxt <- tempfile(fileext=".err")
  status <- system2("bash", c("-o", "pipefail", "-c", shQuote(cmd)), stderr = errTxt)
  err    <- if(file.exists(errTxt)) readLines(errTxt, warn = FALSE) else character(0)
  unlink(errTxt)
  if(!identical(as.integer(status), 0L)) {
    unlink(tmpTxt)
    stop("bcftools failed while reading ", vcfFile, " (exit status ", status, "):\n",
         paste(err[seq_len(min(20, length(err)))], collapse = "\n"), call. = FALSE)
  }
  if(length(err) > 0) message(paste(err, collapse = "\n"))
  
  # 5. Read data into R and delete temporary file
  df <- fread(tmpTxt, 
              sep="\t",
              header=FALSE, 
              col.names=col_names, 
              quote="", 
              fill=TRUE,
              colClasses = "character")  # force all columns as character
  
  df <- as.data.frame(df)
  unlink(tmpTxt)
  
  # Convert numeric columns
  df$POS <- as.integer(df$POS)
  if(has_GQ) df$GQ <- as.numeric(df$GQ)
  if(has_DP) df$DP <- as.numeric(df$DP)
  
  # Basic cleanup
  df <- df[df$GT != "0/0" & df$GT != "0|0", ]
  df$CHROM <- gsub("^chr", "", df$CHROM)
  
  # 6. Parse CSQ (VEP)
  hdr_csq  <- system(paste0("bcftools view -h ", vcfFile, " | grep 'ID=CSQ'"), intern=TRUE)
  meta     <- strsplit(gsub('.*Format: ([^"]+)".*', "\\1", hdr_csq), "[|]")[[1]]
  
  ConsequencePOS <- match("Consequence", meta)
  ImpactPOS      <- match("IMPACT",      meta)
  cDNAPOS        <- match("HGVSc",       meta)
  
  # One row per variant and gene, with the most severe CSQ entry of each gene
  if(nrow(df) > 0) {
    sel     <- .csqPerGene(df$INFO, meta, df$REF, df$ALT)
    df      <- df[sel$row, ]
    df$INFO <- sel$info
  }
  csq_raw   <- stri_replace_first_regex(df$INFO, ".*CSQ=", "")
  csq_first <- stri_replace_first_regex(csq_raw, ",.*", "")
  csq_split <- stri_split_fixed(csq_first, "|")
  
  consequence <- stri_list2matrix(lapply(csq_split, `[`, ConsequencePOS), byrow=TRUE)[,1]
  impact      <- stri_list2matrix(lapply(csq_split, `[`, ImpactPOS), byrow=TRUE)[,1]
  cDNA        <- stri_list2matrix(lapply(csq_split, `[`, cDNAPOS), byrow=TRUE)[,1]
  
  # 7. Splicing and Impact filters
  # Resolve distSplicThreshold: "full" or Inf means no distance limit
  distSplicThreshold <- if(identical(distSplicThreshold, "full") || 
                           (is.numeric(distSplicThreshold) && is.infinite(distSplicThreshold))) {
    Inf
  } else {
    as.numeric(distSplicThreshold)
  }
  
  hgvsc_only <- gsub(".*:(c\\.[^|]+)", "\\1", cDNA)
  m          <- regexpr("[+-][0-9]+", hgvsc_only, perl=TRUE)
  dist_match <- regmatches(hgvsc_only, m, invert=FALSE)
  dist_num   <- rep(NA_real_, length(cDNA))
  dist_num[m > 0] <- as.numeric(dist_match)
  
  pos_intron  <- grepl("intron", consequence)
  keep_intron <- if(is.infinite(distSplicThreshold)) {
    pos_intron  # include all intronic variants regardless of distance
  } else {
    pos_intron & !is.na(dist_num) & abs(dist_num) <= distSplicThreshold
  }
  
  keep <- grepl("HIGH|MODERATE|LOW", impact) | keep_intron
  if(!synonymous) {
    synom_no_splice <- grepl("synonymous", consequence) & !grepl("splice", consequence)
    keep <- keep & !synom_no_splice
  }
  df <- df[keep, ]
  
  # 8. Robust allelic bias filter
  if(has_AD && nrow(df) > 0) {
    ad_split <- strsplit(df$AD, ",")
    
    vaf <- sapply(ad_split, function(a) {
      a <- as.numeric(a)
      if(length(a) < 2) return(NA)
      total <- sum(a)
      if(total == 0) return(NA)
      return(a[2] / total)
    })
    
    bias_fail <- !is.na(vaf) & (vaf < 0.25)
    df <- df[!bias_fail, ]
  }
  
  # 9. Blacklist filter
  if(nrow(df) > 0) {
    variant_id <- paste(df$CHROM, df$POS, df$REF, df$ALT, sep="-")
    df <- df[!variant_id %in% ensemble, ]
    
    var_chrompos <- paste(df$CHROM, df$POS, sep="-")
    pos_match    <- match(var_chrompos, bl_chrompos)
    has_match    <- !is.na(pos_match)
    
    ref_len_diff   <- nchar(df$REF[has_match]) != nchar(bl_ref[pos_match[has_match]])
    alt_len_diff   <- nchar(df$ALT[has_match]) != nchar(bl_alt[pos_match[has_match]])
    
    indel_mismatch <- has_match
    indel_mismatch[has_match] <- ref_len_diff | alt_len_diff
    df <- df[!indel_mismatch, ]
  }
  attr(df, "vcfFile") <- vcfFile
  return(df)
}

# VEP consequence terms from most to least severe
# (https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html)
.vepConsequenceOrder <- c(
  "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
  "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
  "transcript_amplification", "feature_elongation", "feature_truncation",
  "inframe_insertion", "inframe_deletion", "missense_variant",
  "protein_altering_variant", "splice_donor_5th_base_variant",
  "splice_region_variant", "splice_donor_region_variant",
  "splice_polypyrimidine_tract_variant", "incomplete_terminal_codon_variant",
  "start_retained_variant", "stop_retained_variant", "synonymous_variant",
  "coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant",
  "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
  "NMD_transcript_variant", "non_coding_transcript_variant",
  "coding_transcript_variant", "upstream_gene_variant",
  "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification",
  "TF_binding_site_variant", "regulatory_region_ablation",
  "regulatory_region_amplification", "regulatory_region_variant",
  "intergenic_variant", "sequence_variant"
)

# One CSQ entry per variant and gene.
# info, ref, alt: one element per variant (after bcftools norm -m -any).
# meta: CSQ field names from the VCF header.
# Returns list(row = index of the variant for each output row,
#              info = INFO with CSQ reduced to the chosen entry).
.csqPerGene <- function(info, meta, ref, alt) {
  n    <- length(info)
  loc  <- stri_locate_first_regex(info, "CSQ=[^;\t]*")
  has  <- which(!is.na(loc[, 1]))
  if(length(has) == 0) return(list(row = seq_len(n), info = info))
  
  csq     <- stri_sub(info[has], loc[has, 1] + 4L, loc[has, 2])
  entries <- stri_split_fixed(csq, ",")
  row     <- rep.int(has, lengths(entries))
  flat    <- unlist(entries, use.names = FALSE)
  fields  <- stri_split_fixed(flat, "|")
  
  getF <- function(name) {
    p <- match(name, meta)
    if(is.na(p)) return(rep("", length(flat)))
    v <- vapply(fields, `[`, "", p)
    v[is.na(v)] <- ""
    v
  }
  allele  <- getF("Allele")
  conseq  <- getF("Consequence")
  symbol  <- getF("SYMBOL")
  gene    <- getF("Gene")
  biotype <- getF("BIOTYPE")
  mane    <- getF("MANE_SELECT")
  canon   <- getF("CANONICAL")
  
  # Entries of the row's own allele. VEP writes the allele without the base
  # shared by all alleles of the original record ("-" if empty); when they do
  # not share it (e.g. A>T,AT) it keeps the full ALT. If nothing matches,
  # all entries are kept.
  strip   <- substr(ref, 1, 1) == substr(alt, 1, 1) & (nchar(ref) > 1 | nchar(alt) > 1)
  alt_vep <- ifelse(strip, substring(alt, 2), alt)
  alt_vep[alt_vep == ""] <- "-"
  exact    <- allele == alt[row]
  stripped <- allele == alt_vep[row]
  anyRow   <- function(x) (rowsum(as.integer(x), row, reorder = FALSE)[, 1] > 0)[match(row, unique(row))]
  allele_ok <- ifelse(anyRow(exact), exact, ifelse(anyRow(stripped), stripped, TRUE))
  
  # Most severe term of each entry (terms joined with "&")
  terms     <- stri_split_fixed(conseq, "&")
  cons_rank <- vapply(terms, function(t) min(match(t, .vepConsequenceOrder,
                                                   nomatch = length(.vepConsequenceOrder) + 1L)), 1L)
  
  gene_key <- ifelse(symbol != "", symbol, gene)
  keep     <- which(allele_ok)
  ord <- keep[order(row[keep],
                    biotype[keep] != "protein_coding",
                    cons_rank[keep],
                    mane[keep] == "",
                    canon[keep] != "YES",
                    !startsWith(gene[keep], "ENSG"),
                    keep)]
  best <- ord[!duplicated(paste(row[ord], gene_key[ord], sep = "\r"))]
  
  # Rebuild INFO for each chosen entry; variants without CSQ stay as they are
  i        <- row[best]
  new_info <- paste0(stri_sub(info[i], 1L, loc[i, 1] - 1L), "CSQ=", flat[best],
                     stri_sub(info[i], loc[i, 2] + 1L))
  no_csq   <- setdiff(seq_len(n), has)
  out_row  <- c(i, no_csq)
  out_info <- c(new_info, info[no_csq])
  o <- order(out_row, seq_along(out_row))
  list(row = out_row[o], info = out_info[o])
}