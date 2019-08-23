#' @export
#'
#'
bb_readgtf <- function(gtf){

  ## Read gtf
  gtf <- rtracklayer::import(gtf)
  gtf_df <- as.data.frame(gtf)


  ## Split into genes
  by_genes <- split(gtf_df, gtf_df$gene_name)
  #list2env(by_genes)

  #for each one gene:

  plus_genes <- one_gene[which(one_gene$strand == "+"),]
  minus_genes <- one_gene[which(one_gene$strand == "-"),]
  plus_exons <- plus_genes[which(plus_genes$type == "exon"),]
  minus_exons <-  minus_genes[which(minus_genes$type == "exon"),]

  #plus_UTRS <- plus_genes[which(plus_genes$type == "UTR"),]

  plus_new <- data.frame()

  read_exon <- function(df, new_df){

    exon_seq <- df[1]
    exon_start <- df[2]
    exon_end <- df[3]

    exon <- paste(exon_seq,exon_start,exon_end)


  }




}
