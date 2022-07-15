getExons <- function(assembly, chromosome, start, stop) {
    if (is(assembly$TxDb, "TxDb")) {
        tx_db <- assembly$TxDb
    } else {
        tx_db <- eval(parse(text = paste0(as.name(assembly$TxDb),
                                        "::",
                                        as.name(assembly$TxDb))))
    }

    if (length(assembly$OrgDb) == 1){
        
        ## Custom orgDb
        if (is(assembly$OrgDb, "OrgDb")){
            org_db <- assembly$OrgDb
        } else {
            org_db <- eval(parse(text = paste0(as.name(assembly$OrgDb),
                                               "::",
                                               as.name(assembly$OrgDb))))
        }
    }
    
    genes_on_chrom <- suppressMessages(GenomicFeatures::genes(tx_db,
        filter = list(tx_chrom = chromosome)
    )) ## Still GRanges
    genes_in_range <- genes_on_chrom[genes_on_chrom@ranges@start <= stop &
        genes_on_chrom@ranges@start +
            genes_on_chrom@ranges@width >= start]
    txs_wo_symbols <- suppressMessages(AnnotationDbi::select(tx_db,
        keys = genes_in_range$gene_id,
        columns = AnnotationDbi::columns(tx_db),
        keytype = "GENEID"
    ))
    ## Double check filtering of chromosome
    txs_wo_symbols <- txs_wo_symbols[which(txs_wo_symbols$TXCHROM==chromosome),]

    if (assembly$gene.id.column == assembly$display.column) {
        txs_w_symbols <- txs_wo_symbols
    } else if (assembly$gene.id.column == "ORF") {
        txs_wo_symbols <- cbind(txs_wo_symbols,
            "SHORTORF" = sub("-.+", "", txs_wo_symbols$GENEID)
        )
        gene_symbols <- suppressMessages(AnnotationDbi::select(org_db,
            keys = sub("-.+", "", genes_in_range$gene_id),
            columns = assembly$display.column,
            keytype = "ORF"
        ))
        txs_w_symbols <- merge(txs_wo_symbols, gene_symbols,
            by.x = "SHORTORF", by.y = "ORF"
        )
    } else {
        gene_symbols <- suppressMessages(AnnotationDbi::select(org_db,
            keys = genes_in_range$gene_id,
            columns = assembly$display.column,
            keytype = assembly$gene.id.column
        ))
        txs_w_symbols <- merge(txs_wo_symbols, gene_symbols,
            by.x = "GENEID", by.y = assembly$gene.id.column
        )
    }

    return(txs_w_symbols)
}
