#' Makes a bb_assembly object for alternate TxDB, OrgDb,and BSgenome packages
#'
#' @param Genome string of genome assembly name
#' @param TxDb string of the desired TxDb package
#' @param OrgDb string of the desired OrgDb package
#' @param gene.id.column string name of OrgDb column that refers to the given TxDb gene ID's
#' @param display.column string name of OrgDb column that is the type of gene symbol to be displayed in plots
#' @param BSgenome string of desired BSgenome package
#'
#' @export

bb_assembly <- function(Genome, TxDb, OrgDb, gene.id.column = "ENTREZID", display.column = "SYMBOL", BSgenome = NULL){

  object <- structure(list(Genome, TxDb, OrgDb, gene.id.column, display.column), class = "bb_assembly")
  if(!is.null(BSgenome)){
    object$BSgenome <- BSgenome
  }

  return(object)

}
