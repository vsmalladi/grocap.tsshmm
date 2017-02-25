#
# BigWig Datasets
#

#' Load the main GROcap, GROseq, CAGE datasets
#'
#' @param basePath Root directory for bigWig files
#' @param cellLine Directory name for cell specific files
#' @return List of bigWig objects for the given celline
#' @export
mainBwSet.grocap <- function(basePath = "/usr/data/GROseq.parser/hg19/", cellLine = "k562") {
    #
  # load bigWigs
  path.merge <- function(...) {
    paste(basePath, cellLine, "/", ..., sep='')
  }

  file.alt <- function(fname) {
    fname1 = paste(fname, ".bigWig", sep='')
    fname2 = paste(fname, "Rep1.bigWig", sep='')
    if (file.exists(fname1))
      return(fname1)
    return(fname2)
  }

  loadOrNULL.bigWig <- function(path) {
    if (file.exists(path))
      return(load.bigWig(path))
    return(NULL)
  }
  
  #
  cCell = paste(toupper(substring(cellLine, 1, 1)), substring(cellLine, 2), sep='')

  
  list(  
       ## GROseq data
       #GROseq.plus = load.bigWig(path.merge("groseq/groseq_plus.bigWig")),
       #GROseq.minus = load.bigWig(path.merge("groseq/groseq_minus.bigWig")),
       
       # GROseq TSS data
       GROcap.plus = load.bigWig(path.merge("groseq_tss/groseq_tss_wTAP_plus.bigWig")),
       GROcap.minus = load.bigWig(path.merge("groseq_tss/groseq_tss_wTAP_minus.bigWig")),
       
      #  # CAGE Whole Cell TSS data
      #  CAGEWC.plus = load.bigWig(path.merge("cage/wgEncodeRikenCage", cCell, "CellPapPlusSignalRep1.bigWig")),
      #  CAGEWC.minus = load.bigWig(path.merge("cage/wgEncodeRikenCage", cCell, "CellPapMinusSignalRep1.bigWig")),

      #  # CAGE Nucleus TSS data PA+
      #  CAGENuc.plus = load.bigWig(path.merge("cage/wgEncodeRikenCage", cCell, "NucleusPapPlusSignalRep1.bigWig")),
      #  CAGENuc.minus = load.bigWig(path.merge("cage/wgEncodeRikenCage", cCell, "NucleusPapMinusSignalRep1.bigWig")),

      #  # CAGE Nucleus TSS data PA-
      #  CAGENucPAM.plus = load.bigWig(file.alt(path.merge("cage/wgEncodeRikenCage", cCell, "NucleusPamPlusSignal"))),
      #  CAGENucPAM.minus = load.bigWig(file.alt(path.merge("cage/wgEncodeRikenCage", cCell, "NucleusPamMinusSignal"))),
  
       
      #  # CAGE Chromatin TSS data
      #  CAGEChr.plus = loadOrNULL.bigWig(file.alt(path.merge("cage/wgEncodeRikenCage", cCell, "ChromatinTotalPlusSignal"))),
      #  CAGEChr.minus = loadOrNULL.bigWig(file.alt(path.merge("cage/wgEncodeRikenCage", cCell, "ChromatinTotalMinusSignal"))),

       # Revised CAGE data
       CAGENuc.bp.plus = load.bigWig(file.alt(path.merge("cage_bp/riken_cage_nucleus_pap_plus"))),
       CAGENuc.bp.minus = load.bigWig(file.alt(path.merge("cage_bp/riken_cage_nucleus_pap_minus")))
       )
}


#' Load the background GROcap datasets
#'
#' @param basePath Root directory for bigWig files
#' @param cellLine Directory name for cell specific files
#' @return List of bigWig objects for the given celline
#' @export
backBwSet.grocap <- function(basePath = "/usr/data/GROseq.parser/hg19/", cellLine = "k562") {
    #
  # load bigWigs
  path.merge <- function(...) {
    paste(basePath, cellLine, "/", ..., sep='')
  }
  
  #
  cCell = paste(toupper(substring(cellLine, 1, 1)), substring(cellLine, 2), sep='')

  
  list(  
       # GROseq TSS data
       GROcap.plus = load.bigWig(path.merge("groseq_tss/groseq_tss_noTAP_plus.bigWig")),
       GROcap.minus = load.bigWig(path.merge("groseq_tss/groseq_tss_noTAP_minus.bigWig"))
       
       )
}
