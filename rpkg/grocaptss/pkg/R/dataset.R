#
# BigWig Datasets
#

#' Load the main GROcap, GROseq, CAGE datasets
#'
#' @param basePath Root directory for bigWig files
#' @param cellLine Directory name for cell specific files
#' @return List of bigWig objects for the given celline
#' @export
mainBwSet.grocap <- function(basePath = "/Volumes/project/GCRB/Lee_Lab/s163035/Tim_MC7_E2/signal/", cellLine = "E2") {
    #
  # load bigWigs
  path.merge <- function(...) {
    paste(basePath, cellLine, ..., sep='')
  }

  list(
       # PRO-cap TAP data
       GROcap.plus = load.bigWig(path.merge("_PlusTAP.norm.positive.bw")),
       GROcap.minus = load.bigWig(path.merge("_PlusTAP.norm.negative.bw"))
       )
}


#' Load the background GROcap datasets
#'
#' @param basePath Root directory for bigWig files
#' @param cellLine Directory name for cell specific files
#' @return List of bigWig objects for the given celline
#' @export
backBwSet.grocap <- function(basePath = "/Volumes/project/GCRB/Lee_Lab/s163035/Tim_MC7_E2/signal/", cellLine = "E2") {
    #
  # load bigWigs
  path.merge <- function(...) {
    paste(basePath, cellLine, ..., sep='')
  }


  list(
       # PRO-cap minus TAP data
       GROcap.plus = load.bigWig(path.merge("_MinusTap.norm.positive.bw")),
       GROcap.minus = load.bigWig(path.merge("_MinusTap.norm.negative.bw"))

       )
}
