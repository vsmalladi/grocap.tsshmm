
#' Take a vector of read counts (in 10bp steps) and produce the (Y,X) pair
#' needed by the HMM
#'
#' Peaks identification makes use of steps in the local region marked as depleated. If a sufficient number is present, a step is called a peak if the TAP+ read count is above the mean TAP+ reads in neighboring background (depleated) steps.
#'
#' @param reads.in.steps.tap reads in TAP+ condition
#' @param reads.in.steps.notap reads in TAP- condition
#' @param scale.factor TAP- to TAP+ scale factor
#' @param w number of steps to look for background values (pos +/- w)
#' @param k minimum number of background steps to run log2 ratio test
#' @param log2.thresh a step is peaked if ratio to background windows is above this value
#' @return data suitable for the Naive HMMs
#' @export
sequence.to.data <- function(reads.in.steps.tap, reads.in.steps.notap, scale.factor, w = 10, k = 2, log2.thresh = 1) {
  reads.tap = abs(reads.in.steps.tap)
  reads.notap = abs(reads.in.steps.notap) * scale.factor

  ratio = log2(reads.tap / reads.notap)
  
  useful = reads.tap > 0 & reads.notap > 0
  increased = reads.tap > reads.notap

  L = length(reads.tap)

  # create X
  #
  X = rep(1, L) # default no data
  X[increased] = 2 # enriched
  X[useful & ratio < 0] = 3 # depleated (but > 0)

  # create Y
  #
  Y = rep(0, L)
  if (any(increased & useful)) {
    Y[increased & useful] = sapply(which(increased & useful), function(idx) {
      i = max(idx - w, 1)
      j = min(idx + w, L)
      bck.idxs = which(X[i:j] == 3)
      if (length(bck.idxs) >= k) {
        bck = mean(reads.tap[i:j][bck.idxs])
        val = log2(reads.tap[idx] / bck)
        
        if (val > log2.thresh)
          return(1)
      }
      return(0)
    })
  }

  return(rbind(Y, X))
}


#' Write a BED data frame as a prediction track (with colors)
#'
#' @param bed BED data frame
#' @param filename output BED filename
#' @param track track/description string
#' @param color.plus RGB color (comma separated) for plus strand
#' @param color.minus RGB color (comma separated) for minus strand
#' @export
write.track.bed <- function(bed, filename, track="hmm.preds", color.plus = "197,0,11", color.minus = "0,132,209") {
  colors = rep(color.plus, dim(bed)[1])
  colors[bed[,6] == '-'] = color.minus

  bed = cbind(bed[,1], as.integer(bed[,2]), as.integer(bed[,3]), bed[,4:6],
    as.integer(bed[,2]), as.integer(bed[,3]), colors)

  fout = file(filename, "w")
  cat("track name=", track, " description=\"", track, "\"  itemRgb=On\n", sep='', file=fout)
  write.table(bed, file=fout, col.names=F, row.names=F, quote=F, sep='\t')
  close(fout)
}
