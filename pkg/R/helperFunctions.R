#' Get the CODs assigned to each individual in the data
#'
#' @param x An fitted object returned from openVA::codeVA()
#' @param top An integer specifying the number of top causes included
#' in the results.
#'
#' @return A matrix with the death ID, top causes, and associated
#' probability for the top causes.
#'
#' @export
#'
indivCOD <- function (x, top = 3) {

  probs <- openVA::getIndivProb(x)
  cods  <- colnames(probs)

  out <- matrix(NA, nrow = nrow(probs), ncol = top*2)

  for (i in 1:nrow(probs)) {
    probsOrdered <- order(probs[i,], decreasing = TRUE)
    newTop <- top
    if(length(unique(probsOrdered)) <= top) newTop = (top - 1)
    if(newTop < 1){
      cat("Error: not enough unique causes \n")
      next
    }
    for (j in 1:newTop) {
      k <- seq(1, top*2, by = 2)[j]
      out[i, k  ] <- cods[probsOrdered[j]]
      out[i, k+1] <- round(probs[i, probsOrdered[j]], 4)
    }
  }
  out <- cbind(rownames(probs), out)
  colnames(out) <- c("ID", "Most Likely Cause", "Probability",
                     "Second Most Likely Cause", "Probability",
                     "Third Most Likely Cause", "Probability")
  out
}

#' Create CSMF plot with Undetermined for InterVA4 Output
#'
#' @param x An fitted object returned from openVA::codeVA()
#' @param top An integer specifying the number of top causes included
#' in the results.
#'
#' @return
#' @export
#'
CSMF2 <- function (x, top) {

  csmf  <- InterVA4::CSMF(x, InterVA.rule=TRUE, noplot=TRUE, min.prob=.00001)
  csmf2 <- sort(csmf, decreasing=TRUE)[1:top]
  isUndet <- "Undetermined" %in% names(csmf2)
  if(isUndet){
    idUndet <- which(names(csmf2) == "Undetermined")
    newBars <- c(csmf2[-idUndet], csmf2[idUndet])
    par(las = 2)
    par(mar = c(5, 15, 4, 2))
    barplot(newBars[top:1], horiz=TRUE, names.arg=names(newBars)[top:1],
            cex.names = 0.8, xlab = "Probability",
            col = rev(grey.colors(length(newBars))))
  } else{
    InterVA4::CSMF(x, InterVA.rule=TRUE, min.prob=.00001)
  }
}
