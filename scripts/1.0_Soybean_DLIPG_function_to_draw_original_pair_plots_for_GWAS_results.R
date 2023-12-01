##############################################################################################
######  Title: 1.0_Soybean_DLIPG_function_to_draw_original_pair_plots_for_GWAS_results    ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                                ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo       ######
######  Date: 2021/12/16 (Created), 2021/12/16 (Last Updated)                           ######
##############################################################################################




pairsForScores <- function(scoresAll,
                           genoMap = NULL,
                           thresholds = NULL,
                           sig.level = 0.05,
                           method = "BH",
                           verbose = TRUE) {
  allMissing <- apply(X = scoresAll,
                      MARGIN = 2,
                      FUN = function(x) {
                        all(is.na(x))
                      })
  scoresAll <- scoresAll[, !allMissing]
  
  if (is.null(thresholds)) {
    if (verbose) {
      print("Now computing thresholds...")
    }
    if (is.null(genoMap)) {
      stop("You must specify `genoMap` when you don't specify `thresholds`!!")
    }
    
    thresholds <- apply(
      X = scoresAll,
      MARGIN = 2,
      FUN = function(x, genoMap, sig.level, method) {
        input <- cbind(genoMap, x)
        thres <- try(RAINBOWR::CalcThreshold(input = input,
                                             sig.level = sig.level,
                                             method = method),
                     silent = TRUE)
        if (class(thres) %in% "try-error") {
          thres <- NA
        }
        
        return(thres)
      },
      genoMap = genoMap,
      sig.level = sig.level,
      method = method
    )
  } else {
    stopifnot(length(thresholds) == ncol(scoresAll))
  }
  
  
  rangeAll <- range(scoresAll, na.rm = TRUE)
  nMethods <- ncol(scoresAll)
  methodNames <- colnames(scoresAll)
  
  
  op <- par(mfrow = c(nMethods, nMethods),
            mai = rep(0.05, 4),
            omi = rep(0.35, 4))
  
  pb <- utils::txtProgressBar(min = 1,
                              max = nMethods ^ 2,
                              style = 3)
  for (i in 1:nMethods) {
    for (j in 1:nMethods) {
      if (verbose) {
        utils::setTxtProgressBar(pb = pb, value = (i - 1) * nMethods + j)
      }
      
      xaxt <- ifelse(i %in% nMethods, "s", "n")
      yaxt <- ifelse(j %in% 1, "s", "n")
      
      if (i == j) {
        histRes <- hist(
          x = scoresAll[, i],
          main = "",
          cex.main = 1,
          xlab = "",
          ylab = "",
          xlim = rangeAll,
          xaxt = "n",
          yaxt = "n"
        )
        text(
          x = mean(rangeAll),
          y = mean(range(histRes$counts)),
          labels = methodNames[i],
          cex = 2.3,
          font = 2
        )
      } else if (i < j) {
        plot(
          x = scoresAll[, i],
          y = scoresAll[, j],
          xlim = rangeAll,
          ylim = rangeAll,
          xlab = "",
          ylab = "",
          axes = TRUE,
          xaxt = xaxt,
          yaxt = yaxt,
          main = "",
          type = "n"
        )
        if (i == 1) {
          axis(side = 3)
        }
        if (j == nMethods) {
          axis(side = 4)
        }
        
        text(
          x = mean(rangeAll),
          y = mean(rangeAll),
          labels = round(cor(scoresAll[, j],
                             scoresAll[, i],
                             use = "complete.obs"), 3),
          cex = 2.3
        )
      } else {
        plot(
          x = scoresAll[, j],
          y = scoresAll[, i],
          xlim = rangeAll,
          ylim = rangeAll,
          xlab = "",
          ylab = "",
          axes = TRUE,
          xaxt = xaxt,
          yaxt = yaxt,
          main = ""
        )
        abline(0, 1, col = 2, lty = 2)
        abline(h = thresholds[i], col = 3, lty = 2)
        abline(v = thresholds[j], col = 3, lty = 2)
      }
      
      
    }
  }
  
  
  par(mfrow = c(1, 1),
      mai = c(1.02, 0.82, 0.82, 0.42))
  cat("\n")
}
