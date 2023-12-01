#############################################################################################
######  Title: 1.3.1_Soybean_DLIPG_function_to_define_haplotype_block_for_each_marker    ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                               ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo      ######
######  Date: 2022/01/25 (Created), 2022/01/25 (Last Updated)                          ######
#############################################################################################




defineHaploBlock <- function(mrkName,
                             haploBlockListDf,
                             genoMap,
                             resGWAS = NULL) {
  haploBlockIncludeMrk <- haploBlockListDf[mrkName, 1]

  if (is.na(haploBlockIncludeMrk)) {
    whereInMap <- which(genoMap$marker %in% mrkName)
    chrMrk <- genoMap$chr[whereInMap]
    posMrk <- genoMap$pos[whereInMap]

    haploBlockIncludeMrkAdjacents <- rep(NA, 2)
    posAdjacentsInHaplo <- rep(NA, 2)
    finishDefine <- finishBackward <- finishForward <- FALSE
    distanceSNP <- 0

    while (!finishDefine) {
      distanceSNP <- distanceSNP + 1

      if (!finishBackward) {
        if (whereInMap - distanceSNP >= 1) {
          if (chrMrk == genoMap$chr[whereInMap - distanceSNP]) {
            haploBlockIncludeMrkBackward <-
              haploBlockListDf[genoMap$marker[whereInMap - distanceSNP], 1]
            finishBackward <- !is.na(haploBlockIncludeMrkBackward)
            if (finishBackward) {
              posAdjacentsInHaplo[1] <- genoMap$pos[whereInMap - distanceSNP]
              haploBlockIncludeMrkAdjacents[1] <- haploBlockIncludeMrkBackward
            }
          } else {
            finishBackward <- TRUE
          }
        } else {
          finishBackward <- TRUE
        }
      }


      if (!finishForward) {
        if (whereInMap + distanceSNP <= nrow(genoMap)) {
          if (chrMrk == genoMap$chr[whereInMap + distanceSNP]) {
            haploBlockIncludeMrkForward <-
              haploBlockListDf[genoMap$marker[whereInMap + distanceSNP], 1]
            finishForward <- !is.na(haploBlockIncludeMrkForward)
            if (finishForward) {
              posAdjacentsInHaplo[2] <- genoMap$pos[whereInMap + distanceSNP]
              haploBlockIncludeMrkAdjacents[2] <- haploBlockIncludeMrkForward
            }
          } else {
            finishForward <- TRUE
          }
        } else {
          finishForward <- TRUE
        }
      }

      finishDefine <- finishBackward & finishForward
    }

    if (any(is.na(haploBlockIncludeMrkAdjacents))) {
      haploBlockIncludeMrk <- haploBlockIncludeMrkAdjacents[!is.na(haploBlockIncludeMrkAdjacents)]
    } else {
      if (is.null(resGWAS)) {
        haploBlockIncludeMrk <- haploBlockIncludeMrkAdjacents[which.min(abs(posMrk - posAdjacentsInHaplo))]
      } else {
        haploBlockIncludeMrk <- haploBlockIncludeMrkAdjacents[which.max(resGWAS[resGWAS$marker %in% haploBlockIncludeMrkAdjacents, 4])]
      }
    }
  }


  return(haploBlockIncludeMrk)
}
