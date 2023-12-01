##################################################################################################
######  Title: 1.4.2_Soybean_DLIPG_function_to_drop_QTLs_from_GWAS_results_haplotype_based  ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                                    ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo           ######
######  Date: 2022/01/27 (Created), 2022/01/28 (Last Updated)                               ######
##################################################################################################




dropQtlFromData <- function(genoScore,
                            genoMap,
                            resGWAS,
                            phenoSimInfo,
                            idsAroundMrcList = NULL,
                            geneSet = NULL,
                            selectQtlMethod = "GWAS",
                            lenLD = 150000,
                            corThres = 0.35,
                            windowSize = 0) {
  if (is.null(geneSet)) {
    selectQtlMethod <- "LD"
  }

  traitName <- colnames(resGWAS)[4]

  qtlNamesListAll <- phenoSimInfo$mrkNamesList$Qtls
  names(qtlNamesListAll) <- c("Qtls", "HaploQtls")
  qtlNamesListQtls <- qtlNamesListAll$Qtls
  qtlNamesListHaploQtls <- qtlNamesListAll$HaploQtls
  qtlNamesListQtls <- lapply(X = qtlNamesListQtls,
                             FUN = function(qtlNames) {
                               if (!is.null(qtlNames)) {
                                 qtlNames <- stringr::str_replace(string = qtlNames,
                                                                  pattern = "X.", replacement = "X-")
                               }

                               return(qtlNames)
                             })
  qtlNamesListHaploQtls <- lapply(X = qtlNamesListHaploQtls,
                                  FUN = function(qtlNames) {
                                    if (!is.null(qtlNames)) {
                                      qtlNames <- stringr::str_replace(string = qtlNames,
                                                                       pattern = "X.", replacement = "X-")
                                    }

                                    return(qtlNames)
                                  })
  qtlNamesListAll <- list(Qtls = qtlNamesListQtls,
                          HaploQtls = qtlNamesListHaploQtls)



  qtlNamesList <- unlist(x = qtlNamesListAll, recursive = FALSE)
  qtlNamesList <- qtlNamesList[!sapply(X = qtlNamesList, FUN = is.null)]

  qtlNamesList <- lapply(X = qtlNamesList,
                         FUN = function(qtlNames) {
                           stringr::str_replace(string = qtlNames,
                                                pattern = "X.", replacement = "X-")
                         })
  qtlNamesAll <- unlist(qtlNamesList)

  nQtlType <- length(qtlNamesList)
  qtlTypeNames <- names(qtlNamesList)

  nQtlTypeWithTotal <- nQtlType + 1
  qtlTypeNamesWithTotal <- c("Total", qtlTypeNames)

  nTrueQtls <- sapply(X = qtlNamesList, FUN = length,
                      simplify = TRUE)
  nTrueQtlTotal <- sum(nTrueQtls)

  if (is.null(geneSet)) {
    genoMap <- genoMap[genoMap$marker %in% resGWAS$marker, , drop = FALSE]
  }

  qtlCandsList <- lapply(X = qtlNamesList,
                         FUN = function(x) {
                           match(x, genoMap$marker)
                         })
  qtlCandsTotal <- unlist(qtlCandsList)
  qtlCandsListWithTotal <- c(list(Total = qtlCandsTotal),
                             qtlCandsList)

  marker <- as.character(genoMap[, 1])
  chr <- genoMap[, 2]
  pos <- genoMap[, 3]
  nMarkers <- ncol(genoScore)

  cumPos <- RAINBOWR::cumsumPos(map = genoMap[, 1:3])

  positionQtlsList <- lapply(X = qtlCandsList,
                             FUN = function(qtlCands) {
                               cumPos[qtlCands]
                             })



  ###### 1. For normal Results ######
  resGWAS <- resGWAS[order(resGWAS[, 2], resGWAS[, 3]),]

  qtlCandsReplacedList <- sapply(
    X = 1:nQtlType,
    FUN = function(qtlTypeNo) {
      qtlCands <- qtlCandsList[[qtlTypeNo]]
      positionQtls <- positionQtlsList[[qtlTypeNo]]
      nTrueQtl <- nTrueQtls[qtlTypeNo]

      qtlCandsReplaced <- sapply(X = 1:nTrueQtl,
                                 FUN = function(trueQtlNo) {
                                   qtlCand <- qtlCands[trueQtlNo]
                                   positionQtl <- positionQtls[trueQtlNo]

                                   idsLdBlockRaw <- which(cumPos >= (positionQtl - lenLD) &
                                                            cumPos <= (positionQtl + lenLD))
                                   if (length(idsLdBlockRaw) >= 2) {
                                     corLDBlockRaw <- (cor(genoScore[, qtlCand],
                                                           genoScore[, idsLdBlockRaw]) ^ 2)[1, ]
                                   } else {
                                     corLDBlockRaw <- cor(genoScore[, qtlCand],
                                                          genoScore[, idsLdBlockRaw]) ^ 2
                                   }

                                   startBlockId <- max(min(which(x = corLDBlockRaw >= corThres))
                                                       - windowSize, 1)
                                   endBlockId <- min(max(which(x = corLDBlockRaw >= corThres))
                                                     + windowSize, length(idsLdBlockRaw))
                                   idsLdBlock <- idsLdBlockRaw[startBlockId:endBlockId]

                                   if (selectQtlMethod == "GWAS") {
                                     idsLdBlockLogpOrd <- idsLdBlock[order(resGWAS[idsLdBlock, 4], decreasing = TRUE)]
                                     qtlCandReplaced <- (idsLdBlockLogpOrd[!(idsLdBlockLogpOrd %in% qtlCand)])[1]
                                   } else if (selectQtlMethod == "LD")  {
                                     idsLdBlockCorOrd <- idsLdBlock[order(corLDBlockRaw[startBlockId:endBlockId], decreasing = TRUE)]
                                     qtlCandReplaced <- (idsLdBlockCorOrd[!(idsLdBlockCorOrd %in% qtlCand)])[1]
                                   }

                                   return(qtlCandReplaced)
                                 }, simplify = TRUE)

      return(qtlCandsReplaced)
    },
    simplify = FALSE
  )
  names(qtlCandsReplacedList) <- qtlTypeNames


  qtlNamesReplacedList <- lapply(X = qtlCandsReplacedList,
                                 FUN = function(qtlCands) {
                                   qtlNamesReplaced <- marker[qtlCands]

                                   return(qtlNamesReplaced)
                                 })

  qtlNamesReplacedList2 <- lapply(X = qtlNamesList,
                                  FUN = function(qtlNames) {
                                    sapply(X = qtlNames,
                                           FUN = function(qtlName) {
                                             whichMarker <- which(marker %in% qtlName)
                                             chrMarker <- chr[whichMarker]
                                             posMarker <- pos[whichMarker]
                                             whichMarkerInsteadCands <- c(whichMarker - 1, whichMarker + 1)
                                             finishLoop <- FALSE
                                             while (!finishLoop) {
                                               if (all(chrMarker == genoMap[whichMarkerInsteadCands, 2])) {
                                                 # whichNearest <- which.min(abs(genoMap[whichMarkerInsteadCands, 3] - posMarker))
                                                 if (selectQtlMethod == "GWAS") {
                                                   whichNearest <- which.max(resGWAS[whichMarkerInsteadCands, 4])
                                                 } else if (selectQtlMethod == "LD") {
                                                   whichNearest <- which.max(abs(cor(genoScore[, whichMarker],
                                                                                     genoScore[, whichMarkerInsteadCands])))
                                                 }
                                               } else {
                                                 whichNearest <- which(chrMarker == genoMap[whichMarkerInsteadCands, 2])
                                               }
                                               whichMarkerInstead <- whichMarkerInsteadCands[whichNearest]
                                               nearQtlName <- marker[whichMarkerInstead]
                                               finishLoop <- !(nearQtlName %in% qtlNamesAll)
                                               whichMarkerInsteadCands[whichNearest] <- whichMarkerInsteadCands[whichNearest] +
                                                 c(-1, 1)[whichNearest]
                                             }
                                             return(nearQtlName)
                                           })
                                  })
  for (qtlTypeNo in 1:nQtlType) {
    qtlNamesReplaced <- qtlNamesReplacedList[[qtlTypeNo]]
    if (any(is.na(qtlNamesReplaced))) {
      qtlNamesReplaced[is.na(qtlNamesReplaced)] <- qtlNamesReplacedList2[[qtlTypeNo]][is.na(qtlNamesReplaced)]
    }

    qtlNamesReplacedList[[qtlTypeNo]] <- qtlNamesReplaced
  }



  if (is.null(geneSet)) {
    if (!is.null(idsAroundMrcList)) {
      whichQtls <- sort(which(marker %in% qtlNamesAll))
      markerNos <- 1:nrow(genoMap)
      markerNos[whichQtls] <- NA
      for (whichQtl in whichQtls) {
        if (whichQtl != nrow(genoMap)) {
          markerNos[(whichQtl + 1):nrow(genoMap)] <- markerNos[(whichQtl + 1):nrow(genoMap)] - 1
        }
      }
      idsAroundMrcList <- lapply(X = idsAroundMrcList,
                                 FUN = function(x) na.omit(markerNos[x]))
      idsAroundMrcList <- idsAroundMrcList[unlist(lapply(X = idsAroundMrcList,
                                                         FUN = length)) != 0]
    }

    genoScore <- genoScore[, !(colnames(genoScore) %in% qtlNamesAll)]
    genoMap <- genoMap[!(genoMap$marker %in% qtlNamesAll), ]
    resGWAS <- resGWAS[!(resGWAS$marker %in% qtlNamesAll), ]


    qtlsOrHaploQtls <- unlist(lapply(
      X = stringr::str_split(string = qtlTypeNames, pattern = "\\."),
      FUN = function(x) x[[1]]
    ))

    qtlTypeNamesWoHaplo <- unlist(lapply(
      X = stringr::str_split(string = qtlTypeNames, pattern = "\\."),
      FUN = function(x) x[[2]]
    ))

    names(qtlNamesReplacedList) <- qtlTypeNamesWoHaplo

    phenoSimInfo$mrkNamesList$Qtls <- split(x = qtlNamesReplacedList,
                                            f = qtlsOrHaploQtls)[unique(qtlsOrHaploQtls)]
  } else {

    qtlHaploNamesList <- lapply(X = qtlNamesList,
                                FUN = function(qtlNames) {
                                  sapply(X = qtlNames, FUN = defineHaploBlock,
                                         haploBlockListDf = geneSet,
                                         genoMap = genoMap,
                                         resGWAS = resGWAS)
                                })
    qtlHaploNamesAll <- unlist(qtlHaploNamesList)




    oneMarkerInHaploBlockList <- lapply(X = qtlHaploNamesList,
                                        FUN = function(qtlHaploNames) {
                                          sapply(X = qtlHaploNames,
                                                 FUN = function(qtlHaploName) {
                                                   sum(geneSet[, 1] %in% qtlHaploName) == 1
                                                 })
                                        })
    qtlNamesForHaploReplacedList <- qtlNamesList
    for (qtlTypeNo in 1:nQtlType) {
      oneMarkerInHaplo <- oneMarkerInHaploBlockList[[qtlTypeNo]]
      if (any(oneMarkerInHaplo)) {
        qtlNamesForHaploReplaced <- qtlNamesForHaploReplacedList[[qtlTypeNo]]
        qtlNamesReplaced <- qtlNamesReplacedList[[qtlTypeNo]]

        qtlNamesForHaploReplaced[oneMarkerInHaplo] <- qtlNamesReplaced[oneMarkerInHaplo]
        qtlNamesForHaploReplacedList[[qtlTypeNo]] <- qtlNamesForHaploReplaced
      }
    }

    qtlsOrHaploQtls <- unlist(lapply(
      X = stringr::str_split(string = qtlTypeNames, pattern = "\\."),
      FUN = function(x) x[[1]]
    ))

    qtlTypeNamesWoHaplo <- unlist(lapply(
      X = stringr::str_split(string = qtlTypeNames, pattern = "\\."),
      FUN = function(x) x[[2]]
    ))

    names(qtlNamesForHaploReplacedList) <- qtlTypeNamesWoHaplo

    phenoSimInfo$mrkNamesList$Qtls <- split(x = qtlNamesForHaploReplacedList,
                                            f = qtlsOrHaploQtls)[unique(qtlsOrHaploQtls)]
  }


  return(list(genoScore = genoScore,
              genoMap = genoMap,
              resGWAS = resGWAS,
              phenoSimInfo = phenoSimInfo,
              idsAroundMrcList = idsAroundMrcList))
}
