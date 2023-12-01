#############################################################################################
######  Title: 1.4_Soybean_DLIPG_function_to_summary_GWAS_results_haplotype_based      ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                               ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo      ######
######  Date: 2022/01/18 (Created), 2022/01/25 (Last Updated)                          ######
#############################################################################################


source("scripts/1.3.1_Soybean_DLIPG_function_to_define_haplotype_block_for_each_marker.R")

summaryGwasResults <- function(resGWAS,
                               genoScore,
                               genoMap,
                               phenoSimInfo,
                               geneSet = NULL,
                               nTopBlockFalsePositive = 10,
                               sigLevel = c(0.05, 0.01),
                               methodThres = "BH",
                               inflatorPlus = 2,
                               lenLD = 150000,
                               corThres = 0.35,
                               windowSize = 0,
                               saveName = NULL,
                               plotROC = TRUE,
                               returnAucRes = TRUE,
                               returnBlockRes = TRUE,
                               idsAroundMrcList = NULL,
                               nCores = 1) {
  traitName <- colnames(resGWAS)[4]

  qtlNamesListAll <- phenoSimInfo$mrkNamesList$Qtls
  # names(qtlNamesListAll) <- c("Qtls", "HaploQtls")
  qtlNamesListQtls <- qtlNamesListAll$Qtls
  qtlNamesListHaploQtls <- qtlNamesListAll$HaploQtls
  # qtlNamesListQtls <- lapply(X = qtlNamesListQtls,
  #                            FUN = function(qtlNames) {
  #                              if (!is.null(qtlNames)) {
  #                                qtlNames <- stringr::str_replace(string = qtlNames,
  #                                                                 pattern = "X.", replacement = "X-")
  #                              }
  #
  #                              return(qtlNames)
  #                            })
  # qtlNamesListHaploQtls <- lapply(X = qtlNamesListHaploQtls,
  #                                 FUN = function(qtlNames) {
  #                                   if (!is.null(qtlNames)) {
  #                                     qtlNames <- stringr::str_replace(string = qtlNames,
  #                                                                      pattern = "X.", replacement = "X-")
  #                                   }
  #
  #                                   return(qtlNames)
  #                                 })
  # qtlNamesListAll <- list(Qtls = qtlNamesListQtls,
  #                         HaploQtls = qtlNamesListHaploQtls)



  qtlNamesList <- unlist(x = qtlNamesListAll, recursive = FALSE)
  qtlNamesList <- qtlNamesList[!sapply(X = qtlNamesList, FUN = is.null)]


  # qtlNamesList <- lapply(X = qtlNamesList,
  #                        FUN = function(qtlNames) {
  #                          stringr::str_replace(string = qtlNames,
  #                                               pattern = "X.", replacement = "X-")
  #                        })

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

  nullChrNo <- phenoSimInfo$inputParams$nullChrNo


  if (is.null(geneSet)) {
    ###### 1. For normal Results ######
    ##### 1.1. Chromosome and position ######
    resGWAS <- resGWAS[order(resGWAS[, 2], resGWAS[, 3]),]



    ##### 1.2. Calculate theshold #####
    thresholds <- try(RAINBOWR::CalcThreshold(input = resGWAS,
                                              sig.level = sigLevel,
                                              method = methodThres), silent = TRUE)
    if ("try-error" %in% class(thresholds)) {
      thresholds <- rep(NA, length(sigLevel))
    }
    nOverThresholds <- sapply(X = thresholds,
                              FUN = function(threshold) {
                                ifelse(test = is.na(threshold),
                                       yes = 0,
                                       no = sum(resGWAS[, 4] >= threshold))
                              })


    ##### 1.3. The rank of -log10(p) of qtls #####
    logpOrd <- order(resGWAS[, 4], decreasing = TRUE)
    qtlLogpOrdsList <- lapply(X = qtlCandsList,
                              FUN = function(qtlCands) {
                                match(qtlCands, logpOrd)
                              })


    ldBlockResList <- sapply(
      X = 1:nQtlType,
      FUN = function(qtlTypeNo) {
        qtlCands <- qtlCandsList[[qtlTypeNo]]
        positionQtls <- positionQtlsList[[qtlTypeNo]]
        nTrueQtl <- nTrueQtls[qtlTypeNo]

        ldBlockRes <- sapply(X = 1:nTrueQtl,
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
                               idsLdBlockLogpOrd <- idsLdBlock[order(resGWAS[idsLdBlock, 4], decreasing = TRUE)]
                               overQtlInLdBlock <- idsLdBlockLogpOrd[
                                 1:match(qtlCand, idsLdBlockLogpOrd)
                               ]

                               return(list(idsLdBlock = idsLdBlock,
                                           overQtlInLdBlock = overQtlInLdBlock))
                             }, simplify = FALSE)

        return(ldBlockRes)
      },
      simplify = FALSE
    )
    names(ldBlockResList) <- qtlTypeNames
    idsLdBlockList <- lapply(X = ldBlockResList,
                             FUN = function(ldBlockResults) {
                               lapply(X = ldBlockResults,
                                      FUN = function(ldBlockRes) {
                                        ldBlockRes$idsLdBlock
                                      })
                             })

    overQtlInLdBlockList <- lapply(X = ldBlockResList,
                                   FUN = function(ldBlockResults) {
                                     lapply(X = ldBlockResults,
                                            FUN = function(ldBlockRes) {
                                              ldBlockRes$overQtlInLdBlock
                                            })
                                   })
    qtlLogpOrdsVec <- unlist(qtlLogpOrdsList)



    #### 1.4. Calculation of AUC and draw ROC curve ####
    if (returnAucRes) {
      overQtlInLdBlockListUnlist <- unlist(overQtlInLdBlockList, recursive = FALSE)
      overQtlInLdBlockListUnlistOrd <- overQtlInLdBlockListUnlist[order(qtlLogpOrdsVec)]
      nOverQtlOrdInLdBlocksCumsumTotal <- rep(NA, nTrueQtlTotal)
      overQtlOrdInLdBlocksVec <- NULL
      for (trueQtlNoTotal in 1:nTrueQtlTotal) {
        overQtlOrdInLdBlocksVec <- c(overQtlOrdInLdBlocksVec,
                                     overQtlInLdBlockListUnlistOrd[[trueQtlNoTotal]])
        nOverQtlOrdInLdBlocksCumsumTotal[order(qtlLogpOrdsVec)[trueQtlNoTotal]] <- length(unique(overQtlOrdInLdBlocksVec))
      }
      nOverQtlOrdInLdBlocksCumsumTotalSorted <- sort(nOverQtlOrdInLdBlocksCumsumTotal)

      nOverQtlOrdInLdBlocksCumsumList <- split(
        x = nOverQtlOrdInLdBlocksCumsumTotal,
        f = rep(qtlTypeNames, nTrueQtls)
      )[qtlTypeNames]
      nOverQtlOrdInLdBlocksCumsumSotedList <- lapply(
        X = nOverQtlOrdInLdBlocksCumsumList, FUN = sort
      )




      #### 1.4.1. Each result of AUC and draw ROC curve ####
      aucEachRes <- sapply(
        X = 1:nQtlType,
        FUN = function(qtlTypeNo) {
          qtlLogpOrds <- qtlLogpOrdsList[[qtlTypeNo]]
          nTrueQtl <- nTrueQtls[qtlTypeNo]
          nOverQtlOrdInLdBlocksCumsum <-
            nOverQtlOrdInLdBlocksCumsumSotedList[[qtlTypeNo]]

          aucTP <- c(0, 1:nTrueQtl, nTrueQtl) / nTrueQtl
          aucFP <- (c(0, sort(qtlLogpOrds), nMarkers) -
                      c(0, nOverQtlOrdInLdBlocksCumsum,
                        nOverQtlOrdInLdBlocksCumsum[nTrueQtl])) /
            (nMarkers - nOverQtlOrdInLdBlocksCumsum[nTrueQtl])

          aucFP[aucFP <= 0] <- 0


          auc <- aucTP[2] * aucFP[2] / 2

          for (m in 2:(length(aucFP) - 1)) {
            aucSegment <-
              (aucTP[m] + aucTP[m + 1]) * (aucFP[m + 1] - aucFP[m]) / 2
            auc <- auc + aucSegment
          }

          return(list(aucTP = aucTP,
                      aucFP = aucFP,
                      auc = auc))
        }, simplify = FALSE
      )
      names(aucEachRes) <- qtlTypeNames
      aucEach <- sapply(X = aucEachRes,
                        FUN = function(x) x$auc,
                        simplify = TRUE)


      if (plotROC) {
        for (qtlTypeNo in 1:nQtlType) {
          qtlTypeName <- qtlTypeNames[qtlTypeNo]
          aucEachNow <- aucEach[qtlTypeNo]
          aucTP <- aucEachRes[[qtlTypeNo]]$aucTP
          aucFP <- aucEachRes[[qtlTypeNo]]$aucFP

          if (is.null(saveName)) {
            plot(aucFP, aucTP, type = "l",
                 main = qtlTypeName)
            text(0.83, 0.15, paste("AUC = ", round(aucEachNow, 3), sep = ""), cex = 1.5)
          } else{
            png(paste(saveName, traitName, "_", qtlTypeName,
                      "_ROC_curve.png", sep = ""),
                width = 700)
            plot(aucFP, aucTP, type = "l",
                 main = qtlTypeName)
            text(0.83, 0.15, paste("AUC = ", round(aucEachNow, 3), sep = ""), cex = 1.5)
            dev.off()
          }
        }
      }



      #### 1.4.2. Total result of AUC and draw ROC curve ####
      aucTPTotal <- c(0, 1:nTrueQtlTotal, nTrueQtlTotal) / nTrueQtlTotal
      aucFPTotal <- (c(0, sort(qtlLogpOrdsVec), nMarkers) -
                       c(0, nOverQtlOrdInLdBlocksCumsumTotalSorted,
                         nOverQtlOrdInLdBlocksCumsumTotalSorted[nTrueQtlTotal])) /
        (nMarkers - nOverQtlOrdInLdBlocksCumsumTotalSorted[nTrueQtlTotal])

      aucFPTotal[aucFPTotal <= 0] <- 0

      aucTotal <- aucTPTotal[2] * aucFPTotal[2] / 2
      for (m in 2:(length(aucFPTotal) - 1)) {
        aucSegmentTotal <-
          (aucTPTotal[m] + aucTPTotal[m + 1]) *
          (aucFPTotal[m + 1] - aucFPTotal[m]) / 2
        aucTotal <- aucTotal + aucSegmentTotal
      }

      if (plotROC) {
        if (is.null(saveName)) {
          plot(aucFPTotal, aucTPTotal, type = "l",
               main = "Total")
          text(0.83, 0.15, paste("AUC = ", round(aucTotal, 3), sep = ""), cex = 1.5)
        } else{
          png(paste(saveName, traitName, "_Total",
                    "_ROC_curve.png", sep = ""),
              width = 700)
          plot(aucFPTotal, aucTPTotal, type = "l",
               main = "Total")
          text(0.83, 0.15, paste("AUC = ", round(aucTotal, 3), sep = ""), cex = 1.5)
          dev.off()
        }
      }
      names(aucTotal) <- "Total"

      auc <- c(aucTotal, aucEach)
    } else {
      auc <- NULL
    }




    #### 1.5. Calculate FPR, FNR, FDR and Harmonic mean of precision and recall ####
    ### 1.5.1. Markers over thresholds ###
    if (returnBlockRes) {
      overThresholdVals <- sapply(X = thresholds,
                                  FUN = function(threshold) {
                                    resGWAS[which(resGWAS[, 4] >= threshold), 4]
                                  }, simplify = FALSE)

      overThresholdCandsList <- sapply(X = thresholds,
                                       FUN = function(threshold) {
                                         which(resGWAS[, 4] >= threshold)
                                       }, simplify = FALSE)

      # overThresholdCandsLogpOrdList <-
      #   sapply(X = nOverThresholds,
      #          FUN = function(nOverThreshold) {
      #            if (nOverThreshold != 0) {
      #              overThresholdCandsLogpOrd <- logpOrd[1:nOverThreshold]
      #            } else {
      #              overThresholdCandsLogpOrd <- NULL
      #            }
      #            return(overThresholdCandsLogpOrd)
      #          }, simplify = FALSE)



      ### 1.5.2. Preparation for computation of number of TP, FP, FN, and TN ###
      # overThresholdCandsLogpOrdListForHm <- overThresholdCandsLogpOrdList
      logpOrdForHm <- logpOrd




      ### 1.5.3. Start computation of number of TP, FP, FN, and TN ###
      if (is.null(idsAroundMrcList)) {
        while (length(logpOrdForHm) >= 1) {
          ## 1.5.3.1. Detect LD block for the marker of interest ##
          idMrcNow <- logpOrdForHm[1]
          idsAroundMrcRaw <- which(cumPos >= (cumPos[idMrcNow] - lenLD) &
                                     cumPos <= (cumPos[idMrcNow] + lenLD))
          corAroundMrcRaw <- (cor(genoScore[, idMrcNow],
                                  genoScore[, idsAroundMrcRaw, drop = FALSE]) ^ 2)[1, ]

          startBlockId <- max(min(which(x = corAroundMrcRaw >= corThres))
                              - windowSize, 1)
          endBlockId <- min(max(which(x = corAroundMrcRaw >= corThres))
                            + windowSize, length(idsAroundMrcRaw))
          idsAroundMrc <- idsAroundMrcRaw[startBlockId:endBlockId]
          idsAroundMrc <- idsAroundMrc[idsAroundMrc %in% logpOrdForHm]
          idsAroundMrcListNow <- list(idsAroundMrc)
          names(idsAroundMrcListNow) <- idMrcNow
          idsAroundMrcList <- c(idsAroundMrcList, idsAroundMrcListNow)



          ## 1.5.3.1. Detect LD block for the marker of interest ##
          logpOrdForHm <- logpOrdForHm[!(logpOrdForHm %in% idsAroundMrc)]

        }
      } else {
        blockLogPs <- unlist(x = lapply(X = (idsAroundMrcList),
                                        FUN = function(x) {
                                          suppressWarnings(max(resGWAS[x, 4], na.rm = TRUE))
                                        }))
        idsAroundMrcList <- idsAroundMrcList[order(blockLogPs, decreasing = TRUE)]
      }


      whichQtlCandsList <- lapply(
        X = qtlCandsList,
        FUN = function(qtlCands) {
          unlist(parallel::mclapply(
            X = idsAroundMrcList,
            FUN = function(idsAroundMrc) {
              any(idsAroundMrc %in% qtlCands)
            },
            mc.cores = nCores
          ))
        }
      )


      whichQtlCandsListWithTotal <- c(
        list(Total = apply(X = do.call(what = cbind,
                                       args = whichQtlCandsList),
                           MARGIN = 1, FUN = any)),
        whichQtlCandsList
      )

      # whichQtlCandsTotalMat <- do.call(what = rbind,
      #                                  args = parallel::mclapply(
      #   X = idsAroundMrcList,
      #   FUN = function(idsAroundMrc) {
      #     sapply(X = qtlCandsListWithTotal$Total,
      #            FUN = function(qtlCands) {
      #              any(idsAroundMrc %in% qtlCands)
      #            })
      #   },
      #   mc.cores = nCores
      # ))
      #
      # whichNoQtlCandsTotal <- apply(X = whichQtlCandsTotalMat,
      #                             MARGIN = 2, FUN = which)
      # names(whichNoQtlCandsTotal) <- rownames(whichQtlCandsTotalMat)[whichNoQtlCandsTotal]
      # whichNoQtlCandsList <- c(list(Total = whichNoQtlCandsTotal),
      #                        split(x = whichNoQtlCandsTotal,
      #                              f = rep(qtlTypeNames, nTrueQtls)))
      # names(whichNoQtlCandsList) <- qtlTypeNamesWithTotal



      countQtlList <- lapply(X = whichQtlCandsListWithTotal,
                             FUN = sum)
      falseBlockLogps <- unlist(x = lapply(X = (idsAroundMrcList[!whichQtlCandsListWithTotal$Total]),
                                           FUN = function(x) {
                                             suppressWarnings(max(resGWAS[x, 4], na.rm = TRUE))
                                           }))
      falseBlockLogpsList <- rep(list(falseBlockLogps), nQtlTypeWithTotal)
      countFP <- sapply(X = thresholds,
                        FUN = function(threshold) {
                          if (!is.na(threshold)) {
                            return(sum(falseBlockLogps >= threshold))
                          } else {
                            return(0)
                          }
                        }, simplify = TRUE)
      countFPList <- rep(list(countFP), nQtlTypeWithTotal)

      names(countFPList) <- names(falseBlockLogpsList) <-
        qtlTypeNamesWithTotal

      countAllBlocksList <- lapply(X = countQtlList,
                                   FUN = function(countQtl) {
                                     length(idsAroundMrcList) -
                                       (countQtlList$Total - countQtl)
                                   })

      countNonQtlList <- lapply(X = countQtlList,
                                FUN = function(countQtl) {
                                  length(idsAroundMrcList) - countQtlList$Total
                                })

      # qtlCandsOverThresholdList <- lapply(
      #   X = qtlCandsListWithTotal,
      #   FUN = function(qtlCands) {
      #     sapply(X = thresholds,
      #            FUN = function(threshold) {
      #              if (!is.na(threshold)) {
      #                qtlCands[resGWAS[qtlCands, 4] >= threshold]
      #              } else {
      #                NULL
      #              }
      #            }, simplify = FALSE)
      #   }
      # )

      # countTPList <- lapply(
      #   X = qtlCandsOverThresholdList,
      #   FUN = function(qtlCandsOverThreshold) {
      #     unlist(lapply(X = qtlCandsOverThreshold,
      #                   FUN = function(qtlCandsOverThresholdNow) {
      #                     if (length(qtlCandsOverThresholdNow) >= 1) {
      #                       return(sum(unlist(parallel::mclapply(
      #                         X = idsAroundMrcList,
      #                         FUN = function(idsAroundMrc) {
      #                           any(idsAroundMrc %in% qtlCandsOverThresholdNow)
      #                         },
      #                         mc.cores = nCores
      #                       ))))
      #                     } else {
      #                       return(0)
      #                     }
      #                   }))
      #   }
      # )

      trueBlockLogpsList <- lapply(X = whichQtlCandsListWithTotal,
                                  FUN = function(whichQtlCandsEach) {
                                    unlist(x = lapply(X = (idsAroundMrcList[whichQtlCandsEach]),
                                                      FUN = function(x) {
                                                        suppressWarnings(max(resGWAS[x, 4], na.rm = TRUE))
                                                      }))
                                  })
      countTPList <- lapply(X = trueBlockLogpsList,
                            FUN = function(trueBlockLogps) {
                              sapply(X = thresholds,
                                     FUN = function(threshold) {
                                       if (!is.na(threshold)) {
                                         return(sum(trueBlockLogps >= threshold))
                                       } else {
                                         return(0)
                                       }
                                     }, simplify = TRUE)
                            })

      whichNoQtlCandsListWithTotal <- lapply(X = whichQtlCandsListWithTotal,
                                             FUN = which)

      countBlocksForAucList <- lapply(
        X = whichNoQtlCandsListWithTotal,
        FUN = function(whichNoQtlCands) {
          whichNoQtlCandsNames <- readr::parse_number(names(whichNoQtlCands))
          whichNoQtlCandsTotalNames <- readr::parse_number(names(whichNoQtlCandsListWithTotal$Total))
          countBlocksForAucUnSorted <- whichNoQtlCands -
            match(whichNoQtlCandsNames, whichNoQtlCandsTotalNames) +
            1:length(whichNoQtlCands)

          countBlocksForAuc <- countBlocksForAucUnSorted[rank(whichNoQtlCandsNames)]
        }
      )
      countBlocksForAucList$Total <- countBlocksForAucList$Total[
        unlist(sapply(X = countBlocksForAucList[-1],
                      FUN = names, simplify = FALSE))
      ]




      nTPList <- countTPList
      nFPList <- countFPList
      nFNList <- mapply(FUN = "-",
                        countQtlList,
                        countTPList,
                        SIMPLIFY = FALSE)
      nTNList <- mapply(
        FUN = function(countAllBlocks,
                       nTP, nFP, nFN) {
          countAllBlocks - nTP - nFP - nFN
        },
        countAllBlocks = countAllBlocksList,
        nTP = nTPList,
        nFP = nFPList,
        nFN = nFNList,
        SIMPLIFY = FALSE)
      names(nTNList) <- names(nFNList)

      summaryStatResList <- mapply(FUN = function(nTP, nFP, nFN, nTN) {
        fdr <- nFP / (nFP + nTP)
        fpr <- nFP / (nFP + nTN)
        fnr <- nFN / (nTP + nFN)
        recall <- nTP / (nTP + nFN)
        precision <- nTP / (nTP + nFP)
        accuracy <- (nTP + nTN) / (nTP + nFP + nTN + nFN)
        harmonicMean <- 2 * recall * precision / (recall + precision)

        if (any(is.nan(precision))) {
          fdr[is.nan(precision)] <- NA
          harmonicMean[is.nan(precision)] <- NA
          precision[is.nan(precision)] <- NA
        }

        if (any(recall == 0 | precision == 0)) {
          harmonicMean[recall == 0 | precision == 0] <- 0
        }

        return(list(
          fdr = fdr,
          fpr = fpr,
          fnr = fnr,
          recall = recall,
          precision = precision,
          accuracy = accuracy,
          harmonicMean = harmonicMean
        ))
      },
      nTP = nTPList,
      nFP = nFPList,
      nFN = nFNList,
      nTN = nTNList,
      SIMPLIFY = FALSE)
    } else {
      summaryStatResList <- rep(list(list(
        fdr = NULL,
        fpr = NULL,
        fnr = NULL,
        recall = NULL,
        precision = NULL,
        accuracy = NULL,
        harmonicMean = NULL
      )), nQtlTypeWithTotal)
      names(summaryStatResList) <- qtlTypeNamesWithTotal
    }







    #### 1.6. LD-based AUC ####
    if (returnAucRes) {
      aucLdBaseRes <- sapply(X = 1:(nQtlTypeWithTotal),
                             FUN = function(qtlTypeNo) {
                               countQtl <- countQtlList[[qtlTypeNo]]
                               countNonQtl <- countNonQtlList[[qtlTypeNo]]
                               countBlocksForAuc <- countBlocksForAucList[[qtlTypeNo]]
                               aucLdBaseTP <- c(0, 1:countQtl, countQtl) / countQtl
                               aucLdBaseFP <- c(0, (sort(unique(countBlocksForAuc)) -
                                                      1:countQtl), countNonQtl) / countNonQtl
                               aucLdBaseFP[aucLdBaseFP <= 0] <- 0


                               aucLdBase <- aucLdBaseTP[2] * aucLdBaseFP[2] / 2
                               for (m in 2:(length(aucLdBaseFP) - 1)) {
                                 aucLdBaseSegment <- (aucLdBaseTP[m] + aucLdBaseTP[m + 1]) *
                                   (aucLdBaseFP[m + 1] - aucLdBaseFP[m]) / 2
                                 aucLdBase <- aucLdBase + aucLdBaseSegment
                               }


                               return(list(
                                 aucLdBaseTP = aucLdBaseTP,
                                 aucLdBaseFP = aucLdBaseFP,
                                 aucLdBase = aucLdBase
                               ))
                             }, simplify = FALSE)
      names(aucLdBaseRes) <- qtlTypeNamesWithTotal
      aucLdBase <- unlist(lapply(X = aucLdBaseRes,
                                 FUN = function(x) x$aucLdBase))

      if (plotROC) {
        for (qtlTypeNo in 1:(nQtlTypeWithTotal)) {
          qtlTypeName <- qtlTypeNamesWithTotal[qtlTypeNo]
          aucLdBaseTP <- aucLdBaseRes[[qtlTypeNo]]$aucLdBaseTP
          aucLdBaseFP <- aucLdBaseRes[[qtlTypeNo]]$aucLdBaseFP

          if (is.null(saveName)) {
            plot(aucLdBaseFP, aucLdBaseTP, type = "l",
                 main = qtlTypeName)
            text(0.83, 0.15, paste("AUC (LD) = ",
                                   round(aucLdBase[qtlTypeNo], 3), sep = ""),
                 cex = 1.5)
          } else{
            png(paste(saveName, traitName, "_", qtlTypeName,
                      "_LD_based_ROC_curve.png", sep = ""),
                width = 700)
            plot(aucLdBaseFP, aucLdBaseTP, type = "l",
                 main = qtlTypeName)
            text(0.83, 0.15, paste("AUC (LD) = ",
                                   round(aucLdBase[qtlTypeNo], 3), sep = ""),
                 cex = 1.5)
            dev.off()
          }
        }
      }
    } else {
      aucLdBase <- NULL
    }


    #### 1.7. Other summary statistics ####
    if (is.null(idsAroundMrcList)) {
      idsLdBlockListWithTotal <- c(Total = list(unlist(idsLdBlockList, recursive = FALSE)),
                                   idsLdBlockList)
    } else {
      idsLdBlockListWithTotal <- lapply(X = qtlCandsListWithTotal,
                                        FUN = function(qtlCandsListEach) {
                                          sapply(X = qtlCandsListEach,
                                                 FUN = function(qtlCandsEach) {
                                                   idsAroundMrcList[unlist(lapply(X = idsAroundMrcList,
                                                                                  FUN = function(idsAroundMrcEach) {
                                                                                    any(idsAroundMrcEach %in% qtlCandsEach)
                                                                                  }))][[1]]
                                                 }, simplify = FALSE)
                                        })
    }


    ### 1.7.1. Maximum -log10(p) values for each LD block around QTL ###
    maxLogpInLdListWithTotal <- lapply(X = idsLdBlockListWithTotal,
                                       FUN = function(idsLdBlockEach) {
                                         unlist(lapply(X = idsLdBlockEach,
                                                       FUN = function(idsLdBlock) {
                                                         maxLogpInLd <- max(resGWAS[idsLdBlock, 4])

                                                         return(maxLogpInLd)
                                                       }))
                                       })
    recallMaxLogpInLdListWithTotal <- lapply(
      X = maxLogpInLdListWithTotal,
      FUN = function(maxLogpInLd) {
        sapply(X = thresholds,
               FUN = function(threshold) {
                 if (!is.na(threshold)) {
                   recallMaxLogpInLd <- mean(maxLogpInLd >= threshold)
                 } else {
                   recallMaxLogpInLd <- 0
                 }

                 return(recallMaxLogpInLd)
               })
      }
    )


    ### 1.7.2. Inflation factor (mean of -log10(p) values of the top 10 false positive blocks) ###
    if (!returnBlockRes) {
      if (is.null(idsAroundMrcList)) {
        logpOrdForHm <- logpOrd
        qtlCandsTotal <- unlist(qtlCandsList)
        qtlCandsListWithTotal <- c(list(Total = qtlCandsTotal),
                                   qtlCandsList)
        countBlocks <- 0

        while (countBlocks <= nTrueQtlTotal + nTopBlockFalsePositive) {
          idMrcNow <- logpOrdForHm[1]
          idsAroundMrcRaw <- which(cumPos >= (cumPos[idMrcNow] - lenLD) &
                                     cumPos <= (cumPos[idMrcNow] + lenLD))
          corAroundMrcRaw <- (cor(genoScore[, idMrcNow],
                                  genoScore[, idsAroundMrcRaw, drop = FALSE]) ^ 2)[1, ]

          startBlockId <- max(min(which(x = corAroundMrcRaw >= corThres))
                              - windowSize, 1)
          endBlockId <- min(max(which(x = corAroundMrcRaw >= corThres))
                            + windowSize, length(idsAroundMrcRaw))
          idsAroundMrc <- idsAroundMrcRaw[startBlockId:endBlockId]
          idsAroundMrc <- idsAroundMrc[idsAroundMrc %in% logpOrdForHm]
          idsAroundMrcListNow <- list(idsAroundMrc)
          names(idsAroundMrcListNow) <- idMrcNow
          idsAroundMrcList <- c(idsAroundMrcList, idsAroundMrcListNow)



          logpOrdForHm <- logpOrdForHm[!(logpOrdForHm %in% idsAroundMrc)]

          countBlocks <- countBlocks + 1
        }
      }


      whichQtlCandsTotal <- unlist(parallel::mclapply(
        X = idsAroundMrcList,
        FUN = function(idsAroundMrc) {
          any(idsAroundMrc %in% qtlCandsListWithTotal$Total)
        },
        mc.cores = nCores
      ))

      falseBlockLogps <- unlist(x = lapply(X = (idsAroundMrcList[!whichQtlCandsTotal]),
                                           FUN = function(x) {
                                             suppressWarnings(max(resGWAS[x, 4], na.rm = TRUE))
                                           }))
      falseBlockLogpsList <- rep(list(falseBlockLogps), nQtlTypeWithTotal)
      names(falseBlockLogpsList) <- qtlTypeNamesWithTotal
    }


    meanFalseLogpListWithTotal <- lapply(
      X = falseBlockLogpsList,
      FUN = function(falseBlockLogps) {
        meanFalseLogp <- sapply(X = nTopBlockFalsePositive,
                                FUN = function(x) {
                                  mean(falseBlockLogps[1:x], na.rm = TRUE)
                                }, simplify = TRUE)
        names(meanFalseLogp) <- nTopBlockFalsePositive

        return(meanFalseLogp)
      }
    )



    ### 1.7.3. Adjusted -log10(p) values ###
    qtlLogpListWithTotal <- lapply(X = qtlCandsListWithTotal,
                                   FUN = function(qtlCands) {
                                     resGWAS[qtlCands, 4]
                                   })
    qtlAdjustedLogpListWithTotal <- mapply(FUN = "-",
                                           qtlLogpListWithTotal,
                                           meanFalseLogpListWithTotal)
    maxAdjustedLogpInLdListWithTotal <- mapply(FUN = "-",
                                               maxLogpInLdListWithTotal,
                                               meanFalseLogpListWithTotal)


    ### 1.7.4. AUC around causal ###
    if (returnAucRes) {
      qtlCandsListWithTotalList <- lapply(X = qtlCandsListWithTotal,
                                          FUN = function(x) split(x, 1:length(x)))
      aucCausalList <- mapply(
        FUN = function(qtlCandsListList,
                       idsLdBlock) {
          mapply(FUN = function(qtlCands,
                                idsLdBlockEach) {
            whereQtlCandsInBlock <- match(qtlCands, idsLdBlockEach)
            qtlLogpLdBlock <- resGWAS[idsLdBlockEach, 4]
            qtlLogpOrdLdBlock <- rank(-qtlLogpLdBlock)
            qtlLogpOrdInBlock <- qtlLogpOrdLdBlock[whereQtlCandsInBlock]

            aucCausalTP <- c(0, 1:length(qtlCands), length(qtlCands)) / length(qtlCands)
            if (length(idsLdBlockEach) - length(qtlCands) >= 1) {
              aucCausalFP <- (c(0, sort(qtlLogpOrdInBlock), length(idsLdBlockEach)) -
                                c(0, 1:length(qtlCands), length(qtlCands))) /
                (length(idsLdBlockEach) - length(qtlCands))
              aucCausal <- aucCausalTP[2] * aucCausalFP[2] / 2

              for (m in 2:(length(aucCausalFP) - 1)) {
                aucCausalSegment <- (aucCausalTP[m] + aucCausalTP[m + 1]) * (aucCausalFP[m + 1] - aucCausalFP[m]) / 2
                aucCausal <- aucCausal + aucCausalSegment
              }
            } else {
              aucCausal <- 1
            }
            return(aucCausal)
          },
          qtlCandsListList,
          idsLdBlock,
          SIMPLIFY = TRUE)
        },
        qtlCandsListWithTotalList,
        idsLdBlockListWithTotal,
        SIMPLIFY = FALSE)
      aucCausal <- sapply(X = aucCausalList, FUN = mean, simplify = TRUE)
    } else {
      aucCausal <- NULL
    }

    haploName <- NULL




    #### 1.8. Compute false discovery rate using chromosome under null model ####
    if (!is.null(nullChrNo)) {
      resGWASNullChr <- resGWAS[resGWAS[, 2] %in% nullChrNo, , drop = FALSE]

      fdrNullChrs <- sapply(
        X = thresholds,
        FUN = function(threshold) {
          if (!is.na(threshold)) {
            fdrNullChr <- sum(any(resGWASNullChr$Iteration_1 >= threshold))
          } else {
            fdrNullChr <- 0
          }

          return(fdrNullChr)
        }
      )

      propMarkersUnderSigLevels <- sapply(
        X = -log10(sigLevel),
        FUN = function(sigLevelEach) {
          mean(resGWASNullChr[, 4] >= sigLevelEach, na.rm = TRUE)
        }
      )

      names(fdrNullChrs) <- names(propMarkersUnderSigLevels) <- names(thresholds)

      nullChrRes <- list(nullChrNo = nullChrNo,
                         scoresNullChr = resGWASNullChr[, 4],
                         fdrNullChrs = fdrNullChrs,
                         propMarkersUnderSigLevels = propMarkersUnderSigLevels)
    } else {
      nullChrRes <- NULL
    }



    #### 1.9. Put all results into `summaryStatResList` ####
    qtlLogpOrdsListWithTotal <- c(Total = list(unlist(qtlLogpOrdsList)),
                                  qtlLogpOrdsList)
    for (qtlTypeName in qtlTypeNamesWithTotal) {
      summaryStatResList[[qtlTypeName]][["qtlLogp"]] <-
        qtlLogpListWithTotal[[qtlTypeName]]
      summaryStatResList[[qtlTypeName]][["qtlAdjustedLogp"]] <-
        qtlAdjustedLogpListWithTotal[[qtlTypeName]]
      summaryStatResList[[qtlTypeName]][["qtlLogpOrd"]] <-
        qtlLogpOrdsListWithTotal[[qtlTypeName]]
      summaryStatResList[[qtlTypeName]][["auc"]] <-
        auc[qtlTypeName]
      summaryStatResList[[qtlTypeName]][["aucLdBase"]] <-
        aucLdBase[qtlTypeName]
      summaryStatResList[[qtlTypeName]][["aucCausal"]] <-
        aucCausal[qtlTypeName]
      summaryStatResList[[qtlTypeName]][["meanFalseLogp"]] <-
        meanFalseLogpListWithTotal[[qtlTypeName]]
      summaryStatResList[[qtlTypeName]][["maxLogpInLd"]] <-
        maxLogpInLdListWithTotal[[qtlTypeName]]
      summaryStatResList[[qtlTypeName]][["recallMaxLogpInLd"]] <-
        recallMaxLogpInLdListWithTotal[[qtlTypeName]]
      summaryStatResList[[qtlTypeName]][["maxAdjustedLogpInLd"]] <-
        maxAdjustedLogpInLdListWithTotal[[qtlTypeName]]
      summaryStatResList[[qtlTypeName]][["falseBlockLogps"]] <-
        falseBlockLogpsList[[qtlTypeName]]
    }


  } else {
    ###### 2. For HB Results ######
    ##### 2.1. Chromosome and position ######
    resGWAS <- resGWAS[order(resGWAS[, 2], resGWAS[, 3]), ]
    block <- as.character(geneSet[, 1])
    blockNumber <- readr::parse_number(block)
    nHaploBlock <- nrow(resGWAS)

    cumPosHaplo <- RAINBOWR::cumsumPos(map = resGWAS[, 1:3])


    qtlNamesHaploList <- lapply(X = qtlNamesList,
                                FUN = function(qtlNames) {
                                  haploBlocksIncludeQTl <- sapply(
                                    X = qtlNames,
                                    FUN = function(mrkName) {
                                      defineHaploBlock(mrkName = mrkName,
                                                       haploBlockListDf = geneSet,
                                                       genoMap = genoMap,
                                                       resGWAS = resGWAS)
                                    }
                                  )
                                })

    qtlCandsHaploList <- lapply(X = qtlNamesHaploList,
                                FUN = function(qtlNamesHaplo) {
                                  match(x = qtlNamesHaplo,
                                        table = as.character(resGWAS[, 1]))
                                }
    )

    qtlCandsHaploTotal <- unlist(qtlCandsHaploList)
    qtlCandsHaploListWithTotal <- c(list(Total = qtlCandsHaploTotal),
                                    qtlCandsHaploList)




    ##### 2.2. Calculate theshold #####
    thresholds <- try(RAINBOWR::CalcThreshold(input = resGWAS,
                                              sig.level = sigLevel,
                                              method = methodThres), silent = TRUE)
    if ("try-error" %in% class(thresholds)) {
      thresholds <- rep(NA, length(sigLevel))
    }
    nOverThresholds <- sapply(X = thresholds,
                              FUN = function(threshold) {
                                ifelse(test = is.na(threshold),
                                       yes = 0,
                                       no = sum(resGWAS[, 4] >= threshold))
                              })


    ##### 2.3. The rank of -log10(p) of qtls #####
    logpOrd <- order(resGWAS[, 4], decreasing = TRUE)
    qtlHaploLogpOrdsList <- lapply(X = qtlCandsHaploList,
                                   FUN = function(qtlCandsHaplo) {
                                     match(qtlCandsHaplo, logpOrd)
                                   })

    positionQtlsHaploList <- lapply(X = qtlCandsHaploList,
                                    FUN = function(qtlCandsHaplo) {
                                      cumPosHaplo[qtlCandsHaplo]
                                    })


    ldBlockResList <- sapply(
      X = 1:nQtlType,
      FUN = function(qtlTypeNo) {
        qtlCands <- qtlCandsList[[qtlTypeNo]]
        qtlCandsHaplo <- qtlCandsHaploList[[qtlTypeNo]]
        positionQtls <- positionQtlsList[[qtlTypeNo]]
        positionQtlsHaplo <- positionQtlsHaploList[[qtlTypeNo]]
        nTrueQtl <- nTrueQtls[qtlTypeNo]

        ldBlockRes <- sapply(X = 1:nTrueQtl,
                             FUN = function(trueQtlNo) {
                               qtlCand <- qtlCands[trueQtlNo]
                               qtlCandHaplo <- qtlCandsHaplo[trueQtlNo]
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

                               idsLdBlockHaploStEd <- readr::parse_number(unique(sapply(
                                 X = marker[idsLdBlock[c(1, length(idsLdBlock))]],
                                 FUN = defineHaploBlock,
                                 haploBlockListDf = geneSet,
                                 genoMap = genoMap
                               )))
                               if (length(idsLdBlockHaploStEd) >= 2) {
                                 idsLdBlockHaplo <- idsLdBlockHaploStEd[1]:idsLdBlockHaploStEd[2]
                               } else {
                                 idsLdBlockHaplo <- idsLdBlockHaploStEd
                               }

                               idsLdBlockHaploLogpOrd <- idsLdBlockHaplo[order(resGWAS[idsLdBlockHaplo, 4], decreasing = TRUE)]
                               overQtlHaploInLdBlock <- idsLdBlockHaploLogpOrd[
                                 1:match(qtlCandHaplo, idsLdBlockHaploLogpOrd)
                               ]

                               return(list(idsLdBlock = idsLdBlockHaplo,
                                           overQtlInLdBlock = overQtlHaploInLdBlock))

                             }, simplify = FALSE)

        return(ldBlockRes)
      },
      simplify = FALSE
    )
    names(ldBlockResList) <- qtlTypeNames
    idsLdBlockList <- lapply(X = ldBlockResList,
                             FUN = function(ldBlockResults) {
                               lapply(X = ldBlockResults,
                                      FUN = function(ldBlockRes) {
                                        ldBlockRes$idsLdBlock
                                      })
                             })

    overQtlInLdBlockList <- lapply(X = ldBlockResList,
                                   FUN = function(ldBlockResults) {
                                     lapply(X = ldBlockResults,
                                            FUN = function(ldBlockRes) {
                                              ldBlockRes$overQtlInLdBlock
                                            })
                                   })
    qtlHaploLogpOrdsVec <- unlist(qtlHaploLogpOrdsList)



    #### 2.4. Calculation of AUC and draw ROC curve ####
    if (returnAucRes) {
      overQtlInLdBlockListUnlist <- unlist(overQtlInLdBlockList, recursive = FALSE)
      overQtlInLdBlockListUnlistOrd <- overQtlInLdBlockListUnlist[order(qtlHaploLogpOrdsVec)]
      nOverQtlOrdInLdBlocksCumsumTotal <- rep(NA, nTrueQtlTotal)
      overQtlOrdInLdBlocksVec <- NULL
      for (trueQtlNoTotal in 1:nTrueQtlTotal) {
        overQtlOrdInLdBlocksVec <- c(overQtlOrdInLdBlocksVec,
                                     overQtlInLdBlockListUnlistOrd[[trueQtlNoTotal]])
        nOverQtlOrdInLdBlocksCumsumTotal[order(qtlHaploLogpOrdsVec)[trueQtlNoTotal]] <- length(unique(overQtlOrdInLdBlocksVec))
      }
      nOverQtlOrdInLdBlocksCumsumTotalSorted <- sort(nOverQtlOrdInLdBlocksCumsumTotal)

      nOverQtlOrdInLdBlocksCumsumList <- split(
        x = nOverQtlOrdInLdBlocksCumsumTotal,
        f = rep(qtlTypeNames, nTrueQtls)
      )[qtlTypeNames]
      nOverQtlOrdInLdBlocksCumsumSotedList <- lapply(
        X = nOverQtlOrdInLdBlocksCumsumList, FUN = sort
      )




      #### 2.4.1. Each result of AUC and draw ROC curve ####
      aucEachRes <- sapply(
        X = 1:nQtlType,
        FUN = function(qtlTypeNo) {
          qtlHaploLogpOrds <- qtlHaploLogpOrdsList[[qtlTypeNo]]
          nTrueQtl <- nTrueQtls[qtlTypeNo]
          nOverQtlOrdInLdBlocksCumsum <-
            nOverQtlOrdInLdBlocksCumsumSotedList[[qtlTypeNo]]

          aucTP <- c(0, 1:nTrueQtl, nTrueQtl) / nTrueQtl
          aucFP <- (c(0, sort(qtlHaploLogpOrds), nMarkers) -
                      c(0, nOverQtlOrdInLdBlocksCumsum,
                        nOverQtlOrdInLdBlocksCumsum[nTrueQtl])) /
            (nMarkers - nOverQtlOrdInLdBlocksCumsum[nTrueQtl])

          aucFP[aucFP <= 0] <- 0


          auc <- aucTP[2] * aucFP[2] / 2

          for (m in 2:(length(aucFP) - 1)) {
            aucSegment <-
              (aucTP[m] + aucTP[m + 1]) * (aucFP[m + 1] - aucFP[m]) / 2
            auc <- auc + aucSegment
          }

          return(list(aucTP = aucTP,
                      aucFP = aucFP,
                      auc = auc))
        }, simplify = FALSE
      )
      names(aucEachRes) <- qtlTypeNames
      aucEach <- sapply(X = aucEachRes,
                        FUN = function(x) x$auc,
                        simplify = TRUE)


      if (plotROC) {
        for (qtlTypeNo in 1:nQtlType) {
          qtlTypeName <- qtlTypeNames[qtlTypeNo]
          aucEachNow <- aucEach[qtlTypeNo]
          aucTP <- aucEachRes[[qtlTypeNo]]$aucTP
          aucFP <- aucEachRes[[qtlTypeNo]]$aucFP

          if (is.null(saveName)) {
            plot(aucFP, aucTP, type = "l",
                 main = qtlTypeName)
            text(0.83, 0.15, paste("AUC = ", round(aucEachNow, 3), sep = ""), cex = 2.5)
          } else{
            png(paste(saveName, traitName, "_", qtlTypeName,
                      "_ROC_curve.png", sep = ""),
                width = 700)
            plot(aucFP, aucTP, type = "l",
                 main = qtlTypeName)
            text(0.83, 0.15, paste("AUC = ", round(aucEachNow, 3), sep = ""), cex = 2.5)
            dev.off()
          }
        }
      }



      #### 2.4.2. Total result of AUC and draw ROC curve ####
      aucTPTotal <- c(0, 1:nTrueQtlTotal, nTrueQtlTotal) / nTrueQtlTotal
      aucFPTotal <- (c(0, sort(qtlHaploLogpOrdsVec), nHaploBlock) -
                       c(0, nOverQtlOrdInLdBlocksCumsumTotalSorted,
                         nOverQtlOrdInLdBlocksCumsumTotalSorted[nTrueQtlTotal])) /
        (nHaploBlock - nOverQtlOrdInLdBlocksCumsumTotalSorted[nTrueQtlTotal])

      aucFPTotal[aucFPTotal <= 0] <- 0

      aucTotal <- aucTPTotal[2] * aucFPTotal[2] / 2
      for (m in 2:(length(aucFPTotal) - 1)) {
        aucSegmentTotal <-
          (aucTPTotal[m] + aucTPTotal[m + 1]) *
          (aucFPTotal[m + 1] - aucFPTotal[m]) / 2
        aucTotal <- aucTotal + aucSegmentTotal
      }

      if (plotROC) {
        if (is.null(saveName)) {
          plot(aucFPTotal, aucTPTotal, type = "l",
               main = "Total")
          text(0.83, 0.15, paste("AUC = ", round(aucTotal, 3), sep = ""), cex = 2.5)
        } else{
          png(paste(saveName, traitName, "_Total",
                    "_ROC_curve.png", sep = ""),
              width = 700)
          plot(aucFPTotal, aucTPTotal, type = "l",
               main = "Total")
          text(0.83, 0.15, paste("AUC = ", round(aucTotal, 3), sep = ""), cex = 2.5)
          dev.off()
        }
      }
      names(aucTotal) <- "Total"

      auc <- c(aucTotal, aucEach)
    } else {
      auc <- NULL
    }




    #### 2.5. Calculate FPR, FNR, FDR and Harmonic mean of precision and recall ####
    ### 2.5.1. Markers over thresholds ###
    if (returnBlockRes) {
      overThresholdVals <- sapply(X = thresholds,
                                  FUN = function(threshold) {
                                    resGWAS[which(resGWAS[, 4] >= threshold), 4]
                                  }, simplify = FALSE)

      overThresholdCandsList <- sapply(X = thresholds,
                                       FUN = function(threshold) {
                                         which(resGWAS[, 4] >= threshold)
                                       }, simplify = FALSE)

      # overThresholdCandsLogpOrdList <-
      #   sapply(X = nOverThresholds,
      #          FUN = function(nOverThreshold) {
      #            if (nOverThreshold != 0) {
      #              overThresholdCandsLogpOrd <- logpOrd[1:nOverThreshold]
      #            } else {
      #              overThresholdCandsLogpOrd <- NULL
      #            }
      #            return(overThresholdCandsLogpOrd)
      #          }, simplify = FALSE)



      ### 2.5.2. Preparation for computation of number of TP, FP, FN, and TN ###
      # overThresholdCandsLogpOrdListForHm <- overThresholdCandsLogpOrdList
      logpOrdForHm <- logpOrd




      ### 2.5.3. Start computation of number of TP, FP, FN, and TN ###
      if (is.null(idsAroundMrcList)) {
        while (length(logpOrdForHm) >= 1) {
          ## 2.5.3.1. Detect LD block for the marker of interest ##
          idMrcHaploNow <- logpOrdForHm[1]
          markersHaploNow <- geneSet[blockNumber %in% idMrcHaploNow, 2]
          idMrcNow <- which(marker %in% markersHaploNow)
          idsAroundMrcRaw <- which(cumPos >= (cumPos[min(idMrcNow)] - lenLD) &
                                     cumPos <= (cumPos[max(idMrcNow)] + lenLD))
          corAroundMrcRaw <- (cor(genoScore[, idMrcNow],
                                  genoScore[, idsAroundMrcRaw, drop = FALSE]) ^ 2)


          startBlockId <- max(min(apply(X = corAroundMrcRaw,
                                        MARGIN = 1,
                                        FUN = function(x) {
                                          min(which(x = x >= corThres))
                                        }))
                              - windowSize, 1)
          endBlockId <- min(max(apply(X = corAroundMrcRaw,
                                      MARGIN = 1,
                                      FUN = function(x) {
                                        max(which(x = x >= corThres))
                                      }))
                            + windowSize, length(idsAroundMrcRaw))
          idsAroundMrc <- idsAroundMrcRaw[startBlockId:endBlockId]

          idsAroundMrcHaploStEd <- readr::parse_number(unique(sapply(
            X = marker[idsAroundMrc[c(1, length(idsAroundMrc))]],
            FUN = defineHaploBlock,
            haploBlockListDf = geneSet,
            genoMap = genoMap
          )))
          if (length(idsAroundMrcHaploStEd) >= 2) {
            idsAroundMrcHaplo <- idsAroundMrcHaploStEd[1]:idsAroundMrcHaploStEd[2]
          } else {
            idsAroundMrcHaplo <- idsAroundMrcHaploStEd
          }
          idsAroundMrcHaplo <- idsAroundMrcHaplo[idsAroundMrcHaplo %in% logpOrdForHm]

          idsAroundMrcListNow <- list(idsAroundMrcHaplo)
          names(idsAroundMrcListNow) <- unique(block)[idMrcHaploNow]
          idsAroundMrcList <- c(idsAroundMrcList, idsAroundMrcListNow)



          ## 2.5.3.1. Detect LD block for the marker of interest ##
          logpOrdForHm <- logpOrdForHm[!(logpOrdForHm %in% idsAroundMrcHaplo)]

        }
      } else {
        blockLogPs <- unlist(x = lapply(X = (idsAroundMrcList),
                                        FUN = function(x) {
                                          suppressWarnings(max(resGWAS[x, 4], na.rm = TRUE))
                                        }))
        idsAroundMrcList <- idsAroundMrcList[order(blockLogPs, decreasing = TRUE)]
      }


      whichQtlCandsList <- lapply(
        X = qtlCandsHaploList,
        FUN = function(qtlCandsHaplo) {
          unlist(parallel::mclapply(
            X = idsAroundMrcList,
            FUN = function(idsAroundMrc) {
              any(idsAroundMrc %in% qtlCandsHaplo)
            },
            mc.cores = nCores
          ))
        }
      )


      whichQtlCandsListWithTotal <- c(
        list(Total = apply(X = do.call(what = cbind,
                                       args = whichQtlCandsList),
                           MARGIN = 1, FUN = any)),
        whichQtlCandsList
      )

      # whichQtlCandsTotalMat <- do.call(what = rbind,
      #                                  args = parallel::mclapply(
      #   X = idsAroundMrcList,
      #   FUN = function(idsAroundMrc) {
      #     sapply(X = qtlCandsListWithTotal$Total,
      #            FUN = function(qtlCands) {
      #              any(idsAroundMrc %in% qtlCands)
      #            })
      #   },
      #   mc.cores = nCores
      # ))
      #
      # whichNoQtlCandsTotal <- apply(X = whichQtlCandsTotalMat,
      #                             MARGIN = 2, FUN = which)
      # names(whichNoQtlCandsTotal) <- rownames(whichQtlCandsTotalMat)[whichNoQtlCandsTotal]
      # whichNoQtlCandsList <- c(list(Total = whichNoQtlCandsTotal),
      #                        split(x = whichNoQtlCandsTotal,
      #                              f = rep(qtlTypeNames, nTrueQtls)))
      # names(whichNoQtlCandsList) <- qtlTypeNamesWithTotal



      countQtlList <- lapply(X = whichQtlCandsListWithTotal,
                             FUN = sum)
      falseBlockLogps <- unlist(x = lapply(X = (idsAroundMrcList[!whichQtlCandsListWithTotal$Total]),
                                           FUN = function(x) {
                                             suppressWarnings(max(resGWAS[x, 4], na.rm = TRUE))
                                           }))
      falseBlockLogpsList <- rep(list(falseBlockLogps), nQtlTypeWithTotal)
      countFP <- sapply(X = thresholds,
                        FUN = function(threshold) {
                          if (!is.na(threshold)) {
                            return(sum(falseBlockLogps >= threshold))
                          } else {
                            return(0)
                          }
                        }, simplify = TRUE)
      countFPList <- rep(list(countFP), nQtlTypeWithTotal)

      names(countFPList) <- names(falseBlockLogpsList) <-
        qtlTypeNamesWithTotal

      countAllBlocksList <- lapply(X = countQtlList,
                                   FUN = function(countQtl) {
                                     length(idsAroundMrcList) -
                                       (countQtlList$Total - countQtl)
                                   })

      countNonQtlList <- lapply(X = countQtlList,
                                FUN = function(countQtl) {
                                  length(idsAroundMrcList) - countQtlList$Total
                                })

      # qtlCandsOverThresholdList <- lapply(
      #   X = qtlCandsHaploListWithTotal,
      #   FUN = function(qtlCandsHaplo) {
      #     sapply(X = thresholds,
      #            FUN = function(threshold) {
      #              if (!is.na(threshold)) {
      #                qtlCandsHaplo[resGWAS[qtlCandsHaplo, 4] >= threshold]
      #              } else {
      #                NULL
      #              }
      #            }, simplify = FALSE)
      #   }
      # )
      #
      # countTPList <- lapply(
      #   X = qtlCandsOverThresholdList,
      #   FUN = function(qtlCandsOverThreshold) {
      #     unlist(lapply(X = qtlCandsOverThreshold,
      #                   FUN = function(qtlCandsOverThresholdNow) {
      #                     if (length(qtlCandsOverThresholdNow) >= 1) {
      #                       return(sum(unlist(parallel::mclapply(
      #                         X = idsAroundMrcList,
      #                         FUN = function(idsAroundMrc) {
      #                           any(idsAroundMrc %in% qtlCandsOverThresholdNow)
      #                         },
      #                         mc.cores = nCores
      #                       ))))
      #                     } else {
      #                       return(0)
      #                     }
      #                   }))
      #   }
      # )


      trueBlockLogpsList <- lapply(X = whichQtlCandsListWithTotal,
                                   FUN = function(whichQtlCandsEach) {
                                     unlist(x = lapply(X = (idsAroundMrcList[whichQtlCandsEach]),
                                                       FUN = function(x) {
                                                         suppressWarnings(max(resGWAS[x, 4], na.rm = TRUE))
                                                       }))
                                   })
      countTPList <- lapply(X = trueBlockLogpsList,
                            FUN = function(trueBlockLogps) {
                              sapply(X = thresholds,
                                     FUN = function(threshold) {
                                       if (!is.na(threshold)) {
                                         return(sum(trueBlockLogps >= threshold))
                                       } else {
                                         return(0)
                                       }
                                     }, simplify = TRUE)
                            })

      whichNoQtlCandsListWithTotal <- lapply(X = whichQtlCandsListWithTotal,
                                             FUN = which)

      countBlocksForAucList <- lapply(
        X = whichNoQtlCandsListWithTotal,
        FUN = function(whichNoQtlCands) {
          whichNoQtlCandsNames <- readr::parse_number(names(whichNoQtlCands))
          whichNoQtlCandsTotalNames <- readr::parse_number(names(whichNoQtlCandsListWithTotal$Total))
          countBlocksForAucUnSorted <- whichNoQtlCands -
            match(whichNoQtlCandsNames, whichNoQtlCandsTotalNames) +
            1:length(whichNoQtlCands)

          countBlocksForAuc <- countBlocksForAucUnSorted[rank(whichNoQtlCandsNames)]
        }
      )
      countBlocksForAucList$Total <- countBlocksForAucList$Total[
        unlist(sapply(X = countBlocksForAucList[-1],
                      FUN = names, simplify = FALSE))
      ]




      nTPList <- countTPList
      nFPList <- countFPList
      nFNList <- mapply(FUN = "-",
                        countQtlList,
                        countTPList,
                        SIMPLIFY = FALSE)
      nTNList <- mapply(
        FUN = function(countAllBlocks,
                       nTP, nFP, nFN) {
          countAllBlocks - nTP - nFP - nFN
        },
        countAllBlocks = countAllBlocksList,
        nTP = nTPList,
        nFP = nFPList,
        nFN = nFNList,
        SIMPLIFY = FALSE)
      names(nTNList) <- names(nFNList)

      summaryStatResList <- mapply(FUN = function(nTP, nFP, nFN, nTN) {
        fdr <- nFP / (nFP + nTP)
        fpr <- nFP / (nFP + nTN)
        fnr <- nFN / (nTP + nFN)
        recall <- nTP / (nTP + nFN)
        precision <- nTP / (nTP + nFP)
        accuracy <- (nTP + nTN) / (nTP + nFP + nTN + nFN)
        harmonicMean <- 2 * recall * precision / (recall + precision)

        if (any(is.nan(precision))) {
          fdr[is.nan(precision)] <- NA
          harmonicMean[is.nan(precision)] <- NA
          precision[is.nan(precision)] <- NA
        }

        if (any(recall == 0 | precision == 0)) {
          harmonicMean[recall == 0 | precision == 0] <- 0
        }

        return(list(
          fdr = fdr,
          fpr = fpr,
          fnr = fnr,
          recall = recall,
          precision = precision,
          accuracy = accuracy,
          harmonicMean = harmonicMean
        ))
      },
      nTP = nTPList,
      nFP = nFPList,
      nFN = nFNList,
      nTN = nTNList,
      SIMPLIFY = FALSE)
    } else {
      summaryStatResList <- rep(list(list(
        fdr = NULL,
        fpr = NULL,
        fnr = NULL,
        recall = NULL,
        precision = NULL,
        accuracy = NULL,
        harmonicMean = NULL
      )), nQtlTypeWithTotal)
      names(summaryStatResList) <- qtlTypeNamesWithTotal
    }







    #### 2.6. LD-based AUC ####
    if (returnAucRes) {
      aucLdBaseRes <- sapply(X = 1:(nQtlTypeWithTotal),
                             FUN = function(qtlTypeNo) {
                               countQtl <- countQtlList[[qtlTypeNo]]
                               countNonQtl <- countNonQtlList[[qtlTypeNo]]
                               countBlocksForAuc <- countBlocksForAucList[[qtlTypeNo]]
                               aucLdBaseTP <- c(0, 1:countQtl, countQtl) / countQtl
                               aucLdBaseFP <- c(0, (sort(unique(countBlocksForAuc)) -
                                                      1:countQtl), countNonQtl) / countNonQtl
                               aucLdBaseFP[aucLdBaseFP <= 0] <- 0


                               aucLdBase <- aucLdBaseTP[2] * aucLdBaseFP[2] / 2
                               for (m in 2:(length(aucLdBaseFP) - 1)) {
                                 aucLdBaseSegment <- (aucLdBaseTP[m] + aucLdBaseTP[m + 1]) *
                                   (aucLdBaseFP[m + 1] - aucLdBaseFP[m]) / 2
                                 aucLdBase <- aucLdBase + aucLdBaseSegment
                               }


                               return(list(
                                 aucLdBaseTP = aucLdBaseTP,
                                 aucLdBaseFP = aucLdBaseFP,
                                 aucLdBase = aucLdBase
                               ))
                             }, simplify = FALSE)
      names(aucLdBaseRes) <- qtlTypeNamesWithTotal
      aucLdBase <- unlist(lapply(X = aucLdBaseRes,
                                 FUN = function(x) x$aucLdBase))

      if (plotROC) {
        for (qtlTypeNo in 1:(nQtlTypeWithTotal)) {
          qtlTypeName <- qtlTypeNamesWithTotal[qtlTypeNo]
          aucLdBaseTP <- aucLdBaseRes[[qtlTypeNo]]$aucLdBaseTP
          aucLdBaseFP <- aucLdBaseRes[[qtlTypeNo]]$aucLdBaseFP

          if (is.null(saveName)) {
            plot(aucLdBaseFP, aucLdBaseTP, type = "l",
                 main = qtlTypeName)
            text(0.83, 0.15, paste("AUC (LD) = ",
                                   round(aucLdBase[qtlTypeNo], 3), sep = ""),
                 cex = 2.5)
          } else{
            png(paste(saveName, traitName, "_", qtlTypeName,
                      "_LD_based_ROC_curve.png", sep = ""),
                width = 700)
            plot(aucLdBaseFP, aucLdBaseTP, type = "l",
                 main = qtlTypeName)
            text(0.83, 0.15, paste("AUC (LD) = ",
                                   round(aucLdBase[qtlTypeNo], 3), sep = ""),
                 cex = 2.5)
            dev.off()
          }
        }
      }
    } else {
      aucLdBase <- NULL
    }


    #### 2.7. Other summary statistics ####
    if (is.null(idsAroundMrcList)) {
      idsLdBlockListWithTotal <- c(Total = list(unlist(idsLdBlockList, recursive = FALSE)),
                                   idsLdBlockList)
    } else {
      idsLdBlockListWithTotal <- lapply(X = qtlCandsHaploListWithTotal,
                                        FUN = function(qtlCandsHaploListEach) {
                                          sapply(X = qtlCandsHaploListEach,
                                                 FUN = function(qtlCandsHaploEach) {
                                                   idsAroundMrcList[unlist(lapply(X = idsAroundMrcList,
                                                                                  FUN = function(idsAroundMrcEach) {
                                                                                    any(idsAroundMrcEach %in% qtlCandsHaploEach)
                                                                                  }))][[1]]
                                                 }, simplify = FALSE)
                                        })
    }

    ### 2.7.1. Maximum -log10(p) values for each LD block around QTL ###
    maxLogpInLdListWithTotal <- lapply(X = idsLdBlockListWithTotal,
                                       FUN = function(idsLdBlockEach) {
                                         unlist(lapply(X = idsLdBlockEach,
                                                       FUN = function(idsLdBlock) {
                                                         maxLogpInLd <- max(resGWAS[idsLdBlock, 4], na.rm = TRUE)

                                                         return(maxLogpInLd)
                                                       }))
                                       })
    recallMaxLogpInLdListWithTotal <- lapply(
      X = maxLogpInLdListWithTotal,
      FUN = function(maxLogpInLd) {
        sapply(X = thresholds,
               FUN = function(threshold) {
                 if (!is.na(threshold)) {
                   recallMaxLogpInLd <- mean(maxLogpInLd >= threshold)
                 } else {
                   recallMaxLogpInLd <- 0
                 }

                 return(recallMaxLogpInLd)
               })
      }
    )


    ### 2.7.2. Inflation factor (mean of -log10(p) values of the top 10 false positive blocks) ###
    if (!returnBlockRes) {
      if (is.null(idsAroundMrcList)) {
        logpOrdForHm <- logpOrd
        qtlCandsTotal <- unlist(qtlCandsList)
        qtlCandsListWithTotal <- c(list(Total = qtlCandsTotal),
                                   qtlCandsList)
        countBlocks <- 0

        while (countBlocks <= nTrueQtlTotal + nTopBlockFalsePositive) {
          idMrcHaploNow <- logpOrdForHm[1]
          markersHaploNow <- geneSet[blockNumber %in% idMrcHaploNow, 2]
          idMrcNow <- which(marker %in% markersHaploNow)
          idsAroundMrcRaw <- which(cumPos >= (cumPos[min(idMrcNow)] - lenLD) &
                                     cumPos <= (cumPos[max(idMrcNow)] + lenLD))
          corAroundMrcRaw <- (cor(genoScore[, idMrcNow],
                                  genoScore[, idsAroundMrcRaw, drop = FALSE]) ^ 2)


          startBlockId <- max(min(apply(X = corAroundMrcRaw,
                                        MARGIN = 1,
                                        FUN = function(x) {
                                          min(which(x = x >= corThres))
                                        }))
                              - windowSize, 1)
          endBlockId <- min(max(apply(X = corAroundMrcRaw,
                                      MARGIN = 1,
                                      FUN = function(x) {
                                        max(which(x = x >= corThres))
                                      }))
                            + windowSize, length(idsAroundMrcRaw))
          idsAroundMrc <- idsAroundMrcRaw[startBlockId:endBlockId]

          idsAroundMrcHaploStEd <- readr::parse_number(unique(sapply(
            X = marker[idsAroundMrc[c(1, length(idsAroundMrc))]],
            FUN = defineHaploBlock,
            haploBlockListDf = geneSet,
            genoMap = genoMap
          )))
          if (length(idsAroundMrcHaploStEd) >= 2) {
            idsAroundMrcHaplo <- idsAroundMrcHaploStEd[1]:idsAroundMrcHaploStEd[2]
          } else {
            idsAroundMrcHaplo <- idsAroundMrcHaploStEd
          }
          idsAroundMrcHaplo <- idsAroundMrcHaplo[idsAroundMrcHaplo %in% logpOrdForHm]

          idsAroundMrcListNow <- list(idsAroundMrcHaplo)
          names(idsAroundMrcListNow) <- unique(block)[idMrcHaploNow]
          idsAroundMrcList <- c(idsAroundMrcList, idsAroundMrcListNow)

          logpOrdForHm <- logpOrdForHm[!(logpOrdForHm %in% idsAroundMrcHaplo)]

          countBlocks <- countBlocks + 1
        }
      }


      whichQtlCandsTotal <- unlist(parallel::mclapply(
        X = idsAroundMrcList,
        FUN = function(idsAroundMrc) {
          any(idsAroundMrc %in% qtlCandsHaploListWithTotal$Total)
        },
        mc.cores = nCores
      ))

      falseBlockLogps <- unlist(x = lapply(X = (idsAroundMrcList[!whichQtlCandsTotal]),
                                           FUN = function(x) {
                                             suppressWarnings(max(resGWAS[x, 4], na.rm = TRUE))
                                           }))
      falseBlockLogpsList <- rep(list(falseBlockLogps), nQtlTypeWithTotal)
      names(falseBlockLogpsList) <- qtlTypeNamesWithTotal
    }


    meanFalseLogpListWithTotal <- lapply(
      X = falseBlockLogpsList,
      FUN = function(falseBlockLogps) {
        meanFalseLogp <- sapply(X = nTopBlockFalsePositive,
                                FUN = function(x) {
                                  mean(falseBlockLogps[1:x], na.rm = TRUE)
                                }, simplify = TRUE)
        names(meanFalseLogp) <- nTopBlockFalsePositive

        return(meanFalseLogp)
      }
    )



    ### 2.7.3. Adjusted -log10(p) values ###
    qtlHaploLogpListWithTotal <- lapply(X = qtlCandsHaploListWithTotal,
                                        FUN = function(qtlCandsHaplo) {
                                          resGWAS[qtlCandsHaplo, 4]
                                        })
    qtlAdjustedLogpListWithTotal <- mapply(FUN = "-",
                                           qtlHaploLogpListWithTotal,
                                           meanFalseLogpListWithTotal)
    maxAdjustedLogpInLdListWithTotal <- mapply(FUN = "-",
                                               maxLogpInLdListWithTotal,
                                               meanFalseLogpListWithTotal)


    ### 2.7.4. AUC around causal ###
    if (returnAucRes) {
      qtlCandsHaploListWithTotalList <- lapply(X = qtlCandsHaploListWithTotal,
                                               FUN = function(x) split(x, 1:length(x)))
      aucCausalList <- mapply(
        FUN = function(qtlCandsHaploListList,
                       idsLdBlock) {
          mapply(FUN = function(qtlCandsHaplo,
                                idsLdBlockEach) {
            whereQtlCandsInBlock <- match(qtlCandsHaplo, idsLdBlockEach)
            qtlLogpLdBlock <- resGWAS[idsLdBlockEach, 4]
            qtlLogpOrdLdBlock <- rank(-qtlLogpLdBlock)
            qtlLogpOrdInBlock <- qtlLogpOrdLdBlock[whereQtlCandsInBlock]

            aucCausalTP <- c(0, 1:length(qtlCandsHaplo), length(qtlCandsHaplo)) / length(qtlCandsHaplo)
            if (length(idsLdBlockEach) - length(qtlCandsHaplo) >= 1) {
              aucCausalFP <- (c(0, sort(qtlLogpOrdInBlock), length(idsLdBlockEach)) -
                                c(0, 1:length(qtlCandsHaplo), length(qtlCandsHaplo))) /
                (length(idsLdBlockEach) - length(qtlCandsHaplo))
              aucCausal <- aucCausalTP[2] * aucCausalFP[2] / 2

              for (m in 2:(length(aucCausalFP) - 1)) {
                aucCausalSegment <- (aucCausalTP[m] + aucCausalTP[m + 1]) * (aucCausalFP[m + 1] - aucCausalFP[m]) / 2
                aucCausal <- aucCausal + aucCausalSegment
              }
            } else {
              aucCausal <- 1
            }
            return(aucCausal)
          },
          qtlCandsHaploListList,
          idsLdBlock,
          SIMPLIFY = TRUE)
        },
        qtlCandsHaploListWithTotalList,
        idsLdBlockListWithTotal,
        SIMPLIFY = FALSE)
      aucCausal <- sapply(X = aucCausalList, FUN = mean, simplify = TRUE)
    } else {
      aucCausal <- NULL
    }

    haploName <- NULL




    #### 2.8. Compute false discovery rate using chromosome under null model ####
    if (!is.null(nullChrNo)) {
      resGWASNullChr <- resGWAS[resGWAS[, 2] %in% nullChrNo, , drop = FALSE]

      fdrNullChrs <- sapply(
        X = thresholds,
        FUN = function(threshold) {
          if (!is.na(threshold)) {
            fdrNullChr <- sum(any(resGWASNullChr$Iteration_1 >= threshold))
          } else {
            fdrNullChr <- 0
          }

          return(fdrNullChr)
        }
      )

      propMarkersUnderSigLevels <- sapply(
        X = -log10(sigLevel),
        FUN = function(sigLevelEach) {
          mean(resGWASNullChr[, 4] >= sigLevelEach, na.rm = TRUE)
        }
      )

      names(fdrNullChrs) <- names(propMarkersUnderSigLevels) <- names(thresholds)

      nullChrRes <- list(nullChrNo = nullChrNo,
                         scoresNullChr = resGWASNullChr[, 4],
                         fdrNullChrs = fdrNullChrs,
                         propMarkersUnderSigLevels = propMarkersUnderSigLevels)
    } else {
      nullChrRes <- NULL
    }



    #### 2.9. Put all results into `summaryStatResList` ####
    qtlHaploLogpOrdsListWithTotal <- c(Total = list(unlist(qtlHaploLogpOrdsList)),
                                       qtlHaploLogpOrdsList)
    for (qtlTypeName in qtlTypeNamesWithTotal) {
      summaryStatResList[[qtlTypeName]][["qtlLogp"]] <-
        qtlHaploLogpListWithTotal[[qtlTypeName]]
      summaryStatResList[[qtlTypeName]][["qtlAdjustedLogp"]] <-
        qtlAdjustedLogpListWithTotal[[qtlTypeName]]
      summaryStatResList[[qtlTypeName]][["qtlLogpOrd"]] <-
        qtlHaploLogpOrdsListWithTotal[[qtlTypeName]]
      summaryStatResList[[qtlTypeName]][["auc"]] <-
        auc[qtlTypeName]
      summaryStatResList[[qtlTypeName]][["aucLdBase"]] <-
        aucLdBase[qtlTypeName]
      summaryStatResList[[qtlTypeName]][["aucCausal"]] <-
        aucCausal[qtlTypeName]
      summaryStatResList[[qtlTypeName]][["meanFalseLogp"]] <-
        meanFalseLogpListWithTotal[[qtlTypeName]]
      summaryStatResList[[qtlTypeName]][["maxLogpInLd"]] <-
        maxLogpInLdListWithTotal[[qtlTypeName]]
      summaryStatResList[[qtlTypeName]][["recallMaxLogpInLd"]] <-
        recallMaxLogpInLdListWithTotal[[qtlTypeName]]
      summaryStatResList[[qtlTypeName]][["maxAdjustedLogpInLd"]] <-
        maxAdjustedLogpInLdListWithTotal[[qtlTypeName]]
      summaryStatResList[[qtlTypeName]][["falseBlockLogps"]] <-
        falseBlockLogpsList[[qtlTypeName]]
    }


  }


  return(
    list(
      summaryStatResList = summaryStatResList,
      nullChrRes = nullChrRes,
      thresholds = thresholds,
      nOverThresholds = nOverThresholds,
      idsAroundMrcList = idsAroundMrcList,
      haploName = haploName
    )
  )
}
