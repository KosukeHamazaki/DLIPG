############################################################################################
######  Title: 1.1_Soybean_DLIPG_function_to_perform_GWAS_by_Rio_model_with_MM4LMM    ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                              ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo     ######
######  Date: 2022/06/26 (Created), 2022/07/13 (Last Updated)                         ######
############################################################################################




performGWASRio <- function(genoScore,
                           genoMap,
                           phenoNow,
                           addGRM,
                           gbInfo,
                           gbNow,
                           genoScoreWithGenBack,
                           nCores = 1,
                           mafThres = 0.05,
                           ldThres = 0.95,
                           adjustAddGRM = FALSE) {
  if (-1 %in% genoScore) {
    genoScore <- genoScore + 1
  }
  phenoNow <- as.matrix(phenoNow)
  nLines <- nrow(genoScore)
  nSNPs <- ncol(genoScore)

  if (length(gbNow) == 1) {
    ### M1 model
    ## Dataset for subgroup of interest
    genoScoreSubGrRaw <- genoScore[gbInfo %in% gbNow, ]
    mafCutRes <- RAINBOWR::MAF.cut(x.0 = genoScoreSubGrRaw,
                                   map.0 = genoMap,
                                   min.MAF = mafThres)

    gastonDataMaf <- gaston::as.bed.matrix(x = mafCutRes$x)
    gastonDataMaf@snps[, c(2, 1, 4)] <- mafCutRes$map

    gastonDataLd <- gaston::LD.thin(x = gastonDataMaf,
                                    threshold = ldThres)

    genoScoreSubGr <- gaston::as.matrix(gastonDataLd)
    genoMapSubGr <- gastonDataLd@snps[, c(2, 1, 4)]
    colnames(genoMapSubGr) <- colnames(mafCutRes$map)
    phenoSubGr <- phenoNow[gbInfo %in% gbNow, , drop = TRUE]
    addGRMSubGr <- addGRM[gbInfo %in% gbNow, gbInfo %in% gbNow]
    nLinesSubGr <- nrow(genoScoreSubGr)


    ## GWAS - Parameters inference
    varListSubGr <- list(addGRMSubGr, diag(nLinesSubGr))
    modelSubGr <- MM4LMM::MMEst(Y = phenoSubGr, X = genoScoreSubGr,
                                VarList = varListSubGr,
                                NbCores = nCores)

    ## GWAS - Tests
    # deltaSubGr <- sapply(modelSubGr, function(i) i$Beta[2])
    testSubGr <- MM4LMM::AnovaTest(modelSubGr,
                                   TestedCombination = matrix(c(0, 1), 1, 2),
                                   NbCores = nCores)
    # waldSubGr <- sapply(testSubGr, function(i) i[1, 1])
    pvalSubGr <- sapply(testSubGr, function(i) i[1, 2])
    scoresSubGr <- -log10(pvalSubGr)

    resDf <- data.frame(genoMapSubGr,
                        scoresSubGr)
    colnames(resDf)[4] <- colnames(phenoNow)
  } else if (length(gbNow) == 2) {
    ### M2 model
    ## Flint and Dent dataset
    genoScoreSubGrRaw <- genoScore[gbInfo %in% gbNow, ]
    mafCutRes <- RAINBOWR::MAF.cut(x.0 = genoScoreSubGrRaw,
                                   map.0 = genoMap,
                                   min.MAF = mafThres)
    genoScoreSubGr <- mafCutRes$x
    genoMapSubGr <- mafCutRes$map
    genoScoreWithGenBackSubGr <- genoScoreWithGenBack[gbInfo %in% gbNow, genoMapSubGr$marker]
    phenoSubGr <- phenoNow[gbInfo %in% gbNow, , drop = TRUE]
    addGRMSubGr <- addGRM[gbInfo %in% gbNow, gbInfo %in% gbNow]
    nLinesSubGr <- nrow(genoScoreWithGenBackSubGr)
    gbInfoSubGr <- gbInfo[gbInfo %in% gbNow]


    ## Null M2 general polygenic model
    # varListSubGr0 <- list(addGRMSubGr)
    # for (gbNowEach in gbNow) {
    #   addGRMSubGrDifVarNow <- matrix(data = 0,
    #                                  nrow = nLinesSubGr,
    #                                  ncol = nLinesSubGr,
    #                                  dimnames = dimnames(addGRMSubGr))
    #   addGRMSubGrDifVarNow[gbInfoSubGr %in% gbNowEach, gbInfoSubGr %in% gbNowEach] <-
    #     addGRMSubGr[gbInfoSubGr %in% gbNowEach, gbInfoSubGr %in% gbNowEach]
    #
    #   varListSubGr0 <- c(varListSubGr0, list(addGRMSubGrDifVarNow))
    # }
    #
    # varListSubGrNullModel <- c(varListSubGr0, list(diag(nLinesSubGr)))
    #
    #
    # modelNullSubGr <- MM4LMM::MMEst(Y = phenoSubGr,
    #                                 Cofactor = gbInfoSubGr,
    #                                 VarList = varListSubGrNullModel,
    #                                 NbCores = nCores)
    #
    # addGRMSubGrNullModel <- Reduce(f = "+",
    #                                x = lapply(X = 1:(length(varListSubGrNullModel) - 1),
    #                                           FUN = function(i) {
    #                                             varListSubGrNullModel[[i]] *
    #                                               modelNullSubGr[[1]]$Sigma2[i]
    #                                           }))

    # ZETASubGrNullModel <- list(A = list(Z = diag(nrow(addGRMSubGr)), K = addGRMSubGr))
    # for (gbNowEach in gbNow) {
    #   addGRMSubGrDifVarNow <- matrix(data = 0,
    #                                  nrow = nLinesSubGr,
    #                                  ncol = nLinesSubGr,
    #                                  dimnames = dimnames(addGRMSubGr))
    #   addGRMSubGrDifVarNow[gbInfoSubGr %in% gbNowEach, gbInfoSubGr %in% gbNowEach] <-
    #     addGRMSubGr[gbInfoSubGr %in% gbNowEach, gbInfoSubGr %in% gbNowEach]
    #
    #   ZETASubGrNullModel <- c(ZETASubGrNullModel,
    #                           list(list(Z = diag(nrow(addGRMSubGrDifVarNow)),
    #                                     K = addGRMSubGrDifVarNow)))
    # }
    # modelNullSubGr <- RAINBOWR::EM3.general(y = phenoSubGr,
    #                                         X0 = model.matrix(~ factor(gbInfoSubGr)),
    #                                         ZETA = ZETASubGrNullModel,
    #                                         package = "gaston")
    #
    # addGRMSubGrNullModel <- Reduce(f = "+",
    #                                x = lapply(X = 1:length(ZETASubGrNullModel),
    #                                           FUN = function(i) {
    #                                             ZETASubGrNullModel[[i]]$K *
    #                                               modelNullSubGr$weights[i] *
    #                                               modelNullSubGr$Vu
    #                                           }))

    if (adjustAddGRM) {
      adjustGRMRes <- RAINBOWR::adjustGRM(y = phenoSubGr,
                                          X = NULL,
                                          ZETA = list(A = list(Z = diag(nrow(addGRMSubGr)),
                                                               K = addGRMSubGr)),
                                          subpopInfo = gbInfoSubGr,
                                          nSubpop = length(gbNow),
                                          nPcsFindCluster = 10,
                                          include.epistasis = FALSE,
                                          package.MM = "gaston")
      addGRMSubGrNullModel <- adjustGRMRes$ZETAAdjust$Adjust$K
    } else {
      addGRMSubGrNullModel <- addGRMSubGr
    }


    ## GWAS - Parameters inference
    # varListSubGr <- list(addGRMSubGrNullModel / (sum(diag(addGRMSubGrNullModel)) / nLinesSubGr),
    #                      diag(nLinesSubGr))
    varListSubGr <- list(addGRMSubGrNullModel,
                         diag(nLinesSubGr))
    modelSubGr <- MM4LMM::MMEst(Y = phenoSubGr,
                                X = genoScoreWithGenBackSubGr,
                                VarList = varListSubGr,
                                NbCores = nCores)

    ## GWAS - Tests
    gbHeadNow <- stringr::str_sub(string = sort(gbNow), start = 1, end = 1)
    if (identical(gbHeadNow, c("C", "K"))) {
      testNames <- c("All",
                     sort(gbNow),
                     paste(rev(gbHeadNow), collapse = "+"),
                     paste(rev(gbHeadNow), collapse = "-"))
    } else {
      testNames <- c("All",
                     sort(gbNow),
                     paste(gbHeadNow, collapse = "+"),
                     paste(gbHeadNow, collapse = "-"))
    }
    # deltaSubGr <- sapply(X = testList,
    #                     FUN = function(a) {
    #                       sapply(X = modelSubGr,
    #                              FUN = function(i) {
    #                                a %*% i$Beta
    #                              })
    #                     })


    testSubGr <- parallel::mclapply(
      X = modelSubGr,
      FUN = function(modelSubGrEach) {
        effects <- modelSubGrEach$Beta
        effectNames <- names(effects)
        effectNames[1] <- paste0(0, sort(gbHeadNow)[1])
        effectNames <- stringr::str_remove(string = effectNames,
                                           pattern = "Xeffect")

        effectsIncZero <- stringr::str_detect(string = effectNames,
                                              pattern = "0")
        effectsIncTwo <- stringr::str_detect(string = effectNames,
                                             pattern = "2")

        effectsIncGb1 <- stringr::str_detect(string = effectNames,
                                             pattern = sort(gbHeadNow)[1])
        effectsIncGb2 <- stringr::str_detect(string = effectNames,
                                             pattern = sort(gbHeadNow)[2])

        testFirstGrVec <- rep(0, length(effectNames))
        testFirstGrVec[effectsIncTwo & effectsIncGb1] <- 1

        testSecondGrVec <- rep(0, length(effectNames))
        testSecondGrVec[effectsIncZero & effectsIncGb2] <- -1
        testSecondGrVec[effectsIncTwo & effectsIncGb2] <- 1

        testFirstGr <- t(as.matrix(testFirstGrVec))
        testSecondGr <- t(as.matrix(testSecondGrVec))

        existTestFirstGr <- any(testFirstGr != 0)
        existTestSecondGr <- sum(testSecondGr != 0) >= 2
        if (existTestFirstGr & existTestSecondGr) {
          testList <- list(rbind(testFirstGr, testSecondGr),
                           testFirstGr, testSecondGr,
                           testFirstGr + testSecondGr,
                           testFirstGr - testSecondGr)
        } else if (existTestSecondGr) {
          testList <- list(testSecondGr)
        } else if (existTestFirstGr) {
          testList <- list(testFirstGr)
        } else {
          testList <- NULL
        }

        modelSubGrEachList <- list(modelSubGrEach)

        testSubGrEach <- matrix(data = NA,
                                nrow = length(testNames),
                                ncol = 3)
        colnames(testSubGrEach) <- c("Wald", "pval", "df")


        if (!is.null(testList)) {
          testSubGrEach0 <- MM4LMM::AnovaTest(ResMMEst = modelSubGrEachList,
                                              TestedCombination = testList,
                                              NbCores = 1)[[1]]
          if (existTestFirstGr & existTestSecondGr) {
            testSubGrEach <- testSubGrEach0
          } else if (existTestSecondGr) {
            testSubGrEach[3, ] <- testSubGrEach0
          } else if (existTestFirstGr) {
            testSubGrEach[2, ] <- testSubGrEach0
          }
        }

        return(testSubGrEach)
      },
      mc.cores = nCores
    )

    # waldSubGr <- t(sapply(X = testSubGr,
    #                       FUN = function(a) {
    #                         a[, 1]
    #                }))
    pvalSubGr <- t(sapply(X = testSubGr,
                          FUN = function(a) {
                            a[, 2]
                          }))
    # colnames(deltaSubGr) <- colnames(waldSubGr) <- testNames
    colnames(pvalSubGr) <- testNames


    scoresSubGr <- -log10(pvalSubGr)
    resDf <- data.frame(genoMapSubGr,
                        scoresSubGr,
                        check.names = FALSE)
  } else if (length(gbNow) == 3) {
    ### M3 model
    ## All dataset
    genoScoreWithGenBackSubGr <- genoScoreWithGenBack[gbInfo %in% gbNow, ]
    phenoSubGr <- phenoNow[gbInfo %in% gbNow, , drop = TRUE]
    addGRMSubGr <- addGRM[gbInfo %in% gbNow, gbInfo %in% gbNow]
    nLinesSubGr <- nrow(genoScoreWithGenBackSubGr)
    gbInfoSubGr <- gbInfo[gbInfo %in% gbNow]


    ## Null M3 general polygenic model
    gbNowSubsets <- c(lapply(X = gbNow, FUN = function(x) x),
                      list(gbNow[2:3], gbNow[c(1, 3)],
                           gbNow[1:3]))

    # varListSubGr0 <- NULL
    # for (gbNowEach in gbNowSubsets) {
    #   addGRMSubGrDifVarNow <- matrix(data = 0,
    #                                  nrow = nLinesSubGr,
    #                                  ncol = nLinesSubGr,
    #                                  dimnames = dimnames(addGRMSubGr))
    #   addGRMSubGrDifVarNow[gbInfoSubGr %in% gbNowEach, gbInfoSubGr %in% gbNowEach] <-
    #     addGRMSubGr[gbInfoSubGr %in% gbNowEach, gbInfoSubGr %in% gbNowEach]
    #
    #   varListSubGr0 <- c(varListSubGr0, list(addGRMSubGrDifVarNow))
    # }
    #
    # varListSubGrNullModel <- c(varListSubGr0, list(diag(nLinesSubGr)))
    #
    #
    # modelNullSubGr <- MM4LMM::MMEst(Y = phenoSubGr,
    #                                 Cofactor = gbInfoSubGr,
    #                                 VarList = varListSubGrNullModel,
    #                                 NbCores = nCores)
    #
    # addGRMSubGrNullModel <- Reduce(f = "+",
    #                                x = lapply(X = 1:(length(varListSubGrNullModel) - 1),
    #                                           FUN = function(i) {
    #                                             varListSubGrNullModel[[i]] *
    #                                               modelNullSubGr[[1]]$Sigma2[i]
    #                                           }))

    ZETASubGrNullModel <- list(A = list(Z = diag(nrow(addGRMSubGr)), K = addGRMSubGr))
    for (gbNowEach in gbNowSubsets) {
      addGRMSubGrDifVarNow <- matrix(data = 0,
                                     nrow = nLinesSubGr,
                                     ncol = nLinesSubGr,
                                     dimnames = dimnames(addGRMSubGr))
      addGRMSubGrDifVarNow[gbInfoSubGr %in% gbNowEach, gbInfoSubGr %in% gbNowEach] <-
        addGRMSubGr[gbInfoSubGr %in% gbNowEach, gbInfoSubGr %in% gbNowEach]

      ZETASubGrNullModel <- c(ZETASubGrNullModel,
                              list(list(Z = diag(nrow(addGRMSubGrDifVarNow)),
                                        K = addGRMSubGrDifVarNow)))
    }
    modelNullSubGr <- RAINBOWR::EM3.general(y = phenoSubGr,
                                            X0 = model.matrix(~ factor(gbInfoSubGr)),
                                            ZETA = ZETASubGrNullModel,
                                            package = "gaston")

    addGRMSubGrNullModel <- Reduce(f = "+",
                                   x = lapply(X = 1:length(ZETASubGrNullModel),
                                              FUN = function(i) {
                                                ZETASubGrNullModel[[i]]$K *
                                                  modelNullSubGr$weights[i] *
                                                  modelNullSubGr$Vu
                                              }))


    ## GWAS - Parameters inference
    varListSubGr <- list(addGRMSubGrNullModel / (sum(diag(addGRMSubGrNullModel)) / nLinesSubGr),
                         diag(nLinesSubGr))
    modelSubGr <- MM4LMM::MMEst(Y = phenoSubGr,
                                X = genoScoreWithGenBackSubGr,
                                VarList = varListSubGr,
                                NbCores = nCores)

    ## GWAS - Tests
    gbHeadNow <- stringr::str_sub(string = sort(gbNow), start = 1, end = 1)
    testNames <- c("All",
                   sort(gbNow),
                   paste(gbHeadNow[1:2], collapse = "+"),
                   paste(gbHeadNow[1:2], collapse = "-"),
                   paste(gbHeadNow[2:3], collapse = "+"),
                   paste(gbHeadNow[2:3], collapse = "-"),
                   paste(gbHeadNow[c(3, 1)], collapse = "+"),
                   paste(gbHeadNow[c(3, 1)], collapse = "-"),
                   paste0(gbHeadNow[1], "-",
                          gbHeadNow[2], "-",
                          gbHeadNow[3]),
                   paste0(gbHeadNow[1], "-",
                          gbHeadNow[2], "+",
                          gbHeadNow[3]),
                   paste0(gbHeadNow[1], "+",
                          gbHeadNow[2], "-",
                          gbHeadNow[3]),
                   paste0(gbHeadNow[1], "+",
                          gbHeadNow[2], "+",
                          gbHeadNow[3])
    )


    # deltaSubGr <- sapply(X = testList,
    #                     FUN = function(a) {
    #                       sapply(X = modelSubGr,
    #                              FUN = function(i) {
    #                                a %*% i$Beta
    #                              })
    #                     })

    testSubGr <- parallel::mclapply(
      X = modelSubGr,
      FUN = function(modelSubGrEach) {
        effects <- modelSubGrEach$Beta
        effectNames <- names(effects)
        effectNames[1] <- paste0(0, sort(gbHeadNow)[1])
        effectNames <- stringr::str_remove(string = effectNames,
                                           pattern = "Xeffect")

        effectsIncZero <- stringr::str_detect(string = effectNames,
                                              pattern = "0")
        effectsIncTwo <- stringr::str_detect(string = effectNames,
                                             pattern = "2")

        effectsIncGb1 <- stringr::str_detect(string = effectNames,
                                             pattern = sort(gbHeadNow)[1])
        effectsIncGb2 <- stringr::str_detect(string = effectNames,
                                             pattern = sort(gbHeadNow)[2])
        effectsIncGb3 <- stringr::str_detect(string = effectNames,
                                             pattern = sort(gbHeadNow)[3])

        testFirstGrVec <- rep(0, length(effectNames))
        testFirstGrVec[effectsIncTwo & effectsIncGb1] <- 1

        testSecondGrVec <- rep(0, length(effectNames))
        testSecondGrVec[effectsIncZero & effectsIncGb2] <- -1
        testSecondGrVec[effectsIncTwo & effectsIncGb2] <- 1

        testThirdGrVec <- rep(0, length(effectNames))
        testThirdGrVec[effectsIncZero & effectsIncGb3] <- -1
        testThirdGrVec[effectsIncTwo & effectsIncGb3] <- 1

        testFirstGr <- t(as.matrix(testFirstGrVec))
        testSecondGr <- t(as.matrix(testSecondGrVec))
        testThirdGr <- t(as.matrix(testThirdGrVec))

        existTestFirstGr <- any(testFirstGr != 0)
        existTestSecondGr <- sum(testSecondGr != 0) >= 2
        existTestThirdGr <- sum(testThirdGr != 0) >= 2
        if (existTestFirstGr & existTestSecondGr & existTestThirdGr) {
          testList <- list(
            rbind(testFirstGr, testSecondGr, testThirdGr),
            testFirstGr, testSecondGr, testThirdGr,
            testFirstGr + testSecondGr,
            testFirstGr - testSecondGr,
            testSecondGr + testThirdGr,
            testSecondGr - testThirdGr,
            testFirstGr + testThirdGr,
            testFirstGr - testThirdGr,
            testFirstGr - testSecondGr - testThirdGr,
            testFirstGr - testSecondGr + testThirdGr,
            testFirstGr + testSecondGr - testThirdGr,
            testFirstGr + testSecondGr + testThirdGr
          )
        } else if (existTestFirstGr & existTestSecondGr) {
          testList <- list(rbind(testFirstGr, testSecondGr),
                           testFirstGr, testSecondGr,
                           testFirstGr + testSecondGr,
                           testFirstGr - testSecondGr)
        } else if (existTestSecondGr & existTestThirdGr) {
          testList <- list(rbind(testSecondGr, testThirdGr),
                           testSecondGr, testThirdGr,
                           testSecondGr + testThirdGr,
                           testSecondGr - testThirdGr)
        } else if (existTestFirstGr & existTestThirdGr) {
          testList <- list(rbind(testFirstGr, testThirdGr),
                           testFirstGr, testThirdGr,
                           testFirstGr + testThirdGr,
                           testFirstGr - testThirdGr)
        } else if (existTestFirstGr) {
          testList <- list(testFirstGr)
        } else if (existTestSecondGr) {
          testList <- list(testSecondGr)
        } else if (existTestThirdGr) {
          testList <- list(testThirdGr)
        } else {
          testList <- NULL
        }

        modelSubGrEachList <- list(modelSubGrEach)

        testSubGrEach <- matrix(data = NA,
                                nrow = length(testNames),
                                ncol = 3)
        colnames(testSubGrEach) <- c("Wald", "pval", "df")


        if (!is.null(testList)) {
          testSubGrEach0 <- MM4LMM::AnovaTest(ResMMEst = modelSubGrEachList,
                                              TestedCombination = testList,
                                              NbCores = 1)[[1]]


          if (existTestFirstGr & existTestSecondGr & existTestThirdGr) {
            testSubGrEach <- testSubGrEach0
          } else if (existTestFirstGr & existTestSecondGr) {
            testSubGrEach[c(1, 2, 3, 5, 6), ] <- testSubGrEach0
          } else if (existTestSecondGr & existTestThirdGr) {
            testSubGrEach[c(1, 3, 4, 7, 8), ] <- testSubGrEach0
          } else if (existTestFirstGr & existTestThirdGr) {
            testSubGrEach[c(1, 2, 4, 9, 10), ] <- testSubGrEach0
          } else if (existTestFirstGr) {
            testSubGrEach[2, ] <- testSubGrEach0
          } else if (existTestSecondGr) {
            testSubGrEach[3, ] <- testSubGrEach0
          } else if (existTestThirdGr) {
            testSubGrEach[4, ] <- testSubGrEach0
          }
        }

        return(testSubGrEach)
      },
      mc.cores = nCores
    )
    # waldSubGr <- t(sapply(X = testSubGr,
    #                       FUN = function(a) {
    #                         a[, 1]
    #                }))
    pvalSubGr <- t(sapply(X = testSubGr,
                          FUN = function(a) {
                            a[, 2]
                          }))
    # colnames(deltaSubGr) <- colnames(waldSubGr) <- testNames
    colnames(pvalSubGr) <- testNames


    scoresSubGr <- -log10(pvalSubGr)
    resDf <- data.frame(genoMap,
                        scoresSubGr,
                        check.names = FALSE)

  }



  return(resDf)
}
