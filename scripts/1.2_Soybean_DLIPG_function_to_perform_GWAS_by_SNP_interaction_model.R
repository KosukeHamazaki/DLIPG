############################################################################################
######  Title: 1.2_Soybean_DLIPG_function_to_perform_GWAS_by_SNP_interaction_model    ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                              ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo     ######
######  Date: 2022/06/27 (Created), 2022/07/13 (Last Updated)                         ######
############################################################################################




performGWASSNPInt <- function(genoGWAS,
                              phenoGWASNow,
                              addGRM = NULL,
                              addGRMSubGr = NULL,
                              gbInfo,
                              gbNow,
                              interactionMatMethod = "PCA",
                              paramForIntGWAS = 1,
                              nPcsDapc = 10,
                              nCores = 1,
                              mafThres = 0.05,
                              ldThres = 0.95,
                              adjustAddGRM = FALSE) {
  nLines <- ncol(genoGWAS) - 3
  nSNPs <- nrow(genoGWAS)

  phenoGWASSubGr <- phenoGWASNow[gbInfo %in% gbNow, , drop = FALSE]
  genoGWASSubGrRaw <- genoGWAS[, c(rep(TRUE, 3), gbInfo %in% gbNow)]
  mafCutRes <- RAINBOWR::MAF.cut(x.0 = t(genoGWASSubGrRaw[, -c(1:3)]),
                                 map.0 = genoGWASSubGrRaw[, 1:3],
                                 min.MAF = mafThres)

  gastonDataMaf <- gaston::as.bed.matrix(x = mafCutRes$x + 1)
  gastonDataMaf@snps[, c(2, 1, 4)] <- mafCutRes$map

  gastonDataLd <- gaston::LD.thin(x = gastonDataMaf,
                                  threshold = ldThres)
  mrkNamesNow <- gastonDataLd@snps[, 2]

  genoGWASSubGr <- genoGWASSubGrRaw[genoGWASSubGrRaw$marker %in% mrkNamesNow, ]
  nLinesSubGr <- nrow(phenoGWASSubGr)
  phenoSubGr <- phenoGWASSubGr[, 2]
  gbInfoSubGr <- gbInfo[gbInfo %in% gbNow]

  if (is.null(addGRMSubGr)) {
    stopifnot(!is.null(addGRM))
    addGRMSubGr <- addGRM[gbInfo %in% gbNow, gbInfo %in% gbNow]
  } else {
    stopifnot(nrow(addGRMSubGr) == nLinesSubGr)
    stopifnot(ncol(addGRMSubGr) == nLinesSubGr)
  }


  if ((paramForIntGWAS == "true") | (interactionMatMethod == "PCA")) {
    nSubGroup <- length(gbNow)
  } else {
    nSubGroup <- as.numeric(paramForIntGWAS)
  }

  if (nSubGroup >= 2) {
    ZETASubGrNullModel <- list(A = list(Z = design.Z(pheno.labels = phenoGWASSubGr[, 1],
                                                     geno.names = rownames(addGRMSubGr)),
                                        K = addGRMSubGr))
    if (paramForIntGWAS != "true") {
      gbInfoUsed <- adegenet::find.clusters(x = ZETASubGrNullModel$A$Z %*%
                                              tcrossprod(ZETASubGrNullModel$A$K,
                                                         ZETASubGrNullModel$A$Z),
                                            n.pca = nPcsDapc,
                                            n.clust = nSubGroup)$grp
    } else {
      gbInfoUsed <- factor(gbInfoSubGr)
    }


    if (adjustAddGRM) {
      # if (nSubGroup == 2) {
      #   gbNowSubsets <- unique(gbInfoUsed)
      # } else if (nSubGroup == 3) {
      #   gbNowSubsets <- c(lapply(X = unique(gbInfoUsed), FUN = function(x) x),
      #                     list(unique(gbInfoUsed)[2:3], unique(gbInfoUsed)[c(1, 3)]))
      # }
      #
      #
      # for (gbNowEach in gbNowSubsets) {
      #   addGRMSubGrDifVarNow <- matrix(data = 0,
      #                                  nrow = nLinesSubGr,
      #                                  ncol = nLinesSubGr,
      #                                  dimnames = dimnames(addGRMSubGr))
      #   addGRMSubGrDifVarNow[gbInfoUsed %in% gbNowEach, gbInfoUsed %in% gbNowEach] <-
      #     addGRMSubGr[gbInfoUsed %in% gbNowEach, gbInfoUsed %in% gbNowEach]
      #
      #   ZETASubGrNullModel <- c(ZETASubGrNullModel,
      #                           list(list(Z = design.Z(pheno.labels = phenoGWASSubGr[, 1],
      #                                                  geno.names = rownames(addGRMSubGrDifVarNow)),
      #                                     K = addGRMSubGrDifVarNow)))
      # }
      # modelNullSubGr <- RAINBOWR::EM3.general(y = phenoSubGr,
      #                                         X0 = model.matrix(~ factor(gbInfoUsed)),
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


      adjustGRMRes <- RAINBOWR::adjustGRM(y = phenoSubGr,
                                          X = NULL,
                                          ZETA = ZETASubGrNullModel,
                                          subpopInfo = gbInfoUsed,
                                          nSubpop = nSubGroup,
                                          nPcsFindCluster = nPcsDapc,
                                          include.epistasis = FALSE,
                                          package.MM = "gaston")
      addGRMSubGrNullModel <- adjustGRMRes$ZETAAdjust$Adjust$K
    } else {
      addGRMSubGrNullModel <- addGRMSubGr
    }
  } else {
    addGRMSubGrNullModel <- addGRMSubGr
  }

  ZETASubGr <- list(A = list(Z = design.Z(pheno.labels = phenoGWASSubGr[, 1],
                                          geno.names = rownames(addGRMSubGrNullModel)),
                             K = addGRMSubGrNullModel))

  rm(genoGWAS); rm(phenoGWASNow); rm(addGRM)
  gc(reset = TRUE); gc(reset = TRUE)

  if (interactionMatMethod == "NI") {
    if (nSubGroup <= 1) {
      covariateFactor <- NULL
    } else {
      covariateFactor <- gbInfoUsed
    }


    resDf <- RAINBOWR::RGWAS.normal(pheno = phenoGWASSubGr,
                                    geno = genoGWASSubGr,
                                    ZETA = ZETASubGr,
                                    package.MM = "gaston",
                                    n.PC = 0,
                                    P3D = TRUE,
                                    n.core = nCores,
                                    parallel.method = "mclapply",
                                    plot.qq = FALSE,
                                    plot.Manhattan = FALSE,
                                    saveName = NULL,
                                    thres = FALSE,
                                    skip.check = TRUE,
                                    covariate.factor = covariateFactor,
                                    count = FALSE,
                                    verbose = FALSE)
  } else {
    if (interactionMatMethod == "PCA") {
      nInteractionElement <- paramForIntGWAS
      nPcsCovariate <- paramForIntGWAS

      # covariateFactor <- NULL
      if (nSubGroup <= 1) {
        covariateFactor <- NULL
      } else {
        covariateFactor <- gbInfoUsed
      }

      nInteractionGroup <- 1
      interactionGroup <- NULL
      covariate <- interactionWithSNPs <- NULL
    } else {
      nPcsCovariate <- 0

      if (paramForIntGWAS == "true") {
        nInteractionGroup <- length(gbNow)
        # interactionGroup <- factor(gbInfo[gbInfo %in% gbNow])
      } else {
        nInteractionGroup <- as.numeric(paramForIntGWAS)
        # interactionGroup <- adegenet::find.clusters(x = ZETASubGr$A$Z %*%
        #                                               tcrossprod(ZETASubGr$A$K,
        #                                                          ZETASubGr$A$Z),
        #                                             n.pca = nPcsDapc,
        #                                             n.clust = nInteractionGroup)$grp
      }

      if (interactionMatMethod == "DAPC") {
        nInteractionElement <- nInteractionGroup - 1

        dapcRes <- adegenet::dapc(x = addGRMSubGrNullModel,
                                  grp = gbInfoUsed,
                                  n.pca = nPcsDapc,
                                  n.da = nlevels(gbInfoUsed) - 1)
        interactionWithSNPs0 <- dapcRes$ind.coord


        if (ncol(interactionWithSNPs0) > nInteractionElement) {
          interactionWithSNPs <- interactionWithSNPs0[, 1:nInteractionElement, drop = FALSE]
        } else {
          interactionWithSNPs <- interactionWithSNPs0
          nInteractionElement <- ncol(interactionWithSNPs)
        }

        covariate <- interactionWithSNPs
      } else if (interactionMatMethod == "GROUP") {
        nInteractionElement <- nInteractionGroup
        covariate <- interactionWithSNPs <- NULL
      }


      if (nInteractionGroup >= 2) {
        covariateFactor <- interactionGroup <- gbInfoUsed
      } else {
        interactionGroup <- covariateFactor <- NULL
        nInteractionElement <- 0
      }
    }



    resDf <- RAINBOWR::RGWAS.normal.interaction(pheno = phenoGWASSubGr,
                                                geno = genoGWASSubGr,
                                                ZETA = ZETASubGr,
                                                package.MM = "gaston",
                                                interaction.with.SNPs = interactionWithSNPs,
                                                interaction.mat.method = interactionMatMethod,
                                                n.interaction.element = nInteractionElement,
                                                interaction.group = interactionGroup,
                                                n.interaction.group = nInteractionGroup,
                                                interaction.group.method = "find.clusters",
                                                n.PC.dapc = nPcsDapc,
                                                test.method.interaction = "simultaneous",
                                                n.core = nCores,
                                                parallel.method = "mclapply",
                                                n.PC = nPcsCovariate,
                                                covariate = covariate,
                                                covariate.factor = covariateFactor,
                                                P3D = TRUE,
                                                thres = FALSE,
                                                skip.check = TRUE,
                                                plot.Manhattan = FALSE,
                                                plot.qq = FALSE,
                                                count = FALSE,
                                                verbose = FALSE)
    if (nInteractionElement != 0) {
      resDf <- resDf$D$All
    }
    # }
  }


  return(resDf)
}
