##########################################################################################
######  Title: 1.3_Soybean_DLIPG_function_to_perform_GWAS_by_HB_interaction_model     ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                            ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2022/06/29 (Created), 2022/07/13 (Last Updated)                       ######
##########################################################################################




performGWASHBInt <- function(genoGWAS,
                             phenoGWASNow,
                             haploBlockListDf,
                             mapHaploBlock = NULL,
                             interactionKernelOrigin = NULL,
                             addGRM = NULL,
                             addGRMSubGr = NULL,
                             gbInfo,
                             gbInfoUsed = NULL,
                             gbNow,
                             interactionMatMethod = "HBxGB",
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

  if (is.null(interactionKernelOrigin)) {
    interactionKernelOrigin <- addGRMSubGr
  } else {
    stopifnot(nrow(interactionKernelOrigin) == nLinesSubGr)
    stopifnot(ncol(interactionKernelOrigin) == nLinesSubGr)
  }


  nSubGroup <- length(gbNow)

  if (nSubGroup >= 2) {
    ZETASubGrNullModel <- list(A = list(Z = design.Z(pheno.labels = phenoGWASSubGr[, 1],
                                                     geno.names = rownames(addGRMSubGr)),
                                        K = addGRMSubGr))
    if (is.null(gbInfoUsed)) {
      gbInfoUsed <- adegenet::find.clusters(x = ZETASubGrNullModel$A$Z %*%
                                              tcrossprod(ZETASubGrNullModel$A$K,
                                                         ZETASubGrNullModel$A$Z),
                                            n.pca = nPcsDapc,
                                            n.clust = nSubGroup)$grp
    } else {
      stopifnot(length(gbInfoUsed) == nLinesSubGr)
    }


    if (adjustAddGRM) {
      includeEpistasis <- (interactionMatMethod == "HBxGB") & is.null(paramForIntGWAS)
      adjustGRMRes <- RAINBOWR::adjustGRM(y = phenoSubGr,
                                          X = NULL,
                                          ZETA = ZETASubGrNullModel,
                                          subpopInfo = gbInfoUsed,
                                          nSubpop = nSubGroup,
                                          nPcsFindCluster = nPcsDapc,
                                          include.epistasis = includeEpistasis,
                                          package.MM = "gaston")
      addGRMSubGrNullModel <- adjustGRMRes$ZETAAdjust$Adjust$K
    } else {
      addGRMSubGrNullModel <- addGRMSubGr
    }
  } else {
    includeEpistasis <- (interactionMatMethod == "HBxGB") & is.null(paramForIntGWAS)
    if (adjustAddGRM & includeEpistasis) {
      adjustGRMRes <- RAINBOWR::adjustGRM(y = phenoSubGr,
                                          X = NULL,
                                          ZETA = ZETASubGrNullModel,
                                          subpopInfo = gbInfoUsed,
                                          nSubpop = nSubGroup,
                                          nPcsFindCluster = nPcsDapc,
                                          include.epistasis = includeEpistasis,
                                          package.MM = "gaston")
      addGRMSubGrNullModel <- adjustGRMRes$ZETAAdjust$Adjust$K
    } else {
      addGRMSubGrNullModel <- addGRMSubGr
    }
  }

  ZETASubGr <- list(A = list(Z = design.Z(pheno.labels = phenoGWASSubGr[, 1],
                                          geno.names = rownames(addGRMSubGrNullModel)),
                             K = addGRMSubGrNullModel))

  rm(genoGWAS); rm(phenoGWASNow); rm(addGRM)
  gc(reset = TRUE); gc(reset = TRUE)

  if (nSubGroup <= 1) {
    covariateFactor <- NULL
  } else {
    covariateFactor <- gbInfoUsed
  }

  if (interactionMatMethod == "HB") {
    resDf <- RAINBOWR::RGWAS.multisnp(pheno = phenoGWASSubGr,
                                      geno = genoGWASSubGr,
                                      ZETA = ZETASubGr,
                                      gene.set = haploBlockListDf,
                                      map.gene.set = mapHaploBlock,
                                      package.MM = "gaston",
                                      test.method = "LR",
                                      kernel.method = "linear",
                                      test.effect = "additive",
                                      n.PC = 0,
                                      n.core = nCores,
                                      parallel.method = "mclapply",
                                      plot.qq = FALSE,
                                      plot.Manhattan = FALSE,
                                      saveName = NULL,
                                      thres = FALSE,
                                      skip.check = TRUE,
                                      covariate.factor = covariateFactor,
                                      count = TRUE,
                                      verbose = TRUE)
  } else {
    if (!is.null(paramForIntGWAS)) {
      paramForIntGWAS <- as.numeric(paramForIntGWAS)
      pcaResInteractionKernel <- prcomp(x = interactionKernelOrigin)
      interactionKernel <- tcrossprod(pcaResInteractionKernel$x[, 1:paramForIntGWAS],
                                      pcaResInteractionKernel$rotation[, 1:paramForIntGWAS])
    } else {
      interactionKernel <- interactionKernelOrigin
    }

    resDf <- RAINBOWR::RGWAS.multisnp.interaction(pheno = phenoGWASSubGr,
                                                  geno = genoGWASSubGr,
                                                  ZETA = ZETASubGr,
                                                  interaction.kernel = interactionKernel,
                                                  include.interaction.kernel.null = FALSE,
                                                  include.interaction.with.gb.null = FALSE,
                                                  gene.set = haploBlockListDf,
                                                  map.gene.set = mapHaploBlock,
                                                  package.MM = "gaston",
                                                  test.method = "LR",
                                                  kernel.method = "linear",
                                                  test.effect = "additive",
                                                  n.PC = 0,
                                                  n.core = nCores,
                                                  parallel.method = "mclapply",
                                                  plot.qq = FALSE,
                                                  plot.Manhattan = FALSE,
                                                  saveName = NULL,
                                                  thres = FALSE,
                                                  skip.check = TRUE,
                                                  covariate.factor = covariateFactor,
                                                  count = TRUE,
                                                  verbose = TRUE)
  }


  return(resDf)
}
