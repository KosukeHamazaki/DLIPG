##############################################################################################################
######  Title: 0.6_Soybean_DLIPG_Simulate_phenotypic_data_including_haplotype_based_QTLs_with_kmeans    ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                                                ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo                       ######
######  Date: 2023/01/14 (Created), 2023/01/15 (Last Updated)                                           ######
##############################################################################################################



###### 1. Settings ######
##### 1.0. Reset workspace ######
# rm(list=ls())



##### 1.1. Setting working directory to the "DLIGB" directory #####
cropName <- "Soybean"
project <- "DLIPG"
os <- osVersion
isRproject <- function(path = getwd()) {
  files <- list.files(path)
  if (length(grep(".Rproj", files)) >= 1) {
    out <- TRUE
  } else {
    out <-  FALSE
  }
  return(out)
}

if (!isRproject()) {
  if (stringr::str_detect(string = os, pattern = "mac")) {
    dirResearchBase <- "/Volumes/HD-PGF-A/research/"  ### for mac OS
    dirScriptBase <- "~/GitHub/research_secret/"
  } else if (stringr::str_detect(string = os, pattern = "Ubuntu 18")) {
    dirResearchBase <- "/media/hamazaki/HDD2/research/"   ### for Ubuntu SGM
    dirScriptBase <- "~/GitHub/research_secret/"
  } else if (stringr::str_detect(string = os, pattern = "Ubuntu 20.04.2")) {
    dirResearchBase <- "/media/hamazaki-64core/39F16BAC190DFF39/research/"     ### for Ubuntu 64
    dirScriptBase <- "~/GitHub/research_secret/"
  } else if (stringr::str_detect(string = os, pattern = "Ubuntu 20.04.4")) {
    dirResearchBase <- "/media/hamazaki/HDD1/research/"     ### for Ubuntu DC1
    dirScriptBase <- "~/GitHub/research_secret/"
  } else {
    stop("Which type of work-station do you use?")
  }

  setwd(paste0(dirResearchBase, cropName, "/Project/", project))
}

scriptIDGb <- "0.1"
scriptIDHaploBlock <- "0.4"
scriptIDPca <- "0.5"
scriptID <- "0.6"





##### 1.2. Setting some parameters #####
dirMidDLIGBBase <- "midstream/"
nCores <- 10
scenarioNo <- 3
trialNo <- 1
overWrite <- FALSE
nIter <- 100


inputParams <- list(
  nIter = nIter,
  nPolyGenesCommon = 4000,
  nPolyGenesAncestry = 0,
  nPolyGenesGroup = c(500, 500, 500,
                      0, 0, 0),
  nQtlsCommon = 1,
  nQtlsAncestry = 0,
  nQtlsGroup = c(0, 0, 0,
                 0, 0, 0),
  nQtlsGbs = 0,
  nQtlsPg = 2,
  nHaploQtlsCommon = 1,
  nHaploQtlsAncestry = 0,
  nHaploQtlsGroup = c(0, 0, 0,
                      0, 0, 0),
  nHaploQtlsGbs = 0,
  nHaploQtlsPg = 2,
  nMarkersInBlocksMin = 5,
  qtlCandNamesGiven = NULL,
  effectsAncestryType = NULL,
  effectsGroupType = matrix(data = c(1, 0, 0,
                                     0, 1, 0,
                                     0, 0, 1,
                                     1, -1, 0,
                                     0, 1, -1,
                                     -1, 0, 1),
                            nrow = 6,
                            ncol = 3,
                            byrow = TRUE,
                            dimnames = list(Type = c("China", "Japan", "Korea",
                                                     "C-J", "J-K", "K-C"),
                                            Group = c("China", "Japan", "Korea"))),
  giveQtlCandNames = FALSE,
  nullChrNos = rep(1:20, each = nIter / 20),
  effectsQtlsType = "same",
  heritQtls = 0.5,
  herit = 0.75
)



dirSimData <- paste0(dirMidDLIGBBase,
                     scriptID, "_simulated_data/")
if (!dir.exists(dirSimData)) {
  dir.create(dirSimData)
}

dirScenario <- paste0(dirSimData,
                      "Scenario_", scenarioNo, "/")
if (!dir.exists(dirScenario)) {
  dir.create(dirScenario)
}


dirTrial <- paste0(dirScenario,
                   "Trial_", trialNo, "/")
if (!dir.exists(dirTrial)) {
  dir.create(dirTrial)
}

filesInTrial <- list.files(dirTrial)


if ((length(filesInTrial) == 0) | overWrite) {
  fileParamsSimData <- paste0(dirTrial, scriptID,
                              "_", project, "_simulation_parameters.RData")
  saveRDS(object = inputParams,
          file = fileParamsSimData)






  ##### 1.3. Import packages #####
  require(data.table)
  require(MASS)
  require(BGLR)
  require(RAINBOWR)
  require(ggplot2)




  dirScript <- "scripts/"
  # dirScript <- paste0(dirScriptBase, cropName, "/Project/", project, "/")
  # source(paste0(dirScript, "1.0_Soybean_DLIPG_function_to_draw_original_pair_plots_for_GWAS_results.R"))
  # source(paste0(dirScript, "scriptName2"))




  ##### 1.4. Project options #####
  options(stringAsFactors = FALSE,
          check.names = FALSE)






  ###### 2. Modification of data ######
  ##### 2.1. Read genetic background data into R #####
  fileNameGenBack <- paste0("data/extra/", scriptIDGb,
                            "_accession_names_with_country_of_origin_selected_by_origin_based_kmedoids_equal_sampling.csv")
  genBack <- data.frame(fread(input = fileNameGenBack, header = TRUE),
                        row.names = 1)
  genBackVecRaw <- genBack[, 1]
  names(genBackVecRaw) <- rownames(genBack)

  genBackUniq <- sort(unique(genBackVecRaw))



  genBackCombAll <- c(lapply(X = genBackUniq, FUN = function(x) x),
                      list(genBackUniq[1:2], genBackUniq[2:3],
                           genBackUniq[c(1, 3)], genBackUniq[1:3]))
  genBackNames <- lapply(X = genBackCombAll,
                         FUN = function(x) {
                           paste(x, collapse = "_")
                         })
  names(genBackCombAll) <- genBackNames



  ##### 2.2. Read genotype map for all genetic backgrounds into R and extract common markers to all population #####
  fileNameGenoMap <- "data/extra/L900_Soybean_a1_chrnum_maf0.025_ms0.05_ChKoJa_300x3_maf0.025_imputed_maf0.025_SNP_geno_map.csv"
  genoMap <- data.frame(fread(input = fileNameGenoMap))
  colnames(genoMap)[1] <- "marker"



  ##### 2.3. Read marker genotype into R and extract common markers to all population #####
  fileNameGenoScore <- paste0("data/genotype/L900_Soybean_a1_chrnum_maf0.025_ms0.05_ChKoJa_300x3_maf0.025_imputed_maf0.025_SNP_geno_012.tsv")

  genoScoreDf <- data.frame(fread(input = fileNameGenoScore,
                                  check.names = FALSE),
                            check.names = FALSE)
  genoScore <- t(genoScoreDf[, -c(1:2)])
  colnames(genoScore) <- genoMap$marker

  genBackVec <- genBackVecRaw[rownames(genoScore)]

  ancestralHaploMat <- matrix(data = rep(genBackVec, ncol(genoScore)),
                              nrow = nrow(genoScore),
                              ncol = ncol(genoScore),
                              byrow = FALSE,
                              dimnames = dimnames(genoScore))

  ancestralHaploDf <- data.frame(ancestralHaploMat,
                                 check.names = FALSE)
  rm(ancestralHaploMat)
  gc(reset = TRUE); gc(reset = TRUE)

  # See(genoScoreCommon)




  ##### 2.4. Read additive GRM for each subgroup into R and construct ZETA #####
  fileNameAddGrmEach <- paste0("data/extra/", scriptIDPca, "_",
                               "adjusted_additive_genomic_relationship_",
                               paste(genBackCombAll[[length(genBackCombAll)]],
                                     collapse = "_"), ".rds")
  addGrmEach <- readRDS(file = fileNameAddGrmEach)

  pcaRes <- prcomp(x = addGrmEach)
  pcs <- pcaRes$x



  ##### 2.5. Read haplotype block list into R #####
  haploBlockDataFrame <- data.frame(
    data.table::fread(input = paste0("data/extra/", scriptIDHaploBlock,
                                     "_haplotype_block_data_for_RAINBOWR.csv"))
  )

  fileNameMrkNamesCommon <- paste0("data/extra/", scriptIDPca, "_",
                                   "marker_names_common_to_all_population_combinations.rds")
  mrkNamesCommon <- readRDS(file = fileNameMrkNamesCommon)

  if (inputParams$giveQtlCandNames) {
    qtlCandNamesGiven <- Reduce(f = intersect,
                                x = list(haploBlockDataFrame$marker,
                                         mrkNamesCommon))
  } else {
    qtlCandNamesGiven <- mrkNamesCommon
  }



  ###### 3. Function to simulate phenotypic data ######
  simulatePhenVals <- function(genoScore,
                               genoMap,
                               ancestralHaploDf,
                               haploBlockDf,
                               gbInfo,
                               gbsMat,
                               nIter = 1,
                               nPolyGenesCommon = 2800,
                               nPolyGenesAncestry = c(600, 600, 0),
                               nPolyGenesGroup = 0,
                               nQtlsCommon = 2,
                               nQtlsAncestry = c(1, 1, 1),
                               nQtlsGroup = 0,
                               nQtlsGbs = 0,
                               nQtlsPg = 1,
                               nHaploQtlsCommon = 2,
                               nHaploQtlsAncestry = c(1, 1, 1),
                               nHaploQtlsGroup = 0,
                               nHaploQtlsGbs = 0,
                               nHaploQtlsPg = 1,
                               nMarkersInBlocksMin = 5,
                               qtlCandNamesGiven = NULL,
                               nullChrNos = NULL,
                               effectsAncestryType = matrix(data = c(1, 0,
                                                                     0, 1,
                                                                     1, -1),
                                                            nrow = 3,
                                                            ncol = 2,
                                                            byrow = TRUE,
                                                            dimnames = list(Type = c("Dent", "Flint", "Divergent"),
                                                                            Ancestry = c("Dent", "Flint"))),
                               effectsGroupType = NULL,
                               effectsQtlsType = "same",
                               heritQtls = 0.4,
                               herit = 0.8,
                               saveAt = NULL,
                               overWrite = FALSE,
                               iterNameHead = "Iteration_",
                               returnAll = FALSE,
                               nCores = 1,
                               seedsMat = NULL,
                               returnSeeds = TRUE,
                               nIterOneBlock = 100,
                               scoring012 = FALSE) {
    ##### 3.1. Some settings (define names) for simulation #####
    gbUniq <- unique(gbInfo)
    mrkNames <- colnames(genoScore)
    nLines <- nrow(genoScore)
    lineNames <- rownames(genoScore)

    if (!scoring012) {
      genoScore <- genoScore + 1
    }

    iterNames <- paste0(iterNameHead, 1:nIter)

    qtlTypeNamesCommon <- "Common"
    qtlTypeNamesAncestry <- "Ancestry"
    if (!is.null(effectsAncestryType)) {
      qtlTypeNamesAncestry <- paste(qtlTypeNamesAncestry,
                                    rownames(effectsAncestryType), sep = "-")
    }
    qtlTypeNamesGroup <- "Group"
    if (!is.null(effectsGroupType)) {
      qtlTypeNamesGroup <- paste(qtlTypeNamesGroup,
                                 rownames(effectsGroupType), sep = "-")
    }
    qtlTypeNamesGb <- paste0("PC", 1:length(nHaploQtlsGbs))
    qtlTypeNamesPg <- "Polygenes"
    qtlTypeNamesAll <- c(qtlTypeNamesCommon, qtlTypeNamesAncestry,
                         qtlTypeNamesGroup, qtlTypeNamesGb, qtlTypeNamesPg)

    polyGeneTypeNamesCommon <- "Common"
    polyGeneTypeNamesAncestry <- "Ancestry"
    if (!is.null(effectsAncestryType)) {
      polyGeneTypeNamesAncestry <- paste(polyGeneTypeNamesAncestry,
                                         rownames(effectsAncestryType), sep = "-")
    }
    polyGeneTypeNamesGroup <- "Group"
    if (!is.null(effectsGroupType)) {
      polyGeneTypeNamesGroup <- paste(polyGeneTypeNamesGroup,
                                      rownames(effectsGroupType), sep = "-")
    }
    polyGeneTypeNamesAll <- c(polyGeneTypeNamesCommon, polyGeneTypeNamesAncestry,
                              polyGeneTypeNamesGroup)


    nHaploQtlsAll <- c(nHaploQtlsCommon, nHaploQtlsAncestry,
                       nHaploQtlsGroup, nHaploQtlsGbs, nHaploQtlsPg)
    names(nHaploQtlsAll) <- qtlTypeNamesAll
    nHaploQtlsCum <- cumsum(nHaploQtlsAll)
    nHaploQtls <- sum(nHaploQtlsAll)

    nQtlsAll <- c(nQtlsCommon, nQtlsAncestry,
                  nQtlsGroup, nQtlsGbs, nQtlsPg)
    names(nQtlsAll) <- qtlTypeNamesAll
    nQtlsCum <- cumsum(nQtlsAll)
    nQtls <- sum(nQtlsAll)

    nPolyGenesAll <- c(nPolyGenesCommon, nPolyGenesAncestry, nPolyGenesGroup)
    names(nPolyGenesAll) <- polyGeneTypeNamesAll
    nPolyGenesCum <- cumsum(nPolyGenesAll)
    nPolyGenes <- sum(nPolyGenesAll)


    tabHaploBlock <- table(haploBlockDf$block)[unique(haploBlockDf$block)]


    if (is.null(nullChrNos)) {
      nullChrNos <- rep(0, nIter)
    } else if (length(nullChrNos) == 1) {
      nullChrNos <- rep(nullChrNos, nIter)
    } else {
      stopifnot(length(nullChrNos) == nIter)
    }
    stopifnot(all(nullChrNos %in% c(0, unique(genoMap$chr))))
    names(nullChrNos) <- iterNames

    nSeedTypes <- 5
    seedTypeNames <- c("Qtls_Markers", "Polygenes_Markers",
                       "Qtls_Effects", "Polygenes_Effects",
                       "Residuals_Effects")

    if (is.null(seedsMat)) {
      nSeeds <- nIter * nSeedTypes
      seedsMat <- matrix(data = sample(x = 1:1e09, size = nSeeds,
                                       replace = FALSE),
                         nrow = nIter,
                         ncol = nSeedTypes,
                         dimnames = list(iterNames,
                                         seedTypeNames))

    } else {
      stopifnot(all(dim(seedsMat) == c(nIter, nSeedTypes)))
      dimnames(seedsMat) <- list(iterNames,
                                 seedTypeNames)

      if (any(is.na(seedsMat))) {
        nSeedsMissing <- sum(is.na(seedsMat))
        seedsMat[is.na(seedsMat)] <- sample(x = 1:1e09, size = nSeedsMissing,
                                            replace = FALSE)
      }
    }




    ##### 3.2. Start simulation #####
    iterNamesList <- split(x = iterNames,
                           f = factor(((1:nIter) - 1) %/% nIterOneBlock + 1))
    resultsAllList <- list()
    for (iterNamesNow in iterNamesList) {
      # iterNamesNow <- iterNamesList[[1]]
      resultsAllList <- c(
        resultsAllList,
        pbmcapply::pbmclapply(X = iterNamesNow,
                              FUN = function(iterName) {
                                # iterName <- iterNamesList[[1]][1]
                                # iterName <- iterNamesNow[1]
                                #### 3.2.1. Check save name ####
                                if (!is.null(saveAt)) {
                                  if (!dir.exists(saveAt)) {
                                    dir.create(saveAt)
                                  }

                                  dirIterSave <- paste0(saveAt, iterName, "/")
                                  if (!dir.exists(dirIterSave)) {
                                    dir.create(dirIterSave)
                                  }

                                  fileNameIterAllRes <- paste0(dirIterSave,
                                                               "simulated_phenotype_all_results_",
                                                               iterName,".rds")
                                } else {
                                  fileNameIterAllRes <- ""
                                }


                                if (!file.exists(fileNameIterAllRes) | overWrite) {
                                  #### 3.2.2. Select markers for QTLs ####
                                  mrkNamesCand <- mrkNames
                                  seeds <- seedsMat[iterName, ]
                                  nullChrNo <- nullChrNos[iterName]

                                  if (!is.null(qtlCandNamesGiven)) {
                                    mrkNamesQtlsCand <- Reduce(f = intersect,
                                                               x = list(mrkNamesCand,
                                                                        qtlCandNamesGiven))
                                  } else {
                                    mrkNamesQtlsCand <- mrkNames
                                  }

                                  if (nullChrNo != 0) {
                                    mrkNamesQtlsCandChr <- genoMap$marker[!(genoMap$chr %in% nullChrNo)]
                                    mrkNamesQtlsCand <- Reduce(f = intersect,
                                                               x = list(mrkNamesQtlsCand,
                                                                        mrkNamesQtlsCandChr))
                                  }



                                  set.seed(seed = seeds["Qtls_Markers"])



                                  ### 3.2.2.1. Select markers for haplotype-based QTLs ###
                                  if (nHaploQtls > 0) {
                                    resampleMarkers <- TRUE

                                    while(resampleMarkers) {
                                      mrkNamesHaploQtls <- sapply(X = 1:nHaploQtls,
                                                                  FUN = function(haploQtlsNo) {
                                                                    resampleMarker <- TRUE
                                                                    while (resampleMarker) {
                                                                      mrkNameHaploQtls <- sample(x = mrkNamesQtlsCand,
                                                                                                 size = 1,
                                                                                                 replace = FALSE)
                                                                      whereInBlockList <- haploBlockDf$marker %in% mrkNameHaploQtls

                                                                      if (any(whereInBlockList)) {
                                                                        blockNameHaploQtls <- haploBlockDf[whereInBlockList, 1]

                                                                        nMarkerInBlocksHaploQtlsNow <- tabHaploBlock[blockNameHaploQtls]

                                                                        if (nMarkerInBlocksHaploQtlsNow >= nMarkersInBlocksMin) {
                                                                          mrkNamesInBlockNow <- haploBlockDf[haploBlockDf$block %in% blockNameHaploQtls, 2]
                                                                          whichHaploQtlInBlock <- which(mrkNameHaploQtls == mrkNamesInBlockNow)

                                                                          genoScoreInBlockNow <- genoScore[, mrkNamesInBlockNow, drop = FALSE]
                                                                          genoScoreHaploQtlsNow <- genoScore[, mrkNameHaploQtls]
                                                                          tableHaploQtlsNow <- table((genoScoreHaploQtlsNow + 1) / 2)
                                                                          minAlleleHaploQtlsNow <- as.numeric(names(which.min(x = tableHaploQtlsNow)))

                                                                          haplosNow <- apply(X = (genoScoreInBlockNow + 1) / 2,
                                                                                             MARGIN = 1,
                                                                                             FUN = paste, collapse = "")
                                                                          tableHaplosNow <- table(haplosNow)

                                                                          tableHaplosWithMinAlleleNow <-
                                                                            tableHaplosNow[(as.numeric(
                                                                              stringr::str_sub(
                                                                                string = names(tableHaplosNow),
                                                                                start = whichHaploQtlInBlock,
                                                                                end = whichHaploQtlInBlock
                                                                              )
                                                                            ) == minAlleleHaploQtlsNow)]

                                                                          resampleMarker <- length(tableHaplosWithMinAlleleNow) == 1
                                                                        }
                                                                      }
                                                                    }


                                                                    return(mrkNameHaploQtls)
                                                                  })
                                      blockNamesHaploQtls <- haploBlockDataFrame[haploBlockDataFrame$marker %in% mrkNamesHaploQtls, 1]
                                      resampleMarkers <- length(unique(blockNamesHaploQtls)) == 1
                                    }

                                    names(mrkNamesHaploQtls) <- rep(qtlTypeNamesAll, nHaploQtlsAll)
                                    mrkNamesHaploQtlsList <- split(x = mrkNamesHaploQtls, f = names(mrkNamesHaploQtls))[qtlTypeNamesAll]
                                    names(mrkNamesHaploQtlsList) <- qtlTypeNamesAll

                                    mrkNamesHaploQtlsList <- lapply(X = mrkNamesHaploQtlsList,
                                                                    FUN = function(mrkNamesHaploQtlsEach) {
                                                                      if (!is.null(mrkNamesHaploQtlsEach)) {
                                                                        mrkNamesHaploQtlsEach <- mrkNamesQtlsCand[mrkNamesQtlsCand %in% mrkNamesHaploQtlsEach]
                                                                        return(mrkNamesHaploQtlsEach)
                                                                      } else {
                                                                        return(NULL)
                                                                      }
                                                                    })

                                    mrkNamesCand <- mrkNamesCand[!(mrkNamesCand %in% mrkNamesHaploQtls)]

                                    blockNamesHaploQtls <- haploBlockDf[haploBlockDf$marker %in% mrkNamesHaploQtls, 1]
                                    mrkNamesInSameBlockWithHaploQtls <- haploBlockDf[haploBlockDf$block %in% blockNamesHaploQtls, 2]
                                    mrkNamesQtlsCand <- mrkNamesQtlsCand[!(mrkNamesQtlsCand %in% mrkNamesInSameBlockWithHaploQtls)]
                                  } else {
                                    mrkNamesHaploQtlsList <- rep(list(NULL), length(qtlTypeNamesAll))
                                    names(mrkNamesHaploQtlsList) <- qtlTypeNamesAll
                                  }




                                  ### 3.2.2.2. Select markers for normal QTLs ###
                                  if (nQtls != 0) {
                                    mrkNamesQtls <- sample(x = mrkNamesQtlsCand,
                                                           size = nQtls,
                                                           replace = FALSE)
                                    names(mrkNamesQtls) <- rep(qtlTypeNamesAll, nQtlsAll)
                                    mrkNamesQtlsList <- split(x = mrkNamesQtls, f = names(mrkNamesQtls))[qtlTypeNamesAll]
                                    names(mrkNamesQtlsList) <- qtlTypeNamesAll

                                    mrkNamesQtlsList <- lapply(X = mrkNamesQtlsList,
                                                               FUN = function(mrkNamesQtlsEach) {
                                                                 if (!is.null(mrkNamesQtlsEach)) {
                                                                   mrkNamesQtlsEach <- mrkNamesQtlsCand[mrkNamesQtlsCand %in% mrkNamesQtlsEach]
                                                                   return(mrkNamesQtlsEach)
                                                                 } else {
                                                                   return(NULL)
                                                                 }
                                                               })

                                    mrkNamesCand <- mrkNamesCand[!(mrkNamesCand %in% mrkNamesQtls)]
                                    mrkNamesQtlsCand <- mrkNamesQtlsCand[!(mrkNamesQtlsCand %in% mrkNamesQtls)]
                                  } else {
                                    mrkNamesQtlsList <- rep(list(NULL), length(qtlTypeNamesAll))
                                    names(mrkNamesQtlsList) <- qtlTypeNamesAll
                                  }







                                  #### 3.2.3. Select markers for polygenes ####
                                  set.seed(seed = seeds["Polygenes_Markers"])



                                  ### 3.2.3.1. Select markers for common polygenes ###
                                  if (nPolyGenes != 0) {
                                    mrkNamesPolyGenes <- sample(x = mrkNamesCand,
                                                                size = nPolyGenes,
                                                                replace = FALSE)
                                    names(mrkNamesPolyGenes) <- rep(polyGeneTypeNamesAll, nPolyGenesAll)
                                    mrkNamesPolyGenesList <- split(x = mrkNamesPolyGenes, f = names(mrkNamesPolyGenes))[polyGeneTypeNamesAll]
                                    names(mrkNamesPolyGenesList) <- polyGeneTypeNamesAll

                                    mrkNamesPolyGenesList <- lapply(X = mrkNamesPolyGenesList,
                                                                    FUN = function(mrkNamesPolyGenesEach) {
                                                                      if (!is.null(mrkNamesPolyGenesEach)) {
                                                                        mrkNamesPolyGenesEach <- mrkNamesCand[mrkNamesCand %in% mrkNamesPolyGenesEach]
                                                                        return(mrkNamesPolyGenesEach)
                                                                      } else {
                                                                        return(NULL)
                                                                      }
                                                                    })

                                    mrkNamesCand <- mrkNamesCand[!(mrkNamesCand %in% mrkNamesPolyGenes)]
                                  } else {
                                    mrkNamesPolyGenesList <- rep(list(NULL), length(polyGeneTypeNamesAll))
                                    names(mrkNamesPolyGenesList) <- polyGeneTypeNamesAll
                                  }





                                  #### 3.2.4. Summary marker names ####
                                  mrkNamesList <-
                                    list(
                                      Qtls = list(Qtls = mrkNamesQtlsList,
                                                  HaploQtls = mrkNamesHaploQtlsList),
                                      Polygenes = mrkNamesPolyGenesList,
                                      All = list(
                                        Original = mrkNames,
                                        Remain = list(All = mrkNamesCand,
                                                      Qtl = mrkNamesQtlsCand)
                                      )
                                    )





                                  #### 3.2.5. Determine QTL effects ####
                                  set.seed(seed = seeds["Qtls_Effects"])



                                  ### 3.2.5.1. Determine haplotype-based QTL effects ###
                                  if (nHaploQtlsCommon != 0) {
                                    if (effectsQtlsType == "same") {
                                      effectsHaploQtlsCommonOri <- rep(1, nHaploQtlsCommon)
                                    } else {
                                      effectsHaploQtlsCommonOri <- rgamma(n = nHaploQtlsCommon,
                                                                          shape = 2)
                                    }
                                    effectsHaploQtlsCommon <- effectsHaploQtlsCommonOri * rep(1, nHaploQtlsCommon)
                                    directionHaploQtlsCommon <- sample(x = c(-1, 1), size = nHaploQtlsCommon, replace = TRUE)
                                    effectsHaploQtlsCommon <- effectsHaploQtlsCommon * directionHaploQtlsCommon

                                    names(effectsHaploQtlsCommon) <- mrkNamesHaploQtlsList[["Common"]]
                                  } else {
                                    effectsHaploQtlsCommon <- NULL
                                  }

                                  if (!is.null(effectsAncestryType)) {
                                    if (effectsQtlsType == "same") {
                                      effectsHaploQtlsAncestryOri <- rep(1, sum(nHaploQtlsAncestry))
                                    } else {
                                      effectsHaploQtlsAncestryOri <- rgamma(n = sum(nHaploQtlsAncestry),
                                                                            shape = 2)
                                    }

                                    effectsHaploQtlsAncestryEach <- sqrt(ncol(effectsAncestryType)) / apply(X = effectsAncestryType, MARGIN = 1,
                                                                                                            FUN = function(x) sqrt(sum(x ^ 2)))
                                    effectsHaploQtlsAncestryVec <- effectsHaploQtlsAncestryOri *
                                      rep(effectsHaploQtlsAncestryEach, nHaploQtlsAncestry)
                                    directionHaploQtlsAncestry <- sample(x = c(-1, 1), size = sum(nHaploQtlsAncestry), replace = TRUE)
                                    effectsHaploQtlsAncestryVec <- effectsHaploQtlsAncestryVec * directionHaploQtlsAncestry
                                    effectsHaploQtlsAncestry <- split(x = effectsHaploQtlsAncestryVec,
                                                                      f = rep(qtlTypeNamesAncestry, nHaploQtlsAncestry))[qtlTypeNamesAncestry]
                                    effectsHaploQtlsAncestry <- lapply(X = 1:length(qtlTypeNamesAncestry),
                                                                       FUN = function(qtlTypeNoAncestryEach) {
                                                                         effectsHaploQtlsAncestryEach <- effectsHaploQtlsAncestry[[qtlTypeNoAncestryEach]]
                                                                         if (!is.null(qtlTypeNoAncestryEach)) {
                                                                           names(effectsHaploQtlsAncestryEach) <- (mrkNamesHaploQtlsList[qtlTypeNamesAncestry])[[qtlTypeNoAncestryEach]]
                                                                         }

                                                                         return(effectsHaploQtlsAncestryEach)
                                                                       })
                                    names(effectsHaploQtlsAncestry) <- qtlTypeNamesAncestry

                                  } else {
                                    effectsHaploQtlsAncestry <- rep(list(NULL), length(qtlTypeNamesAncestry))
                                    names(effectsHaploQtlsAncestry) <- qtlTypeNamesAncestry
                                  }


                                  if (!is.null(effectsGroupType)) {
                                    if (effectsQtlsType == "same") {
                                      effectsHaploQtlsGroupOri <- rep(1, sum(nHaploQtlsGroup))
                                    } else {
                                      effectsHaploQtlsGroupOri <- rgamma(n = sum(nHaploQtlsGroup),
                                                                         shape = 2)
                                    }

                                    effectsHaploQtlsGroupEach <- sqrt(ncol(effectsGroupType)) / apply(X = effectsGroupType, MARGIN = 1,
                                                                                                      FUN = function(x) sqrt(sum(x ^ 2)))
                                    effectsHaploQtlsGroupVec <- effectsHaploQtlsGroupOri *
                                      rep(effectsHaploQtlsGroupEach, nHaploQtlsGroup)
                                    directionHaploQtlsGroup <- sample(x = c(-1, 1), size = sum(nHaploQtlsGroup), replace = TRUE)
                                    effectsHaploQtlsGroupVec <- effectsHaploQtlsGroupVec * directionHaploQtlsGroup
                                    effectsHaploQtlsGroup <- split(x = effectsHaploQtlsGroupVec,
                                                                   f = rep(qtlTypeNamesGroup, nHaploQtlsGroup))[qtlTypeNamesGroup]
                                    effectsHaploQtlsGroup <- lapply(X = 1:length(qtlTypeNamesGroup),
                                                                    FUN = function(qtlTypeNoGroupEach) {
                                                                      effectsHaploQtlsGroupEach <- effectsHaploQtlsGroup[[qtlTypeNoGroupEach]]
                                                                      if (!is.null(qtlTypeNoGroupEach)) {
                                                                        names(effectsHaploQtlsGroupEach) <- (mrkNamesHaploQtlsList[qtlTypeNamesGroup])[[qtlTypeNoGroupEach]]
                                                                      }

                                                                      return(effectsHaploQtlsGroupEach)
                                                                    })
                                    names(effectsHaploQtlsGroup) <- qtlTypeNamesGroup
                                  } else {
                                    effectsHaploQtlsGroup <- rep(list(NULL), length(qtlTypeNamesGroup))
                                    names(effectsHaploQtlsGroup) <- qtlTypeNamesGroup
                                  }

                                  nGbs <- length(nQtlsGbs)
                                  if (all(nHaploQtlsGbs == 0)) {
                                    effectsHaploQtlsGbs <- rep(list(NULL), length(nHaploQtlsGbs))
                                    names(effectsHaploQtlsGbs) <- qtlTypeNamesGb
                                  } else {
                                    effectsHaploQtlsGbs <- sapply(X = 1:nGbs,
                                                                  FUN = function(gbNo) {
                                                                    nHaploQtlsGbNow <- nHaploQtlsGbs[gbNo]
                                                                    gbName <- qtlTypeNamesGb[gbNo]

                                                                    if (nHaploQtlsGbNow != 0) {
                                                                      if (effectsQtlsType == "same") {
                                                                        effectsHaploQtlsGbNow <- rep(1, nHaploQtlsGbNow)
                                                                      } else {
                                                                        effectsHaploQtlsGbNow <- rgamma(n = nHaploQtlsGbNow,
                                                                                                        shape = 2)
                                                                      }
                                                                      directionHaploQtlsGbNow <- sample(x = c(-1, 1), size = nHaploQtlsGbNow, replace = TRUE)
                                                                      effectsHaploQtlsGbNow <- effectsHaploQtlsGbNow * directionHaploQtlsGbNow
                                                                      names(effectsHaploQtlsGbNow) <- mrkNamesHaploQtlsList[[gbName]]
                                                                    } else {
                                                                      effectsHaploQtlsGbNow <- NULL
                                                                    }

                                                                    return(effectsHaploQtlsGbNow)
                                                                  })
                                    names(effectsHaploQtlsGbs) <- qtlTypeNamesGb
                                  }


                                  if (nHaploQtlsPg != 0) {
                                    if (effectsQtlsType == "same") {
                                      effectsHaploQtlsPg <- rep(1, nHaploQtlsPg)
                                    } else {
                                      effectsHaploQtlsPg <- rgamma(n = nHaploQtlsPg,
                                                                   shape = 2)
                                    }
                                    directionHaploQtlsPg <- sample(x = c(-1, 1), size = nHaploQtlsPg, replace = TRUE)
                                    effectsHaploQtlsPg <- effectsHaploQtlsPg * directionHaploQtlsPg

                                    names(effectsHaploQtlsPg) <- mrkNamesHaploQtlsList[["Polygenes"]]
                                  } else {
                                    effectsHaploQtlsPg <- NULL
                                  }




                                  ### 3.2.5.2. Determine normal QTL effects ###
                                  if (nQtlsCommon != 0) {
                                    if (effectsQtlsType == "same") {
                                      effectsQtlsCommonOri <- rep(1, nQtlsCommon)
                                    } else {
                                      effectsQtlsCommonOri <- rgamma(n = nQtlsCommon,
                                                                     shape = 2)
                                    }
                                    effectsQtlsCommon <- effectsQtlsCommonOri * rep(1, nQtlsCommon)
                                    directionQtlsCommon <- sample(x = c(-1, 1), size = nQtlsCommon, replace = TRUE)
                                    effectsQtlsCommon <- effectsQtlsCommon * directionQtlsCommon
                                    names(effectsQtlsCommon) <- mrkNamesQtlsList[["Common"]]
                                  } else {
                                    effectsQtlsCommon <- NULL
                                  }

                                  if (!is.null(effectsAncestryType)) {
                                    if (effectsQtlsType == "same") {
                                      effectsQtlsAncestryOri <- rep(1, sum(nQtlsAncestry))
                                    } else {
                                      effectsQtlsAncestryOri <- rgamma(n = sum(nQtlsAncestry),
                                                                       shape = 2)
                                    }

                                    effectsQtlsAncestryEach <- sqrt(ncol(effectsAncestryType)) / apply(X = effectsAncestryType, MARGIN = 1,
                                                                                                       FUN = function(x) sqrt(sum(x ^ 2)))
                                    effectsQtlsAncestryVec <- effectsQtlsAncestryOri *
                                      rep(effectsQtlsAncestryEach, nQtlsAncestry)
                                    directionQtlsAncestry <- sample(x = c(-1, 1), size = sum(nQtlsAncestry), replace = TRUE)
                                    effectsQtlsAncestryVec <- effectsQtlsAncestryVec * directionQtlsAncestry
                                    effectsQtlsAncestry <- split(x = effectsQtlsAncestryVec,
                                                                 f = rep(qtlTypeNamesAncestry, nQtlsAncestry))[qtlTypeNamesAncestry]
                                    effectsQtlsAncestry <- lapply(X = 1:length(qtlTypeNamesAncestry),
                                                                  FUN = function(qtlTypeNoAncestryEach) {
                                                                    effectsQtlsAncestryEach <- effectsQtlsAncestry[[qtlTypeNoAncestryEach]]
                                                                    if (!is.null(qtlTypeNoAncestryEach)) {
                                                                      names(effectsQtlsAncestryEach) <- (mrkNamesQtlsList[qtlTypeNamesAncestry])[[qtlTypeNoAncestryEach]]
                                                                    }

                                                                    return(effectsQtlsAncestryEach)
                                                                  })
                                    names(effectsQtlsAncestry) <- qtlTypeNamesAncestry

                                  } else {
                                    effectsQtlsAncestry <- rep(list(NULL), length(qtlTypeNamesAncestry))
                                    names(effectsQtlsAncestry) <- qtlTypeNamesAncestry
                                  }


                                  if (!is.null(effectsGroupType)) {
                                    if (effectsQtlsType == "same") {
                                      effectsQtlsGroupOri <- rep(1, sum(nQtlsGroup))
                                    } else {
                                      effectsQtlsGroupOri <- rgamma(n = sum(nQtlsGroup),
                                                                    shape = 2)
                                    }

                                    effectsQtlsGroupEach <- sqrt(ncol(effectsGroupType)) / apply(X = effectsGroupType, MARGIN = 1,
                                                                                                 FUN = function(x) sqrt(sum(x ^ 2)))
                                    effectsQtlsGroupVec <- effectsQtlsGroupOri * rep(effectsQtlsGroupEach, nQtlsGroup)
                                    directionQtlsGroup <- sample(x = c(-1, 1), size = sum(nQtlsGroup), replace = TRUE)
                                    effectsQtlsGroupVec <- effectsQtlsGroupVec * directionQtlsGroup
                                    effectsQtlsGroup <- split(x = effectsQtlsGroupVec,
                                                              f = rep(qtlTypeNamesGroup, nQtlsGroup))[qtlTypeNamesGroup]
                                    effectsQtlsGroup <- lapply(X = 1:length(qtlTypeNamesGroup),
                                                               FUN = function(qtlTypeNoGroupEach) {
                                                                 effectsQtlsGroupEach <- effectsQtlsGroup[[qtlTypeNoGroupEach]]
                                                                 if (!is.null(qtlTypeNoGroupEach)) {
                                                                   names(effectsQtlsGroupEach) <- (mrkNamesQtlsList[qtlTypeNamesGroup])[[qtlTypeNoGroupEach]]
                                                                 }

                                                                 return(effectsQtlsGroupEach)
                                                               })
                                    names(effectsQtlsGroup) <- qtlTypeNamesGroup
                                  } else {
                                    effectsQtlsGroup <- rep(list(NULL), length(qtlTypeNamesGroup))
                                    names(effectsQtlsGroup) <- qtlTypeNamesGroup
                                  }


                                  if (all(nQtlsGbs == 0)) {
                                    effectsQtlsGbs <- rep(list(NULL), nGbs)
                                    names(effectsQtlsGbs) <- qtlTypeNamesGb
                                  } else {
                                    effectsQtlsGbs <- sapply(X = 1:nGbs,
                                                             FUN = function(gbNo) {
                                                               nQtlsGbNow <- nQtlsGbs[gbNo]
                                                               gbName <- qtlTypeNamesGb[gbNo]

                                                               if (nQtlsGbNow != 0) {
                                                                 if (effectsQtlsType == "same") {
                                                                   effectsQtlsGbNow <- rep(1, nQtlsGbNow)
                                                                 } else {
                                                                   effectsQtlsGbNow <- rgamma(n = nQtlsGbNow,
                                                                                              shape = 2)
                                                                 }

                                                                 directionQtlsGbNow <- sample(x = c(-1, 1), size = nQtlsGbNow, replace = TRUE)
                                                                 effectsQtlsGbNow <- effectsQtlsGbNow * directionQtlsGbNow
                                                                 names(effectsQtlsGbNow) <- mrkNamesQtlsList[[gbName]]
                                                               } else {
                                                                 effectsQtlsGbNow <- NULL
                                                               }

                                                               return(effectsQtlsGbNow)
                                                             })
                                    names(effectsQtlsGbs) <- qtlTypeNamesGb
                                  }


                                  if (nQtlsPg != 0) {
                                    if (effectsQtlsType == "same") {
                                      effectsQtlsPg <- rep(1, nQtlsPg)
                                    } else {
                                      effectsQtlsPg <- rgamma(n = nQtlsPg,
                                                              shape = 2)
                                    }
                                    directionQtlsPg <- sample(x = c(-1, 1), size = nQtlsPg, replace = TRUE)
                                    effectsQtlsPg <- effectsQtlsPg * directionQtlsPg
                                    names(effectsQtlsPg) <- mrkNamesQtlsList[["Polygenes"]]
                                  } else {
                                    effectsQtlsPg <- NULL
                                  }








                                  #### 3.2.6. Determine polygenetic effects and compute genotypic values ####
                                  set.seed(seed = seeds["Polygenes_Effects"])


                                  ### 3.2.6.1. Determine common polygenetic effects and compute genotypic values ###
                                  if (nPolyGenesCommon != 0) {
                                    effectsPolyGenesCommon <- rnorm(n = nPolyGenesCommon,
                                                                    mean = 0, sd = 1)
                                    names(effectsPolyGenesCommon) <- mrkNamesPolyGenesList[["Common"]]
                                    genValsPolyGenesCommon <- genoScore[, mrkNamesPolyGenesList[["Common"]], drop = FALSE] %*%
                                      as.matrix(effectsPolyGenesCommon)
                                  } else {
                                    effectsPolyGenesCommon <- NULL
                                    genValsPolyGenesCommon <- matrix(data = 0,
                                                                     nrow = nLines,
                                                                     ncol = 1,
                                                                     dimnames = list(lineNames,
                                                                                     "Common"))
                                  }


                                  ### 3.2.6.2. Determine ancestral polygenetic effects and compute genotypic values ###
                                  if ((!is.null(effectsAncestryType)) & (sum(nPolyGenesAncestry) != 0)) {
                                    effectsPolyGenesAncestryVec <- rnorm(n = sum(nPolyGenesAncestry),
                                                                         mean = 0, sd = 1)
                                    effectsPolyGenesAncestry <- split(x = effectsPolyGenesAncestryVec,
                                                                      f = rep(polyGeneTypeNamesAncestry, nPolyGenesAncestry))[polyGeneTypeNamesAncestry]
                                    effectsPolyGenesAncestry <- lapply(X = 1:length(polyGeneTypeNamesAncestry),
                                                                       FUN = function(polyGeneTypeNoAncestryEach) {
                                                                         effectsPolyGenesAncestryEach <- effectsPolyGenesAncestry[[polyGeneTypeNoAncestryEach]]

                                                                         if (!is.null(effectsPolyGenesAncestryEach)) {
                                                                           effectsPolyGenesAncestryEach <- effectsPolyGenesAncestryEach * effectsQtlsAncestryEach[polyGeneTypeNoAncestryEach]
                                                                           names(effectsPolyGenesAncestryEach) <- (mrkNamesPolyGenesList[polyGeneTypeNamesAncestry])[[polyGeneTypeNoAncestryEach]]
                                                                         }

                                                                         return(effectsPolyGenesAncestryEach)
                                                                       })
                                    names(effectsPolyGenesAncestry) <- polyGeneTypeNamesAncestry

                                    genValsPolyGenesAncestry <-
                                      do.call(
                                        what = cbind,
                                        args = lapply(
                                          X = 1:length(polyGeneTypeNamesAncestry),
                                          FUN = function(polyGeneTypeNoAncestry) {
                                            mrkNamesPolyGenesAncestryEach <- mrkNamesPolyGenesList[[polyGeneTypeNamesAncestry[polyGeneTypeNoAncestry]]]

                                            if (!is.null(mrkNamesPolyGenesAncestryEach)) {
                                              genoScorePolyGenesAncestry <- genoScore[, mrkNamesPolyGenesAncestryEach, drop = FALSE]
                                              ancestralHaploDfPolyGenesAncestry <- ancestralHaploDf[, mrkNamesPolyGenesAncestryEach, drop = FALSE]
                                              effectsAncestryTypeEach <- effectsAncestryType[polyGeneTypeNoAncestry, ]


                                              ancestralHaploPolyGenesAncestry <- as.matrix(ancestralHaploDfPolyGenesAncestry)

                                              for (ancestralTypeNo in 1:length(effectsAncestryTypeEach)) {
                                                ancestralHaploPolyGenesAncestry[
                                                  ancestralHaploPolyGenesAncestry ==
                                                    names(effectsAncestryTypeEach)[ancestralTypeNo]
                                                ] <- effectsAncestryTypeEach[ancestralTypeNo]
                                              }
                                              ancestralHaploPolyGenesAncestry <- apply(X = ancestralHaploPolyGenesAncestry,
                                                                                       MARGIN = 2,
                                                                                       FUN = as.numeric)
                                              rownames(ancestralHaploPolyGenesAncestry) <-
                                                rownames(ancestralHaploDfPolyGenesAncestry)

                                              genoScorePolyGenesAncestryForGenVals <- genoScorePolyGenesAncestry * ancestralHaploPolyGenesAncestry

                                              genValsPolyGenesAncestryEach <- genoScorePolyGenesAncestryForGenVals %*%
                                                as.matrix(effectsPolyGenesAncestry[[polyGeneTypeNoAncestry]])
                                              colnames(genValsPolyGenesAncestryEach) <- polyGeneTypeNamesAncestry[polyGeneTypeNoAncestry]
                                            } else {
                                              genValsPolyGenesAncestryEach <- matrix(data = 0,
                                                                                     nrow = nLines,
                                                                                     ncol = 1,
                                                                                     dimnames = list(lineNames,
                                                                                                     polyGeneTypeNamesAncestry[polyGeneTypeNoAncestry]))
                                            }

                                            return(genValsPolyGenesAncestryEach)
                                          }
                                        )
                                      )


                                  } else {
                                    effectsPolyGenesAncestry <- rep(list(NULL), length(polyGeneTypeNamesAncestry))
                                    names(effectsPolyGenesAncestry) <- polyGeneTypeNamesAncestry

                                    genValsPolyGenesAncestry <- matrix(data = 0,
                                                                       nrow = nLines,
                                                                       ncol = length(polyGeneTypeNamesAncestry),
                                                                       dimnames = list(lineNames,
                                                                                       polyGeneTypeNamesAncestry))
                                  }



                                  ### 3.2.6.3. Determine group polygenetic effects and compute genotypic values ###
                                  if ((!is.null(effectsGroupType)) & (sum(nPolyGenesGroup) != 0)) {
                                    effectsPolyGenesGroupVec <- rnorm(n = sum(nPolyGenesGroup),
                                                                      mean = 0, sd = 1)
                                    effectsPolyGenesGroup <- split(x = effectsPolyGenesGroupVec,
                                                                   f = rep(polyGeneTypeNamesGroup, nPolyGenesGroup))[polyGeneTypeNamesGroup]
                                    effectsPolyGenesGroup <- lapply(X = 1:length(polyGeneTypeNamesGroup),
                                                                    FUN = function(polyGeneTypeNoGroupEach) {
                                                                      effectsPolyGenesGroupEach <- effectsPolyGenesGroup[[polyGeneTypeNoGroupEach]]
                                                                      if (!is.null(effectsPolyGenesGroupEach)) {
                                                                        effectsPolyGenesGroupEach <- effectsPolyGenesGroupEach * effectsQtlsGroupEach[polyGeneTypeNoGroupEach]
                                                                        names(effectsPolyGenesGroupEach) <- (mrkNamesPolyGenesList[polyGeneTypeNamesGroup])[[polyGeneTypeNoGroupEach]]
                                                                      }

                                                                      return(effectsPolyGenesGroupEach)
                                                                    })
                                    names(effectsPolyGenesGroup) <- polyGeneTypeNamesGroup

                                    genValsPolyGenesGroup <-
                                      do.call(
                                        what = cbind,
                                        args = lapply(
                                          X = 1:length(polyGeneTypeNamesGroup),
                                          FUN = function(polyGeneTypeNoGroup) {
                                            mrkNamesPolyGenesGroupEach <- mrkNamesPolyGenesList[[polyGeneTypeNamesGroup[polyGeneTypeNoGroup]]]

                                            if (!is.null(mrkNamesPolyGenesGroupEach)) {
                                              genoScorePolyGenesGroup <- genoScore[, mrkNamesPolyGenesGroupEach, drop = FALSE]
                                              gbInfoNumericEach <- rep(0, nLines)
                                              names(gbInfoNumericEach) <- gbInfo
                                              effectsGroupTypeEach <- effectsGroupType[polyGeneTypeNoGroup, ]

                                              for (groupTypeNo in 1:length(effectsGroupTypeEach)) {
                                                gbInfoNumericEach[gbInfo == names(effectsGroupTypeEach)[groupTypeNo]] <-
                                                  effectsGroupTypeEach[groupTypeNo]
                                              }


                                              genValsPolyGenesGroupEachRaw <- genoScorePolyGenesGroup %*%
                                                as.matrix(effectsPolyGenesGroup[[polyGeneTypeNoGroup]])
                                              colnames(genValsPolyGenesGroupEachRaw) <- polyGeneTypeNamesGroup[polyGeneTypeNoGroup]

                                              genValsPolyGenesGroupEach <- genValsPolyGenesGroupEachRaw *
                                                as.matrix(gbInfoNumericEach)
                                            } else {
                                              genValsPolyGenesGroupEach <- matrix(data = 0,
                                                                                  nrow = nLines,
                                                                                  ncol = 1,
                                                                                  dimnames = list(lineNames,
                                                                                                  polyGeneTypeNamesGroup[polyGeneTypeNoGroup]))
                                            }

                                            return(genValsPolyGenesGroupEach)
                                          }
                                        )
                                      )


                                  } else {
                                    effectsPolyGenesGroup <- rep(list(NULL), length(polyGeneTypeNamesGroup))
                                    names(effectsPolyGenesGroup) <- polyGeneTypeNamesGroup

                                    genValsPolyGenesGroup <- matrix(data = 0,
                                                                    nrow = nLines,
                                                                    ncol = length(polyGeneTypeNamesGroup),
                                                                    dimnames = list(lineNames,
                                                                                    polyGeneTypeNamesGroup))
                                  }





                                  ### 3.2.6.4. Summary genotypic values for polygenetic effects ###
                                  genValsPolyGenesRawMat <- cbind(genValsPolyGenesCommon,
                                                                  genValsPolyGenesAncestry,
                                                                  genValsPolyGenesGroup)
                                  rownames(genValsPolyGenesRawMat) <- lineNames
                                  colnames(genValsPolyGenesRawMat) <- polyGeneTypeNamesAll
                                  genValsPolyGenesRaw <- apply(X = genValsPolyGenesRawMat,
                                                               MARGIN = 1, FUN = sum)





                                  #### 3.2.7. Compute genotypic values for QTL effects ####
                                  whichGenoHasHaploQtlEffects <- function(mrkNameHaploQtls) {
                                    blockNameIncHaploQTls <- haploBlockDf$block[haploBlockDf$marker == mrkNameHaploQtls]
                                    haploBlockIncHaploQtls <- haploBlockDf[haploBlockDf$block %in% blockNameIncHaploQTls, ]

                                    genoScoreIncHaploQtls <- genoScore[, haploBlockIncHaploQtls$marker]

                                    freqHaploQtls <- mean(x = genoScoreIncHaploQtls[, mrkNameHaploQtls]) / 2
                                    minAlleleHaploQtls <- ifelse(test = freqHaploQtls >= 0.5,
                                                                 yes = 0, no = 1)

                                    whichHasMinAlleleHaploQtls <- genoScoreIncHaploQtls[, mrkNameHaploQtls] %in% c(1, minAlleleHaploQtls * 2)
                                    decreasingInOrder <- minAlleleHaploQtls == 0


                                    clustRes <- kmeans(x = genoScoreIncHaploQtls[whichHasMinAlleleHaploQtls, ],
                                                       centers = 2,
                                                       nstart = 100)
                                    scoreHaploCenterIncHaploQtls <- apply(X = clustRes$centers,
                                                                          MARGIN = 1,
                                                                          FUN = sum)
                                    scoreHaploCenterIncHaploQtlsSorted <- sort(scoreHaploCenterIncHaploQtls,
                                                                               decreasing = !decreasingInOrder)
                                    clusterNoWithQtlEffects <- as.numeric(names(scoreHaploCenterIncHaploQtlsSorted)[1])
                                    indNamesWithQtlsEffects <- names(clustRes$cluster)[clustRes$cluster == clusterNoWithQtlEffects]

                                    whichGenoHasHaploQtlsEffects <- rownames(genoScoreIncHaploQtls) %in% indNamesWithQtlsEffects

                                    return(whichGenoHasHaploQtlsEffects)
                                  }

                                  convertGenoScoreForGenVals <- function(genoScoreNow) {
                                    freqHaploQtls <- mean(x = genoScoreNow) / 2
                                    minAlleleHaploQtls <- ifelse(test = freqHaploQtls >= 0.5,
                                                                 yes = 0, no = 1)
                                    if (minAlleleHaploQtls == 0) {
                                      genoScoreNow <- 2 - genoScoreNow
                                    }

                                    return(genoScoreNow)
                                  }


                                  ### 3.2.7.1. Compute genotypic values for common QTL effects ###
                                  if (nHaploQtlsCommon != 0) {
                                    whichGenoHasHaploQtlsCommonEffectsMat <- sapply(
                                      X = mrkNamesHaploQtlsList[["Common"]],
                                      FUN = whichGenoHasHaploQtlEffects, simplify = TRUE)
                                    rownames(whichGenoHasHaploQtlsCommonEffectsMat) <- lineNames


                                    genoScoreHaploQtlsCommon <- genoScore[, mrkNamesHaploQtlsList[["Common"]],
                                                                          drop = FALSE]
                                    sdHaploQtlsCommon <- apply(X = genoScoreHaploQtlsCommon,
                                                               MARGIN = 2, FUN = sd)
                                    effectsHaploQtlsCommonReal <- effectsHaploQtlsCommon / sdHaploQtlsCommon

                                    genoScoreHaploQtlsCommonForGenVals <- apply(
                                      X = genoScoreHaploQtlsCommon,
                                      MARGIN = 2,
                                      FUN = convertGenoScoreForGenVals
                                    )
                                    genoScoreHaploQtlsCommonForGenVals <- genoScoreHaploQtlsCommonForGenVals *
                                      whichGenoHasHaploQtlsCommonEffectsMat
                                    genValsHaploQtlsCommon <- genoScoreHaploQtlsCommonForGenVals %*%
                                      as.matrix(effectsHaploQtlsCommonReal)
                                  } else {
                                    effectsHaploQtlsCommonReal <- NULL
                                    genValsHaploQtlsCommon <- matrix(data = 0,
                                                                     nrow = nLines,
                                                                     ncol = 1,
                                                                     dimnames = list(lineNames,
                                                                                     "Common"))
                                  }



                                  if (nQtlsCommon != 0) {
                                    genoScoreQtlsCommon <- genoScore[, mrkNamesQtlsList[["Common"]],
                                                                     drop = FALSE]
                                    sdQtlsCommon <- apply(X = genoScoreQtlsCommon,
                                                          MARGIN = 2, FUN = sd)
                                    effectsQtlsCommonReal <- effectsQtlsCommon / sdQtlsCommon
                                    genoScoreQtlsCommonForGenVals <- apply(
                                      X = genoScoreQtlsCommon,
                                      MARGIN = 2,
                                      FUN = convertGenoScoreForGenVals
                                    )
                                    genValsQtlsCommon <- genoScoreQtlsCommonForGenVals %*%
                                      as.matrix(effectsQtlsCommonReal)
                                  } else {
                                    effectsQtlsCommonReal <- NULL
                                    genValsQtlsCommon <- matrix(data = 0,
                                                                nrow = nLines,
                                                                ncol = 1,
                                                                dimnames = list(lineNames,
                                                                                "Common"))
                                  }



                                  ### 3.2.7.2. Compute genotypic values for ancestral QTL effects ###
                                  if ((!is.null(effectsAncestryType)) & (sum(nHaploQtlsAncestry) != 0)) {
                                    genValsHaploQtlsAncestryRes <-
                                      lapply(
                                        X = 1:length(qtlTypeNamesAncestry),
                                        FUN = function(qtlTypeNoAncestry) {
                                          mrkNamesHaploQtlsAncestryEach <- mrkNamesHaploQtlsList[[qtlTypeNamesAncestry[qtlTypeNoAncestry]]]

                                          if (!is.null(mrkNamesHaploQtlsAncestryEach)) {
                                            whichGenoHasHaploQtlsAncestryEffectsMat <- sapply(
                                              X = mrkNamesHaploQtlsAncestryEach,
                                              FUN = whichGenoHasHaploQtlEffects, simplify = TRUE)
                                            rownames(whichGenoHasHaploQtlsAncestryEffectsMat) <- lineNames

                                            genoScoreHaploQtlsAncestry <- genoScore[, mrkNamesHaploQtlsAncestryEach, drop = FALSE]
                                            ancestralHaploDfHaploQtlsAncestry <- ancestralHaploDf[, mrkNamesHaploQtlsAncestryEach, drop = FALSE]
                                            effectsAncestryTypeEach <- effectsAncestryType[qtlTypeNoAncestry, ]

                                            ancestralHaploHaploQtlsAncestry <- as.matrix(ancestralHaploDfHaploQtlsAncestry)

                                            for (ancestralTypeNo in 1:length(effectsAncestryTypeEach)) {
                                              ancestralHaploHaploQtlsAncestry[
                                                ancestralHaploHaploQtlsAncestry ==
                                                  names(effectsAncestryTypeEach)[ancestralTypeNo]
                                              ] <- effectsAncestryTypeEach[ancestralTypeNo]
                                            }
                                            ancestralHaploHaploQtlsAncestry <- apply(X = ancestralHaploHaploQtlsAncestry,
                                                                                     MARGIN = 2,
                                                                                     FUN = as.numeric)
                                            rownames(ancestralHaploHaploQtlsAncestry) <-
                                              rownames(ancestralHaploDfHaploQtlsAncestry)

                                            genoScoreHaploQtlsAncestryForGenVals <- apply(
                                              X = genoScoreHaploQtlsAncestry,
                                              MARGIN = 2,
                                              FUN = convertGenoScoreForGenVals
                                            )
                                            genoScoreHaploQtlsAncestryForGenVals <- genoScoreHaploQtlsAncestryForGenVals *
                                              ancestralHaploHaploQtlsAncestry *
                                              whichGenoHasHaploQtlsAncestryEffectsMat

                                            sdHaploQtlsAncestryEach <- apply(X = genoScoreHaploQtlsAncestry,
                                                                             MARGIN = 2,
                                                                             FUN = sd)
                                            effectsHaploQtlsAncestryEachReal <- effectsHaploQtlsAncestry[[qtlTypeNoAncestry]] / sdHaploQtlsAncestryEach
                                            genValsHaploQtlsAncestryEach <- genoScoreHaploQtlsAncestryForGenVals %*%
                                              as.matrix(effectsHaploQtlsAncestryEachReal)
                                            colnames(genValsHaploQtlsAncestryEach) <- qtlTypeNamesAncestry[qtlTypeNoAncestry]
                                          } else {
                                            effectsHaploQtlsAncestryEachReal <- NULL
                                            genValsHaploQtlsAncestryEach <- matrix(data = 0,
                                                                                   nrow = nLines,
                                                                                   ncol = 1,
                                                                                   dimnames = list(lineNames,
                                                                                                   qtlTypeNamesAncestry[qtlTypeNoAncestry]))
                                          }

                                          return(list(effectsHaploQtlsAncestryEachReal = effectsHaploQtlsAncestryEachReal,
                                                      genValsHaploQtlsAncestryEach = genValsHaploQtlsAncestryEach))
                                        }
                                      )





                                    effectsHaploQtlsAncestryReal <- lapply(X = genValsHaploQtlsAncestryRes,
                                                                           FUN = function(x) x$effectsHaploQtlsAncestryEachReal)
                                    names(effectsHaploQtlsAncestryReal) <- qtlTypeNamesAncestry

                                    genValsHaploQtlsAncestry <- do.call(
                                      what = cbind,
                                      args = lapply(X = genValsHaploQtlsAncestryRes,
                                                    FUN = function(x) x$genValsHaploQtlsAncestryEach)
                                    )
                                  } else {
                                    effectsHaploQtlsAncestryReal <- rep(list(NULL), length(qtlTypeNamesAncestry))
                                    names(effectsHaploQtlsAncestryReal) <- qtlTypeNamesAncestry
                                    genValsHaploQtlsAncestry <- matrix(data = 0,
                                                                       nrow = nLines,
                                                                       ncol = length(qtlTypeNamesAncestry),
                                                                       dimnames = list(lineNames,
                                                                                       qtlTypeNamesAncestry))
                                  }





                                  if ((!is.null(effectsAncestryType)) & (sum(nQtlsAncestry) != 0)) {
                                    genValsQtlsAncestryRes <- lapply(
                                      X = 1:length(qtlTypeNamesAncestry),
                                      FUN = function(qtlTypeNoAncestry) {
                                        mrkNamesQtlsAncestryEach <- mrkNamesQtlsList[[qtlTypeNamesAncestry[qtlTypeNoAncestry]]]

                                        if (!is.null(mrkNamesQtlsAncestryEach)) {
                                          genoScoreQtlsAncestry <- genoScore[, mrkNamesQtlsAncestryEach, drop = FALSE]
                                          ancestralHaploDfQtlsAncestry <- ancestralHaploDf[, mrkNamesQtlsAncestryEach, drop = FALSE]
                                          effectsAncestryTypeEach <- effectsAncestryType[qtlTypeNoAncestry, ]


                                          ancestralHaploQtlsAncestry <- as.matrix(ancestralHaploDfQtlsAncestry)

                                          for (ancestralTypeNo in 1:length(effectsAncestryTypeEach)) {
                                            ancestralHaploQtlsAncestry[
                                              ancestralHaploQtlsAncestry ==
                                                names(effectsAncestryTypeEach)[ancestralTypeNo]
                                            ] <- effectsAncestryTypeEach[ancestralTypeNo]
                                          }
                                          ancestralHaploQtlsAncestry <- apply(X = ancestralHaploQtlsAncestry,
                                                                              MARGIN = 2,
                                                                              FUN = as.numeric)
                                          rownames(ancestralHaploQtlsAncestry) <-
                                            rownames(ancestralHaploDfQtlsAncestry)

                                          genoScoreQtlsAncestryForGenVals <- apply(
                                            X = genoScoreQtlsAncestry,
                                            MARGIN = 2,
                                            FUN = convertGenoScoreForGenVals
                                          )
                                          genoScoreQtlsAncestryForGenVals <- genoScoreQtlsAncestryForGenVals *
                                            ancestralHaploQtlsAncestry

                                          sdQtlsAncestryEach <- apply(X = genoScoreQtlsAncestry,
                                                                      MARGIN = 2,
                                                                      FUN = sd)
                                          effectsQtlsAncestryEachReal <- effectsQtlsAncestry[[qtlTypeNoAncestry]] / sdQtlsAncestryEach
                                          genValsQtlsAncestryEach <- genoScoreQtlsAncestryForGenVals %*%
                                            as.matrix(effectsQtlsAncestryEachReal)
                                          colnames(genValsQtlsAncestryEach) <- qtlTypeNamesAncestry[qtlTypeNoAncestry]
                                        } else {
                                          effectsQtlsAncestryEachReal <- NULL
                                          genValsQtlsAncestryEach <- matrix(data = 0,
                                                                            nrow = nLines,
                                                                            ncol = 1,
                                                                            dimnames = list(lineNames,
                                                                                            qtlTypeNamesAncestry[qtlTypeNoAncestry]))
                                        }

                                        return(list(effectsQtlsAncestryEachReal = effectsQtlsAncestryEachReal,
                                                    genValsQtlsAncestryEach = genValsQtlsAncestryEach))
                                      }
                                    )

                                    effectsQtlsAncestryReal <- lapply(X = genValsQtlsAncestryRes,
                                                                      FUN = function(x) x$effectsQtlsAncestryEachReal)
                                    names(effectsQtlsAncestryReal) <- qtlTypeNamesAncestry

                                    genValsQtlsAncestry <- do.call(
                                      what = cbind,
                                      args = lapply(X = genValsQtlsAncestryRes,
                                                    FUN = function(x) x$genValsQtlsAncestryEach)
                                    )
                                  } else {
                                    effectsQtlsAncestryReal <- rep(list(NULL), length(qtlTypeNamesAncestry))
                                    names(effectsQtlsAncestryReal) <- qtlTypeNamesAncestry
                                    genValsQtlsAncestry <- matrix(data = 0,
                                                                  nrow = nLines,
                                                                  ncol = length(qtlTypeNamesAncestry),
                                                                  dimnames = list(lineNames,
                                                                                  qtlTypeNamesAncestry))
                                  }





                                  ### 3.2.7.3. Compute genotypic values for group QTL effects ###
                                  if ((!is.null(effectsGroupType)) & (sum(nHaploQtlsGroup) != 0)) {
                                    genValsHaploQtlsGroupRes <-
                                      lapply(
                                        X = 1:length(qtlTypeNamesGroup),
                                        FUN = function(qtlTypeNoGroup) {
                                          mrkNamesHaploQtlsGroupEach <- mrkNamesHaploQtlsList[[qtlTypeNamesGroup[qtlTypeNoGroup]]]

                                          if (!is.null(mrkNamesHaploQtlsGroupEach)) {
                                            whichGenoHasHaploQtlsGroupEffectsMat <- sapply(
                                              X = mrkNamesHaploQtlsGroupEach,
                                              FUN = whichGenoHasHaploQtlEffects, simplify = TRUE)
                                            rownames(whichGenoHasHaploQtlsGroupEffectsMat) <- lineNames


                                            genoScoreHaploQtlsGroup <- genoScore[, mrkNamesHaploQtlsGroupEach, drop = FALSE]
                                            gbInfoNumericEach <- rep(0, nLines)
                                            names(gbInfoNumericEach) <- gbInfo
                                            effectsGroupTypeEach <- effectsGroupType[qtlTypeNoGroup, ]

                                            for (groupTypeNo in 1:length(effectsGroupTypeEach)) {
                                              gbInfoNumericEach[gbInfo == names(effectsGroupTypeEach)[groupTypeNo]] <-
                                                effectsGroupTypeEach[groupTypeNo]
                                            }

                                            genoScoreHaploQtlsGroupForGenVals <- apply(
                                              X = genoScoreHaploQtlsGroup,
                                              MARGIN = 2,
                                              FUN = convertGenoScoreForGenVals
                                            )
                                            genoScoreHaploQtlsGroupForGenVals <- genoScoreHaploQtlsGroupForGenVals *
                                              whichGenoHasHaploQtlsGroupEffectsMat

                                            sdHaploQtlsGroupEach <- apply(X = genoScoreHaploQtlsGroup,
                                                                          MARGIN = 2,
                                                                          FUN = sd)
                                            effectsHaploQtlsGroupEachReal <- effectsHaploQtlsGroup[[qtlTypeNoGroup]] / sdHaploQtlsGroupEach

                                            genValsHaploQtlsGroupEachRaw <- genoScoreHaploQtlsGroupForGenVals %*%
                                              as.matrix(effectsHaploQtlsGroupEachReal)
                                            colnames(genValsHaploQtlsGroupEachRaw) <- qtlTypeNamesGroup[qtlTypeNoGroup]

                                            genValsHaploQtlsGroupEach <- genValsHaploQtlsGroupEachRaw * as.matrix(gbInfoNumericEach)
                                          } else {
                                            effectsHaploQtlsGroupEachReal <- NULL
                                            genValsHaploQtlsGroupEach <- matrix(data = 0,
                                                                                nrow = nLines,
                                                                                ncol = 1,
                                                                                dimnames = list(lineNames,
                                                                                                qtlTypeNamesGroup[qtlTypeNoGroup]))
                                          }

                                          return(list(effectsHaploQtlsGroupEachReal = effectsHaploQtlsGroupEachReal,
                                                      genValsHaploQtlsGroupEach = genValsHaploQtlsGroupEach))
                                        }
                                      )

                                    effectsHaploQtlsGroupReal <- lapply(X = genValsHaploQtlsGroupRes,
                                                                        FUN = function(x) x$effectsHaploQtlsGroupEachReal)
                                    names(effectsHaploQtlsGroupReal) <- qtlTypeNamesGroup

                                    genValsHaploQtlsGroup <- do.call(
                                      what = cbind,
                                      args = lapply(X = genValsHaploQtlsGroupRes,
                                                    FUN = function(x) x$genValsHaploQtlsGroupEach)
                                    )



                                  } else {
                                    effectsHaploQtlsGroupReal <- rep(list(NULL), length(qtlTypeNamesGroup))
                                    names(effectsHaploQtlsGroupReal) <- qtlTypeNamesGroup
                                    genValsHaploQtlsGroup <- matrix(data = 0,
                                                                    nrow = nLines,
                                                                    ncol = length(qtlTypeNamesGroup),
                                                                    dimnames = list(lineNames,
                                                                                    qtlTypeNamesGroup))
                                  }






                                  if ((!is.null(effectsGroupType)) & (sum(nQtlsGroup) != 0)) {
                                    genValsQtlsGroupRes <-
                                      lapply(
                                        X = 1:length(qtlTypeNamesGroup),
                                        FUN = function(qtlTypeNoGroup) {
                                          mrkNamesQtlsGroupEach <- mrkNamesQtlsList[[qtlTypeNamesGroup[qtlTypeNoGroup]]]

                                          if (!is.null(mrkNamesQtlsGroupEach)) {
                                            genoScoreQtlsGroup <- genoScore[, mrkNamesQtlsGroupEach, drop = FALSE]
                                            gbInfoNumericEach <- rep(0, nLines)
                                            names(gbInfoNumericEach) <- gbInfo
                                            effectsGroupTypeEach <- effectsGroupType[qtlTypeNoGroup, ]

                                            for (groupTypeNo in 1:length(effectsGroupTypeEach)) {
                                              gbInfoNumericEach[gbInfo == names(effectsGroupTypeEach)[groupTypeNo]] <-
                                                effectsGroupTypeEach[groupTypeNo]
                                            }

                                            genoScoreQtlsGroupForGenVals <- apply(
                                              X = genoScoreQtlsGroup,
                                              MARGIN = 2,
                                              FUN = convertGenoScoreForGenVals
                                            )

                                            sdQtlsGroupEach <- apply(X = genoScoreQtlsGroup,
                                                                     MARGIN = 2,
                                                                     FUN = sd)
                                            effectsQtlsGroupEachReal <- effectsQtlsGroup[[qtlTypeNoGroup]] / sdQtlsGroupEach

                                            genValsQtlsGroupEachRaw <- genoScoreQtlsGroupForGenVals %*%
                                              as.matrix(effectsQtlsGroupEachReal)
                                            colnames(genValsQtlsGroupEachRaw) <- qtlTypeNamesGroup[qtlTypeNoGroup]

                                            genValsQtlsGroupEach <- genValsQtlsGroupEachRaw * as.matrix(gbInfoNumericEach)
                                          } else {
                                            effectsQtlsGroupEachReal <- NULL
                                            genValsQtlsGroupEach <- matrix(data = 0,
                                                                           nrow = nLines,
                                                                           ncol = 1,
                                                                           dimnames = list(lineNames,
                                                                                           qtlTypeNamesGroup[qtlTypeNoGroup]))
                                          }

                                          return(list(effectsQtlsGroupEachReal = effectsQtlsGroupEachReal,
                                                      genValsQtlsGroupEach = genValsQtlsGroupEach))
                                        }
                                      )




                                    effectsQtlsGroupReal <- lapply(X = genValsQtlsGroupRes,
                                                                   FUN = function(x) x$effectsQtlsGroupEachReal)
                                    names(effectsQtlsGroupReal) <- qtlTypeNamesGroup

                                    genValsQtlsGroup <- do.call(
                                      what = cbind,
                                      args = lapply(X = genValsQtlsGroupRes,
                                                    FUN = function(x) x$genValsQtlsGroupEach)
                                    )


                                  } else {
                                    effectsQtlsGroupReal <- rep(list(NULL), length(qtlTypeNamesGroup))
                                    names(effectsQtlsGroupReal) <- qtlTypeNamesGroup
                                    genValsQtlsGroup <- matrix(data = 0,
                                                               nrow = nLines,
                                                               ncol = length(qtlTypeNamesGroup),
                                                               dimnames = list(lineNames,
                                                                               qtlTypeNamesGroup))
                                  }




                                  ### 3.2.7.4. Compute genotypic values for QTL effects depending on genetic backgrounds (PCs of kinship) ###
                                  if (all(nHaploQtlsGbs == 0)) {
                                    effectsHaploQtlsGbsReal <- list()
                                    genValsHaploQtlsGbs <- matrix(data = 0,
                                                                  nrow = nLines,
                                                                  ncol = nGbs,
                                                                  dimnames = list(lineNames,
                                                                                  qtlTypeNamesGb))
                                  } else {
                                    genValsHaploQtlsGbsCompList <-
                                      sapply(X = 1:nGbs,
                                             FUN = function(gbNo) {
                                               nHaploQtlsGbNow <- nHaploQtlsGbs[gbNo]
                                               gbName <- qtlTypeNamesGb[gbNo]
                                               if (nHaploQtlsGbNow != 0) {
                                                 mrkNamesHaploQtlsGbNow <- mrkNamesHaploQtlsList[[gbName]]
                                                 effectsHaploQtlsGbNow <- effectsHaploQtlsGbs[[gbName]]
                                                 gbsMatNow <- gbsMat[, rep(gbNo, nHaploQtlsGbNow), drop = FALSE]
                                                 gbsMatNow <- scale(gbsMatNow)
                                                 genoScoreHaploQtlsGbNow <- genoScore[, mrkNamesHaploQtlsGbNow, drop = FALSE]

                                                 whichGenoHasHaploQtlsGbNowEffectsMat <- sapply(
                                                   X = mrkNamesHaploQtlsGbNow,
                                                   FUN = whichGenoHasHaploQtlEffects, simplify = TRUE
                                                 )
                                                 genoScoreHaploQtlsGbNowForGenVals <- apply(
                                                   X = genoScoreHaploQtlsGbNow,
                                                   MARGIN = 2,
                                                   FUN = convertGenoScoreForGenVals
                                                 )
                                                 genoScoreHaploQtlsGbNowForGenVals <- genoScoreHaploQtlsGbNowForGenVals *
                                                   whichGenoHasHaploQtlsGbNowEffectsMat
                                                 designMatGbNow <- genoScoreHaploQtlsGbNowForGenVals * gbsMatNow
                                                 sdHaploQtlsGbNow <- apply(X = designMatGbNow,
                                                                           MARGIN = 2, FUN = sd)
                                                 effectsHaploQtlsGbRealNow <- effectsHaploQtlsGbNow / sdHaploQtlsGbNow
                                                 genValsHaploQtlsGbNow <- designMatGbNow %*%
                                                   as.matrix(effectsHaploQtlsGbRealNow)
                                                 colnames(genValsHaploQtlsGbNow) <- gbName
                                               } else {
                                                 effectsHaploQtlsGbRealNow <- NULL
                                                 genValsHaploQtlsGbNow <- matrix(data = 0,
                                                                                 nrow = nLines,
                                                                                 ncol = 1,
                                                                                 dimnames = list(lineNames,
                                                                                                 gbName))
                                               }

                                               return(list(effectsHaploQtlsGbReal = effectsHaploQtlsGbRealNow,
                                                           genValsHaploQtlsGb = genValsHaploQtlsGbNow))
                                             }, simplify = FALSE)
                                    names(genValsHaploQtlsGbsCompList) <- qtlTypeNamesGb


                                    effectsHaploQtlsGbsReal <- lapply(X = genValsHaploQtlsGbsCompList,
                                                                      FUN = function(x) x$effectsHaploQtlsGbReal)

                                    genValsHaploQtlsGbs <-
                                      do.call(
                                        what = cbind,
                                        args = lapply(X = genValsHaploQtlsGbsCompList,
                                                      FUN = function(x) x$genValsHaploQtlsGb)
                                      )
                                  }




                                  if (all(nQtlsGbs == 0)) {
                                    effectsQtlsGbsReal <- list()
                                    genValsQtlsGbs <- matrix(data = 0,
                                                             nrow = nLines,
                                                             ncol = nGbs,
                                                             dimnames = list(lineNames,
                                                                             qtlTypeNamesGb))
                                  } else {
                                    genValsQtlsGbsCompList <-
                                      sapply(X = 1:nGbs,
                                             FUN = function(gbNo) {
                                               nQtlsGbNow <- nQtlsGbs[gbNo]
                                               gbName <- qtlTypeNamesGb[gbNo]
                                               if (nQtlsGbNow != 0) {
                                                 mrkNamesQtlsGbNow <- mrkNamesQtlsGbs[[gbName]]
                                                 effectsQtlsGbNow <- effectsQtlsGbs[[gbName]]
                                                 gbsMatNow <- gbsMat[, rep(gbNo, nQtlsGbNow), drop = FALSE]
                                                 gbsMatNow <- scale(gbsMatNow)
                                                 genoScoreQtlsGbNow <- genoScore[, mrkNamesQtlsGbNow, drop = FALSE]

                                                 genoScoreQtlsGbNowForGenVals <- apply(
                                                   X = genoScoreQtlsGbNow,
                                                   MARGIN = 2,
                                                   FUN = convertGenoScoreForGenVals
                                                 )
                                                 designMatGbNow <- genoScoreQtlsGbNowForGenVals * gbsMatNow
                                                 sdQtlsGbNow <- apply(X = designMatGbNow,
                                                                      MARGIN = 2, FUN = sd)
                                                 effectsQtlsGbRealNow <- effectsQtlsGbNow / sdQtlsGbNow
                                                 genValsQtlsGbNow <- designMatGbNow %*%
                                                   as.matrix(effectsQtlsGbRealNow)
                                                 colnames(genValsQtlsGbNow) <- gbName
                                               } else {
                                                 effectsQtlsGbRealNow <- NULL
                                                 genValsQtlsGbNow <- matrix(data = 0,
                                                                            nrow = nLines,
                                                                            ncol = 1,
                                                                            dimnames = list(lineNames,
                                                                                            gbName))
                                               }

                                               return(list(effectsQtlsGbReal = effectsQtlsGbRealNow,
                                                           genValsQtlsGb = genValsQtlsGbNow))
                                             }, simplify = FALSE)
                                    names(genValsQtlsGbsCompList) <- qtlTypeNamesGb


                                    effectsQtlsGbsReal <- lapply(X = genValsQtlsGbsCompList,
                                                                 FUN = function(x) x$effectsQtlsGbReal)

                                    genValsQtlsGbs <-
                                      do.call(
                                        what = cbind,
                                        args = lapply(X = genValsQtlsGbsCompList,
                                                      FUN = function(x) x$genValsQtlsGb)
                                      )
                                  }




                                  ### 3.2.7.5. Compute genotypic values for QTL effects depending on polygenes ###
                                  if (nHaploQtlsPg != 0) {
                                    genoScoreHaploQtlsPg <- genoScore[, mrkNamesHaploQtlsList[[qtlTypeNamesPg]], drop = FALSE]
                                    genValsPolyGenesRawForInteraction <- as.matrix(x = genValsPolyGenesRaw)[, rep(1, nHaploQtlsPg), drop = FALSE]
                                    genValsPolyGenesRawForInteraction <- scale(genValsPolyGenesRawForInteraction)

                                    genoScoreHaploQtlsPgForGenVals <- apply(
                                      X = genoScoreHaploQtlsPg,
                                      MARGIN = 2,
                                      FUN = convertGenoScoreForGenVals
                                    )
                                    designMatPg <- genoScoreHaploQtlsPgForGenVals * genValsPolyGenesRawForInteraction
                                    sdHaploQtlsPg <- apply(X = designMatPg,
                                                           MARGIN = 2, FUN = sd)
                                    effectsHaploQtlsPgReal <- effectsHaploQtlsPg / sdHaploQtlsPg
                                    genValsHaploQtlsPg <- designMatPg %*%
                                      as.matrix(effectsHaploQtlsPgReal)
                                    colnames(genValsHaploQtlsPg) <- qtlTypeNamesPg
                                  } else {
                                    effectsHaploQtlsPgReal <- NULL
                                    genValsHaploQtlsPg <- matrix(data = 0,
                                                                 nrow = nLines,
                                                                 ncol = 1,
                                                                 dimnames = list(lineNames,
                                                                                 qtlTypeNamesPg))
                                  }



                                  if (nQtlsPg != 0) {
                                    genoScoreQtlsPg <- genoScore[, mrkNamesQtlsList[[qtlTypeNamesPg]], drop = FALSE]
                                    genValsPolyGenesRawForInteraction <- as.matrix(x = genValsPolyGenesRaw)[, rep(1, nQtlsPg), drop = FALSE]
                                    genValsPolyGenesRawForInteraction <- scale(genValsPolyGenesRawForInteraction)

                                    genoScoreQtlsPgForGenVals <- apply(
                                      X = genoScoreQtlsPg,
                                      MARGIN = 2,
                                      FUN = convertGenoScoreForGenVals
                                    )
                                    designMatPg <- genoScoreQtlsPgForGenVals * genValsPolyGenesRawForInteraction
                                    sdQtlsPg <- apply(X = designMatPg,
                                                      MARGIN = 2, FUN = sd)
                                    effectsQtlsPgReal <- effectsQtlsPg / sdQtlsPg
                                    genValsQtlsPg <- designMatPg %*%
                                      as.matrix(effectsQtlsPgReal)
                                    colnames(genValsQtlsPg) <- qtlTypeNamesPg
                                  } else {
                                    effectsQtlsPgReal <- NULL
                                    genValsQtlsPg <- matrix(data = 0,
                                                            nrow = nLines,
                                                            ncol = 1,
                                                            dimnames = list(lineNames,
                                                                            qtlTypeNamesPg))
                                  }





                                  ### 3.2.7.6. Summary genotypic values for QTL and haplotype-based QTL effects ###
                                  genValsHaploQtlsMat <- cbind(genValsHaploQtlsCommon,
                                                               genValsHaploQtlsAncestry,
                                                               genValsHaploQtlsGroup,
                                                               genValsHaploQtlsGbs,
                                                               genValsHaploQtlsPg)
                                  colnames(genValsHaploQtlsMat) <- qtlTypeNamesAll
                                  genValsHaploQtls <- apply(X = genValsHaploQtlsMat,
                                                            MARGIN = 1, FUN = sum)



                                  genValsQtlsMat <- cbind(genValsQtlsCommon,
                                                          genValsQtlsAncestry,
                                                          genValsQtlsGroup,
                                                          genValsQtlsGbs,
                                                          genValsQtlsPg)
                                  colnames(genValsQtlsMat) <- qtlTypeNamesAll
                                  genValsQtls <- apply(X = genValsQtlsMat,
                                                       MARGIN = 1, FUN = sum)


                                  genValsQtlsAllMat <- cbind(genValsHaploQtlsMat,
                                                             genValsQtlsMat)
                                  colnames(genValsQtlsAllMat) <-
                                    c(paste0("HaploQtl-", colnames(genValsHaploQtlsMat)),
                                      paste0("Qtl-", colnames(genValsQtlsMat)))
                                  genValsQtlsAll <- genValsHaploQtls + genValsQtls



                                  #### 3.2.8. Define true genotypic values from variance components and compute total genotypic values ####
                                  heritPolyGenes <- herit - heritQtls
                                  varPropPolyGenes <- heritPolyGenes / heritQtls
                                  correctionForGenValsPolyGenes <- sqrt(varPropPolyGenes) *
                                    (sd(genValsQtlsAll) / sd(genValsPolyGenesRaw))
                                  genValsPolyGenesMat <- genValsPolyGenesRawMat * correctionForGenValsPolyGenes

                                  genValsPolyGenes <- apply(X = genValsPolyGenesMat,
                                                            MARGIN = 1, FUN = sum)
                                  genVals <- genValsPolyGenes + genValsQtlsAll




                                  #### 3.2.9. Summary effects and genotypic values ####
                                  effectsList <- list(
                                    Qtls = list(
                                      HaploQtls = list(
                                        Param = c(list(
                                          Common = effectsHaploQtlsCommon
                                        ),
                                        effectsHaploQtlsAncestry,
                                        effectsHaploQtlsGroup,
                                        effectsHaploQtlsGbs,
                                        list(Polygenes = effectsHaploQtlsPg)),
                                        True = c(list(
                                          Common = effectsHaploQtlsCommonReal
                                        ),
                                        effectsHaploQtlsAncestryReal,
                                        effectsHaploQtlsGroupReal,
                                        effectsHaploQtlsGbsReal,
                                        list(Polygenes = effectsHaploQtlsPgReal))
                                      ),
                                      Qtls = list(
                                        Param = c(list(
                                          Common = effectsQtlsCommon
                                        ),
                                        effectsQtlsAncestry,
                                        effectsQtlsGroup,
                                        effectsQtlsGbs,
                                        list(Polygenes = effectsQtlsPg)),
                                        True = c(list(
                                          Common = effectsQtlsCommonReal
                                        ),
                                        effectsQtlsAncestryReal,
                                        effectsQtlsGroupReal,
                                        effectsQtlsGbsReal,
                                        list(Polygenes = effectsQtlsPgReal))
                                      )
                                    ),
                                    Polygenes = c(list(
                                      Common = effectsPolyGenesCommon * correctionForGenValsPolyGenes
                                    ),
                                    lapply(X = effectsPolyGenesAncestry,
                                           FUN = function(x) {
                                             if (!is.null(x)) {
                                               x <- x * correctionForGenValsPolyGenes
                                             }

                                             return(x)
                                           }),
                                    lapply(X = effectsPolyGenesGroup,
                                           FUN = function(x) {
                                             if (!is.null(x)) {
                                               x <- x * correctionForGenValsPolyGenes
                                             }

                                             return(x)
                                           }))
                                  )


                                  genValsList <- list(
                                    True = genVals,
                                    Total = cbind(HaploQtls = genValsHaploQtls,
                                                  Qtls = genValsQtls,
                                                  QtlsAll = genValsQtlsAll,
                                                  Polygenes = genValsPolyGenes),
                                    Each = list(
                                      HaploQtls = genValsHaploQtlsMat,
                                      Qtls = genValsQtlsMat,
                                      QtlsAll = genValsQtlsAllMat,
                                      Polygenes = genValsPolyGenesMat
                                    )
                                  )



                                  #### 3.2.10. Generate residuals from heritability ####
                                  varPropResid <- (1 - herit) / herit
                                  varResid <- var(genVals) * varPropResid

                                  set.seed(seed = seeds["Residuals_Effects"])
                                  resids <- rnorm(n = nLines,
                                                  mean = 0, sd = sqrt(varResid))
                                  names(resids) <- names(genVals)



                                  #### 3.2.11. Compute phenotypic values ####
                                  phenVals <- genVals + resids



                                  #### 3.2.12. Compute weights for variance components and summary ####
                                  weightsParam <- c(heritQtls, heritPolyGenes, 1 - herit)
                                  weightsParam <- weightsParam / sum(weightsParam)
                                  names(weightsParam) <- c("QtlsAll", "Polygenes", "Residuals")

                                  varHaploQtls <- var(genValsHaploQtls)
                                  varQtls <- var(genValsQtls)
                                  varQtlsAll <- var(genValsQtlsAll)
                                  varPolyGenes <- var(genValsPolyGenes)
                                  varGenotypicValues <- var(genVals)
                                  varResids <- var(resids)
                                  varPhens <- var(phenVals)
                                  weightsTrue <- c(varQtlsAll, varPolyGenes, varResids)
                                  weightsTrue <- weightsTrue / sum(weightsTrue)
                                  names(weightsTrue) <- c("QtlsAll", "Polygenes", "Residuals")


                                  varInfo <- list(
                                    Param = list(
                                      heritability = c(
                                        Qtls = heritQtls,
                                        PolyGenes = heritPolyGenes,
                                        herit = herit
                                      ),
                                      weights = weightsParam
                                    ),
                                    True = list(
                                      variance = c(
                                        HaploQtls = varHaploQtls,
                                        Qtls = varQtls,
                                        QtlsAll = varQtlsAll,
                                        Polygenes = varPolyGenes,
                                        GenotypicValues = varGenotypicValues,
                                        Residuals = varResids,
                                        PhenotypeValues = varPhens
                                      ),
                                      heritability = c(
                                        HaploQtls = varHaploQtls / varPhens,
                                        Qtls = varQtls / varPhens,
                                        QtlsAll = varQtlsAll / varPhens,
                                        PolyGenes = varPolyGenes / varPhens,
                                        herit = varGenotypicValues / varPhens
                                      ),
                                      weights = weightsTrue
                                    )
                                  )



                                  #### 3.2.13. Summary inputted parameters and save/return results ####
                                  inputParams <- list(
                                    gbInfo = gbInfo,
                                    gbsMat = gbsMat,
                                    nPolyGenesCommon = nPolyGenesCommon,
                                    nPolyGenesAncestry = nPolyGenesAncestry,
                                    nPolyGenesGroup = nPolyGenesGroup,
                                    nQtlsCommon = nQtlsCommon,
                                    nQtlsAncestry = nQtlsAncestry,
                                    nQtlsGroup = nQtlsGroup,
                                    nQtlsGbs = nQtlsGbs,
                                    nQtlsPg = nQtlsPg,
                                    nHaploQtlsCommon = nHaploQtlsCommon,
                                    nHaploQtlsAncestry = nHaploQtlsAncestry,
                                    nHaploQtlsGroup = nHaploQtlsGroup,
                                    nHaploQtlsGbs = nHaploQtlsGbs,
                                    nHaploQtlsPg = nHaploQtlsPg,
                                    qtlCandNamesGiven = qtlCandNamesGiven,
                                    nullChrNo = nullChrNo,
                                    effectsAncestryType = effectsAncestryType,
                                    effectsGroupType = effectsGroupType,
                                    effectsQtlsType = effectsQtlsType,
                                    heritQtls = heritQtls,
                                    herit = herit,
                                    seeds = seeds
                                  )

                                  results <- list(
                                    phenVals = phenVals,
                                    resids = resids,
                                    genVals = genVals,
                                    genValsList = genValsList,
                                    mrkNamesList = mrkNamesList,
                                    effectsList = effectsList,
                                    varInfo = varInfo,
                                    inputParams = inputParams
                                  )


                                  if (!is.null(saveAt)) {
                                    saveRDS(object = results,
                                            file = fileNameIterAllRes)
                                  }
                                } else {
                                  results <- readRDS(file = fileNameIterAllRes)
                                }

                                # if (returnAll) {
                                #   returnRes <- results
                                # } else {
                                #   returnRes <- list(phenVals = results$phenVals,
                                #                     seeds = results$inputParams$seeds)
                                #
                                # }
                                # iterNo <- as.numeric(stringr::str_remove(string = iterName,
                                #                               pattern = iterNameHead))
                                # if (iterNo %% 100 == 0) {
                                #   gc(reset = TRUE); gc(reset = TRUE)
                                # }


                                return(results)
                              },
                              mc.cores = nCores)
      )
      cat("\n")
      gc(reset = TRUE); gc(reset = TRUE)
    }




    ##### 3.3. Summary results to be returned as the output of function #####
    if (returnAll) {
      results <- resultsAllList
      names(results) <- iterNames
    } else {
      results <- do.call(what = cbind,
                         args = lapply(X = resultsAllList,
                                       FUN = function(results) {
                                         return(results$phenVals)
                                       }))
      dimnames(results) <- list(lineNames,
                                iterNames)
    }

    if (returnSeeds) {
      if (!(is.null(saveAt) | overWrite)) {
        seedsMat <- do.call(what = rbind,
                            args = lapply(X = resultsAllList,
                                          FUN = function(results) {
                                            return(results$inputParams$seeds)
                                          }))
        dimnames(seedsMat) <- list(iterNames, seedTypeNames)
      }
      resultsAll <- list(results = results,
                         seeds = seedsMat)
    } else {
      resultsAll <- results
    }




    return(resultsAll)
  }












  ###### 4. Simulation of phenotypic data ######
  simulatedRes <- simulatePhenVals(genoScore = genoScore,
                                   genoMap = genoMap,
                                   ancestralHaploDf = ancestralHaploDf,
                                   haploBlockDf = haploBlockDataFrame,
                                   gbInfo = genBackVec,
                                   gbsMat = pcs,
                                   nIter = inputParams$nIter,
                                   nPolyGenesCommon = inputParams$nPolyGenesCommon,
                                   nPolyGenesAncestry = inputParams$nPolyGenesAncestry,
                                   nPolyGenesGroup = inputParams$nPolyGenesGroup,
                                   nQtlsCommon = inputParams$nQtlsCommon,
                                   nQtlsAncestry = inputParams$nQtlsAncestry,
                                   nQtlsGroup = inputParams$nQtlsGroup,
                                   nQtlsGbs = inputParams$nQtlsGbs,
                                   nQtlsPg = inputParams$nQtlsPg,
                                   nHaploQtlsCommon = inputParams$nHaploQtlsCommon,
                                   nHaploQtlsAncestry = inputParams$nHaploQtlsAncestry,
                                   nHaploQtlsGroup = inputParams$nHaploQtlsGroup,
                                   nHaploQtlsGbs = inputParams$nHaploQtlsGbs,
                                   nHaploQtlsPg = inputParams$nHaploQtlsPg,
                                   nMarkersInBlocksMin = inputParams$nMarkersInBlocksMin,
                                   qtlCandNamesGiven = qtlCandNamesGiven,
                                   nullChrNos = inputParams$nullChrNos,
                                   effectsAncestryType = inputParams$effectsAncestryType,
                                   effectsGroupType = inputParams$effectsGroupType,
                                   effectsQtlsType = inputParams$effectsQtlsType,
                                   heritQtls = inputParams$heritQtls,
                                   herit = inputParams$herit,
                                   saveAt = dirTrial,
                                   overWrite = overWrite,
                                   iterNameHead = "Iteration_",
                                   returnAll = FALSE,
                                   nCores = nCores,
                                   seedsMat = NULL,
                                   returnSeeds = TRUE,
                                   nIterOneBlock = 100,
                                   scoring012 = TRUE)
  phenoMat <- simulatedRes$results
  seedsMat <- simulatedRes$seeds


  fileNameSimPhenoMatData <- paste0("data/phenotype/",
                                    scriptID, "_simulated_phenotype_matrix_Scenario_",
                                    scenarioNo, "_Trial_", trialNo, ".csv")
  fileNameSimPhenoMatMid <- paste0(dirTrial,
                                   scriptID, "_simulated_phenotype_matrix_Scenario_",
                                   scenarioNo, "_Trial_", trialNo, ".csv")
  fileNameSimSeedMatMid <- paste0(dirTrial,
                                  scriptID, "_seeds_for_simulation_Scenario_",
                                  scenarioNo, "_Trial_", trialNo, ".csv")


  fwrite(x = data.frame(phenoMat),
         file = fileNameSimPhenoMatData,
         row.names = TRUE)

  fwrite(x = data.frame(phenoMat),
         file = fileNameSimPhenoMatMid,
         row.names = TRUE)

  fwrite(x = data.frame(seedsMat),
         file = fileNameSimSeedMatMid,
         row.names = TRUE)
}
