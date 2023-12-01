######################################################################################################################
######  Title: 2.2_Soybean_DLIPG_Perform_HB_based_GWAS_for_simulated_phenotypes_including_haplotype_based_QTLs  ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                                                        ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo                               ######
######  Date: 2022/06/28 (Created), 2022/07/13 (Last Updated)                                                   ######
######################################################################################################################



###### 1. Settings ######
##### 1.0. Reset workspace ######
# rm(list=ls())



##### 1.1. Setting working directory to the "DLIPG" directory #####
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
scriptIDSimData <- "0.6"

scriptID <- "2.2"







##### 1.2. Setting some parameters #####
dirMidDLIPGBase <- "midstream/"
nCores <- 55
# nCores <- 10
scenarioNo <- 4
trialNo <- 1

nPcsDapc <- 10
mafThres <- 0.05
ldThres <- 0.95
useImputedBlock <- TRUE
includeEpistasisInHBxGB <- TRUE


modelMethodNow <- c("HB", "HBxGB")
methodParamList <- list(
  list(method = "HB", param = NULL),
  list(method = "HBxGB", param = NULL)
  # ,
  # list(method = "HBxGB", param = nPcsDapc),
  # list(method = "HBxGB", param = 970)
)


colQtls <- c("black", "grey",
             "red", "blue", "green",
             "purple", "lightgreen", "orange",
             "grey", "brown")
names(colQtls) <- c("Common", "Ancestry",
                    "Group-China", "Group-Japan", "Group-Korea",
                    "Group-C-J", "Group-J-K", "Group-K-C",
                    "PC1", "Polygenes")
pchQtls <- c(16, 17)
cexQtls <- c(2.5, 2.5)
names(pchQtls) <- names(cexQtls) <- c("Qtls", "HaploQtls")


dirSimData <- paste0(dirMidDLIPGBase,
                     scriptIDSimData, "_simulated_data/")
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



# fileParamsSimData <- paste0(dirTrial, scriptID,
#                             "_", project, "_simulation_parameters.RData")
# saveRDS(object = inputParams,
#         file = fileParamsSimData)






##### 1.3. Import packages #####
require(data.table)
require(MASS)
require(BGLR)
require(RAINBOWR)
require(ggplot2)




dirScript <- "scripts/"
# dirScript <- paste0(dirScriptBase, cropName, "/Project/", project, "/")
# source(paste0(dirScript, "1.0_Soybean_DLIPG_function_to_draw_original_pair_plots_for_GWAS_results.R"))
# source(paste0(dirScript, "1.1_Soybean_DLIPG_function_to_perform_GWAS_by_Rio_model_with_MM4LMM.R"))
# source(paste0(dirScript, "1.2_Soybean_DLIPG_function_to_perform_GWAS_by_SNP_interaction_model.R"))
source(paste0(dirScript, "1.3_Soybean_DLIPG_function_to_perform_GWAS_by_HB_interaction_model.R"))
source(paste0(dirScript, "1.3.1_Soybean_DLIPG_function_to_define_haplotype_block_for_each_marker.R"))
# source(paste0(dirScript, "scriptName2"))








##### 1.4. Project options #####
options(stringAsFactors = FALSE)






###### 2. Read data into R and preparation for GWAS ######
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



##### 2.2. Read genotype map into R and compute cumulative position #####
fileNameGenoMap <- "data/extra/L900_Soybean_a1_chrnum_maf0.025_ms0.05_ChKoJa_300x3_maf0.025_imputed_maf0.025_SNP_geno_map.csv"
genoMap <- data.frame(fread(input = fileNameGenoMap))
colnames(genoMap)[1] <- "marker"

cumPos <- cumsumPos(map = genoMap)
names(cumPos) <- genoMap$marker


##### 2.3. Read marker genotype into R #####
fileNameGenoScore <- paste0("data/genotype/L900_Soybean_a1_chrnum_maf0.025_ms0.05_ChKoJa_300x3_maf0.025_imputed_maf0.025_SNP_geno.tsv")

genoScoreDf <- data.frame(fread(input = fileNameGenoScore,
                                check.names = FALSE),
                          check.names = FALSE)
genoScore <- t(genoScoreDf[, -c(1:2)])
colnames(genoScore) <- genoMap$marker

genBackVec <- genBackVecRaw[rownames(genoScore)]



##### 2.4. Read additive genomic relationship matrix into R #####
fileNameAddGrm <- paste0("data/extra/", scriptIDPca, "_",
                         "additive_genomic_relationship_",
                         paste(genBackCombAll[[length(genBackCombAll)]],
                               collapse = "_"), ".rds")
fileNameAddGrmAdj <- paste0("data/extra/", scriptIDPca, "_",
                            "adjusted_additive_genomic_relationship_",
                            paste(genBackCombAll[[length(genBackCombAll)]],
                                  collapse = "_"), ".rds")
addGrm <- readRDS(file = fileNameAddGrm)
addGrmAdj <- readRDS(file = fileNameAddGrmAdj)



##### 2.5. Read genotype matrix including ancestry information into R #####
# if ("Rio" %in% modelMethodNow) {
#   ancestralHaploMatHead <- apply(X = ancestralHaploMat,
#                                  MARGIN = 2,
#                                  FUN = stringr::str_sub,
#                                  start = 1,
#                                  end = 1)
#   rownames(ancestralHaploMatHead) <- rownames(genoScore)
#   genoScoreWithGenBack <- matrix(data = paste0(genoScore + 1,
#                                                ancestralHaploMatHead),
#                                  nrow = nrow(genoScore),
#                                  ncol = ncol(genoScore),
#                                  byrow = FALSE,
#                                  dimnames = dimnames(genoScore))
# }



##### 2.6. Read simulated phenotypic data into R #####
fileNameSimPhenoMatData <- paste0("data/phenotype/",
                                  scriptIDSimData, "_simulated_phenotype_matrix_Scenario_",
                                  scenarioNo, "_Trial_", trialNo, ".csv")
phenoMat <- data.frame(fread(input = fileNameSimPhenoMatData), row.names = 1)
nTraits <- ncol(phenoMat)
traitNames <- colnames(phenoMat)



##### 2.7. Modify marker genotype and phenotype data for RAINBOWR #####
modifyDataRes <- RAINBOWR::modify.data(pheno.mat = phenoMat,
                                       geno.mat = genoScore,
                                       map = genoMap,
                                       return.ZETA = FALSE,
                                       return.GWAS.format = TRUE)
phenoGWAS <- modifyDataRes$pheno.GWAS
genoGWAS <- modifyDataRes$geno.GWAS



##### 2.8. Read haplotype block list for RAINBOWR into R #####
if (!useImputedBlock) {
  haploBlockListDf <- data.frame(data.table::fread(input = paste0(
    "data/extra/", scriptIDHaploBlock,
    "_haplotype_block_data_for_RAINBOWR.csv"
  )))
  rownames(haploBlockListDf) <- haploBlockListDf$marker

  fileNameHaploBlockMapWithCumPos <- paste0("data/extra/", scriptID,
                                            "_genotype_map_with_cumulative_position_",
                                            "for_haplotype_block_data.csv")
} else {
  haploBlockListDf <- data.frame(data.table::fread(input = paste0(
    "data/extra/", scriptIDHaploBlock,
    "_haplotype_block_data_imputed_for_RAINBOWR.csv"
  )))
  rownames(haploBlockListDf) <- haploBlockListDf$marker

  fileNameHaploBlockMapWithCumPos <- paste0("data/extra/", scriptID,
                                            "_genotype_map_with_cumulative_position_",
                                            "for_haplotype_block_data_imputed.csv")
}

if (!file.exists(fileNameHaploBlockMapWithCumPos)) {
  haploBlockMapWithCumPos <- genesetmap(map = genoMap,
                                        gene.set = haploBlockListDf,
                                        cumulative = TRUE)
  data.table::fwrite(x = haploBlockMapWithCumPos, file = fileNameHaploBlockMapWithCumPos)
} else {
  haploBlockMapWithCumPos <- data.frame(data.table::fread(
    input = fileNameHaploBlockMapWithCumPos
  ))
}




###### 3. Perform GWAS for simulated dataset ######
##### 3.1. Perform GWAS for each trait and prepare directory #####
gbInfoUsedList <- rep(list(NULL), length(genBackCombAll))
names(gbInfoUsedList) <- names(genBackCombAll)
for (traitNo in 1:nTraits) {
  # traitNo <- 1
  startTrait <- Sys.time()
  traitName <- traitNames[traitNo]

  dirIterSave <- paste0(dirTrial, traitName, "/")
  if (!dir.exists(dirIterSave)) {
    dir.create(dirIterSave)
  }



  ##### 3.2. Read simulation information to generate phenotypes into R #####
  fileNameIterAllRes <- paste0(dirIterSave,
                               "simulated_phenotype_all_results_",
                               traitName,".rds")

  phenoSimInfo <- readRDS(file = fileNameIterAllRes)


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


  ##### 3.3. Perform GWAS with each genetic background and prepape for some parameters  #####
  for (genBackEach in genBackCombAll) {
    startGb <- Sys.time()
    # genBackEach <- genBackCombAll[[4]]
    genBackNameEach <- paste(genBackEach, collapse = "_")

    addGrmSubGr <- readRDS(file = paste0("data/extra/", scriptIDPca,
                                         "_additive_genomic_relationship_",
                                         genBackNameEach, ".rds"))
    addGrmAdjSubGr <- readRDS(file = paste0("data/extra/", scriptIDPca,
                                            "_adjusted_additive_genomic_relationship_",
                                            genBackNameEach, ".rds"))


    if ((length(genBackEach) >= 2) & is.null(gbInfoUsedList[[genBackNameEach]])) {
      fcRes <- adegenet::find.clusters(x = addGrmSubGr,
                                       n.pca = nPcsDapc,
                                       n.clust = length(genBackEach))
      gbInfoUsedList[[genBackNameEach]] <- gbInfoUsed <- fcRes$grp
    } else {
      gbInfoUsed <- gbInfoUsedList[[genBackNameEach]]
    }

    # if (length(genBackEach) == 1) {
    #   methodParamListNow <- methodParamList[sapply(methodParamList,
    #                                                function(x) x$method) %in% c("Rio", "PCA")]
    # } else {
    #   methodParamListNow <- lapply(X = methodParamList,
    #                                FUN = function(x) {
    #                                  if (is.null(x$param)) {
    #                                    x$param <- length(genBackEach)
    #                                  }
    #
    #                                  return(x)
    #                                })
    # }
    methodParamListNow <- lapply(X = methodParamList,
                                 FUN = function(x) {
                                   if (!is.null(x$param)) {
                                     x$param <- min(nrow(addGrmSubGr),
                                                    x$param)
                                   }
                                   return(x)
                                 })
    methodParamListPerform <- unlist(lapply(X = methodParamListNow,
                                            FUN = function(x) {
                                              x$method %in% modelMethodNow
                                            }))
    methodParamListNow <- methodParamListNow[methodParamListPerform]

    if (includeEpistasisInHBxGB) {
      includeEpiMethods <- sapply(X = methodParamListNow,
                                  FUN = function(x) (x$method == "HBxGB") & is.null(x$param))
    } else {
      includeEpiMethods <- rep(FALSE, length(methodParamListNow))
    }

    if ((length(genBackEach) >= 2) | any(includeEpiMethods)) {
      if (any(includeEpiMethods)) {
        fileNameAddGrmSubGrEpi <- paste0(dirIterSave, scriptID, "_",
                                         "additive_genomic_relationship_adjusted_by_estimated_subgroups_",
                                         genBackNameEach, "_include_epistasis.rds")

        if (!file.exists(fileNameAddGrmSubGrEpi)) {
          adjustGrmResEpi <- RAINBOWR::adjustGRM(y = phenoMat[genBackVec %in% genBackEach, traitNo, drop = FALSE],
                                                 X = NULL,
                                                 ZETA = list(A = list(Z = design.Z(pheno.labels = rownames(phenoMat)[genBackVec %in% genBackEach],
                                                                                   geno.names = rownames(addGrmSubGr)),
                                                                      K = addGrmSubGr)),
                                                 subpopInfo = gbInfoUsed,
                                                 nSubpop = length(genBackEach),
                                                 nPcsFindCluster = nPcsDapc,
                                                 include.epistasis = TRUE,
                                                 package.MM = "gaston")

          saveRDS(object = adjustGrmResEpi, file = fileNameAddGrmSubGrEpi)
        } else {
          adjustGrmResEpi <- readRDS(file = fileNameAddGrmSubGrEpi)
        }
        addGrmSubGrEpi <- adjustGrmResEpi$ZETAAdjust$Adjust$K
      }
    }


    if ((length(genBackEach) >= 2) & any(!includeEpiMethods)) {
      fileNameAddGrmSubGrNonEpi <- paste0(dirIterSave, scriptID, "_",
                                          "additive_genomic_relationship_adjusted_by_estimated_subgroups_",
                                          genBackNameEach, "_without_epistasis.rds")

      if (!file.exists(fileNameAddGrmSubGrNonEpi)) {
        adjustGrmResNonEpi <- RAINBOWR::adjustGRM(y = phenoMat[genBackVec %in% genBackEach, traitNo, drop = FALSE],
                                                  X = NULL,
                                                  ZETA = list(A = list(Z = design.Z(pheno.labels = rownames(phenoMat)[genBackVec %in% genBackEach],
                                                                                    geno.names = rownames(addGrmSubGr)),
                                                                       K = addGrmSubGr)),
                                                  subpopInfo = gbInfoUsed,
                                                  nSubpop = length(genBackEach),
                                                  nPcsFindCluster = nPcsDapc,
                                                  include.epistasis = FALSE,
                                                  package.MM = "gaston")
        saveRDS(object = adjustGrmResNonEpi, file = fileNameAddGrmSubGrNonEpi)
      } else {
        adjustGrmResNonEpi <- readRDS(file = fileNameAddGrmSubGrNonEpi)
      }
      addGrmSubGrNonEpi <- adjustGrmResNonEpi$ZETAAdjust$Adjust$K
    } else {
      addGrmSubGrNonEpi <- addGrmSubGr
    }





    ##### 3.4. Perform HB GWAS using RAINBOW model with different parameters  #####
    for (methodParams in methodParamListNow) {
      # methodParams <- methodParamListNow[[2]]
      if (methodParams$method == "Rio") {
        modelName <- paste0(methodParams$method,
                            "-M", length(genBackEach))
      } else {
        if (!is.null(methodParams$param)) {
          modelName <- paste0(methodParams$method,
                              "-", methodParams$param)
        } else {
          modelName <- methodParams$method
        }
      }

      print(paste0("Trait: ", traitName,
                   "   Genetic background: ", genBackNameEach,
                   "   Model: ", modelName))

      if (!useImputedBlock) {
        fileNameResGWAS <- paste0(dirIterSave, scriptID, "_",
                                  "GWAS_results_for_", traitName,
                                  "_with_", genBackNameEach,
                                  "_by_", modelName, "_model.rds")
      } else {
        fileNameResGWAS <- paste0(dirIterSave, scriptID, "_",
                                  "GWAS_results_for_", traitName,
                                  "_with_", genBackNameEach,
                                  "_by_", modelName, "_model_imputed.rds")
      }


      #### 3.4.1. Perform HB GWAS ####
      if (includeEpistasisInHBxGB) {
        includeEpistasis <- (methodParams$method == "HBxGB") & is.null(methodParams$param)
      } else {
        includeEpistasis <- FALSE
      }

      if (!file.exists(fileNameResGWAS)) {
        start <- Sys.time()
        if (includeEpistasis) {
          addGrmSubGrNow <- addGrmSubGrEpi
        } else {
          addGrmSubGrNow <- addGrmSubGrNonEpi
        }

        resGWAS <- performGWASHBInt(genoGWAS = genoGWAS,
                                    phenoGWASNow = phenoGWAS[, c(1, traitNo + 1), drop = FALSE],
                                    haploBlockListDf = haploBlockListDf,
                                    mapHaploBlock = haploBlockMapWithCumPos,
                                    interactionKernelOrigin = addGrmSubGr,
                                    addGRM = addGrm,
                                    addGRMSubGr = addGrmSubGrNow,
                                    gbInfo = genBackVec,
                                    gbInfoUsed = gbInfoUsed,
                                    gbNow = genBackEach,
                                    interactionMatMethod = methodParams$method,
                                    paramForIntGWAS = methodParams$param,
                                    nPcsDapc = nPcsDapc,
                                    nCores = nCores,
                                    mafThres = mafThres,
                                    ldThres = ldThres,
                                    adjustAddGRM = FALSE)


        #### 3.4.2. Save the GWAS results ####
        saveRDS(object = resGWAS,
                file = fileNameResGWAS)




        ##### 3.5. Preparation fordrawing manhattan plots and compute thresholds #####
        nTargets <- ncol(resGWAS) - 3
        if (nTargets >= 2) {
          targetNames <- colnames(resGWAS)[-c(1:3)]
        } else {
          targetNames <- genBackNameEach
        }

        thresholds <- sapply(
          X = resGWAS[, -c(1:3), drop = FALSE],
          FUN = function(x) {
            thres <- try(RAINBOWR::CalcThreshold(input = cbind(haploBlockMapWithCumPos[, 1:3], x),
                                                 sig.level = 0.05,
                                                 method = "BH"), silent = TRUE)
            if ("try-error" %in% class(thres)) {
              thres <- NA
            }
            return(thres)
          }
        )
        names(thresholds) <- targetNames



        ##### 3.6. Save the results for real qtls in simulation #####
        if (!useImputedBlock) {
          fileNameResQtl <- paste0(dirIterSave, scriptID, "_",
                                   "GWAS_results_of_QTLs_for_", traitName,
                                   "_with_", genBackNameEach,
                                   "_by_", modelName, "_model.rds")
        } else {
          fileNameResQtl <- paste0(dirIterSave, scriptID, "_",
                                   "GWAS_results_of_QTLs_for_", traitName,
                                   "_with_", genBackNameEach,
                                   "_by_", modelName, "_model_imputed.rds")
        }

        haploBlocksIncludeQtlsListQtls <- lapply(X = qtlNamesListQtls,
                                                 FUN = function(qtlNames) {
                                                   if (!is.null(qtlNames)) {
                                                     haploBlocksIncludeQtl <- sapply(
                                                       X = qtlNames,
                                                       FUN = function(mrkName) {
                                                         defineHaploBlock(mrkName = mrkName,
                                                                          haploBlockListDf = haploBlockListDf,
                                                                          genoMap = genoMapCommon,
                                                                          resGWAS = resGWAS)
                                                       },
                                                       simplify = TRUE
                                                     )
                                                   } else {
                                                     haploBlocksIncludeQtl <- NULL
                                                   }

                                                   return(haploBlocksIncludeQtl)
                                                 })
        qtlResListQtls <- lapply(X = haploBlocksIncludeQtlsListQtls,
                                 FUN = function(haploBlocks) {
                                   if (!is.null(haploBlocks)) {
                                     rownames(resGWAS) <- resGWAS$marker
                                     qtlRes <- resGWAS[haploBlocks, , drop = FALSE]
                                   } else {
                                     qtlRes <- NULL
                                   }
                                   return(qtlRes)
                                 })


        haploBlocksIncludeQtlsListHaploQtls <- lapply(X = qtlNamesListHaploQtls,
                                                      FUN = function(qtlNames) {
                                                        if (!is.null(qtlNames)) {
                                                          haploBlocksIncludeQtl <- sapply(
                                                            X = qtlNames,
                                                            FUN = function(mrkName) {
                                                              defineHaploBlock(mrkName = mrkName,
                                                                               haploBlockListDf = haploBlockListDf,
                                                                               genoMap = genoMapCommon,
                                                                               resGWAS = resGWAS)
                                                            },
                                                            simplify = TRUE
                                                          )
                                                        } else {
                                                          haploBlocksIncludeQtl <- NULL
                                                        }

                                                        return(haploBlocksIncludeQtl)
                                                      })
        qtlResListHaploQtls <- lapply(X = haploBlocksIncludeQtlsListHaploQtls,
                                      FUN = function(haploBlocks) {
                                        if (!is.null(haploBlocks)) {
                                          rownames(resGWAS) <- resGWAS$marker
                                          qtlRes <- resGWAS[haploBlocks, , drop = FALSE]
                                        } else {
                                          qtlRes <- NULL
                                        }
                                        return(qtlRes)
                                      })




        qtlResList <- c(list(HaploQtls = qtlResListHaploQtls),
                        list(Qtls = qtlResListQtls),
                        list(Thres = thresholds))
        saveRDS(object = qtlResList, file = fileNameResQtl)


        resGWAS <- cbind(resGWAS,
                         cumpos = haploBlockMapWithCumPos$cum.pos[haploBlockMapWithCumPos$marker %in% resGWAS$marker])


        ##### 3.7. Draw qq and manhattan plots with true QTL positions for each target #####
        for (targetNo in 1:nTargets) {
          targetName <- targetNames[targetNo]
          resGWASTarget <- resGWAS[, c(1:3, targetNo + 3, ncol(resGWAS))]

          if (!useImputedBlock) {
            # fileNameResGWASOrd <- paste0(dirIterSave, scriptID, "_",
            #                              "GWAS_results_for_", traitName,
            #                              "_with_", genBackNameEach,
            #                              "_by_", modelName, "_model_",
            #                              targetName, "-targeted_ordered.csv")
            fileNameResGWASManhattan <- paste0(dirIterSave, scriptID, "_",
                                               "GWAS_results_for_", traitName,
                                               "_with_", genBackNameEach,
                                               "_by_", modelName, "_model_",
                                               targetName, "-targeted_manhattan.png")
            fileNameResGWASQq <- paste0(dirIterSave, scriptID, "_",
                                        "GWAS_results_for_", traitName,
                                        "_with_", genBackNameEach,
                                        "_by_", modelName, "_model_",
                                        targetName, "-targeted_qq.png")
          } else {
            # fileNameResGWASOrd <- paste0(dirIterSave, scriptID, "_",
            #                              "GWAS_results_for_", traitName,
            #                              "_with_", genBackNameEach,
            #                              "_by_", modelName, "_model_",
            #                              targetName, "-targeted_ordered_imputed.csv")
            fileNameResGWASManhattan <- paste0(dirIterSave, scriptID, "_",
                                               "GWAS_results_for_", traitName,
                                               "_with_", genBackNameEach,
                                               "_by_", modelName, "_model_",
                                               targetName, "-targeted_manhattan_imputed.png")
            fileNameResGWASQq <- paste0(dirIterSave, scriptID, "_",
                                        "GWAS_results_for_", traitName,
                                        "_with_", genBackNameEach,
                                        "_by_", modelName, "_model_",
                                        targetName, "-targeted_qq_imputed.png")
          }

          # targetNo <- 1
          # resGWASTargetOrd <- resGWASTarget[order(resGWASTarget[, 4], decreasing = TRUE), ]
          # fwrite(x = resGWASTargetOrd, file = fileNameResGWASOrd, row.names = TRUE)

          png(filename = fileNameResGWASQq, width = 800, height = 800)
          RAINBOWR::qq(scores = resGWASTarget[, 4])
          title(main = paste0(traitName, " - ", modelName,
                              " - ", targetName))
          dev.off()


          # png(filename = fileNameResGWASManhattan, width = 1450, height = 900)
          # par(mar = c(5.1, 4.1, 4.1, 7.1))
          # RAINBOWR::manhattan(input = resGWASTarget)
          # title(main = paste0(traitName, " - ", modelName,
          #                     " - ", targetName))
          # for (qtlNamesNo in 1:length(qtlNamesList)) {
          #   # qtlNamesNo <- 1
          #   qtlNames <- qtlNamesList[[qtlNamesNo]]
          #   haploBlockQtls <- haploBlocksIncludeQtlsList[[qtlNamesNo]]
          #   colNow <- colQtls[names(qtlNamesList)[qtlNamesNo]]
          #
          #   # qtlResNow <- cbind(resGWASTarget[resGWAS$marker %in% haploBlockQtls, 4],
          #   #                    cumPosCommon[qtlNames])
          #   qtlResNow <- cbind(resGWASTarget[resGWAS$marker %in% haploBlockQtls, 4:5])
          #   for (qtlNo in 1:nrow(qtlResNow)) {
          #     abline(v = qtlResNow[qtlNo, 2],
          #            col = colNow,
          #            lty = 2,
          #            lwd = 1.8)
          #     points(x = qtlResNow[qtlNo, 2],
          #            y = qtlResNow[qtlNo, 1],
          #            col = colNow,
          #            pch = 16,
          #            cex = 2.5)
          #   }
          # }
          # par(xpd = TRUE)
          # legend(par()$usr[2], par()$usr[4],
          #        legend = names(qtlNamesList),
          #        col = colQtls[names(qtlNamesList)],
          #        pch = 16)
          # # legend("topleft",
          # #        legend = names(qtlNamesList),
          # #        col = colQtls,
          # #        pch = 16)
          # dev.off()
          # par(xpd = FALSE, mar = c(5.1, 4.1, 4.1, 2.1))



          png(filename = fileNameResGWASManhattan, width = 1500, height = 950)
          par(mar = c(7.1, 6.1, 4.1, 13.1))
          RAINBOWR::manhattan(input = resGWASTarget,
                              cex.lab = 2,
                              cex.axis.x = 1.8,
                              cex.axis.y = 1.8)
          title(main = paste0(traitName, " - ", modelName,
                              " - ", targetName),
                cex.main = 3)


          for (qtlTypeNo in 1:length(qtlNamesListAll)) {
            # qtlTypeNo <- 1
            qtlNamesListNow <- qtlNamesListAll[[qtlTypeNo]]
            qtlTypeName <- names(qtlNamesListAll)[qtlTypeNo]
            pchNow <- pchQtls[qtlTypeName]
            if (!is.na(as.numeric(pchNow))) {
              pchNow <- as.numeric(pchNow)
            }
            cexNow <- cexQtls[qtlTypeName]

            if (qtlTypeNo == 1) {
              haploBlocksIncludeQtlsList <- haploBlocksIncludeQtlsListQtls
            } else {
              haploBlocksIncludeQtlsList <- haploBlocksIncludeQtlsListHaploQtls
            }

            for (qtlNamesNo in 1:length(qtlNamesListNow)) {
              # qtlNamesNo <- 1
              qtlNames <- qtlNamesListNow[[qtlNamesNo]]
              haploBlockQtls <- haploBlocksIncludeQtlsList[[qtlNamesNo]]
              colNow <- colQtls[names(qtlNamesListNow)[qtlNamesNo]]

              # qtlResNow <- resGWAS[resGWAS$marker %in% qtlNames,
              #                      c(targetNo + 3, ncol(resGWAS))]
              qtlResNow <- cbind(resGWASTarget[resGWAS$marker %in% haploBlockQtls, 4:5])
              for (qtlNo in 1:nrow(qtlResNow)) {
                abline(v = qtlResNow[qtlNo, 2],
                       col = colNow,
                       lty = 2,
                       lwd = 1.8)
                points(x = qtlResNow[qtlNo, 2],
                       y = qtlResNow[qtlNo, 1],
                       col = colNow,
                       pch = pchNow,
                       cex = cexNow)
              }
            }
          }
          par(xpd = TRUE)

          whichNull <- unlist(lapply(X = qtlNamesListAll,
                                     FUN = function(qtlNamesListNow) {
                                       unlist(lapply(X = qtlNamesListNow, FUN = is.null))
                                     }))
          qtlTypeNamesAll <- paste0(rep(names(qtlNamesListAll), each = length(qtlNamesListNow)),
                                    "-", rep(names(qtlNamesListNow), length(qtlNamesListAll)))
          legend(par()$usr[2], par()$usr[4],
                 legend = qtlTypeNamesAll[!whichNull],
                 col = rep(colQtls[names(qtlNamesListNow)], length(qtlNamesListAll))[!whichNull],
                 pch = rep(pchQtls, each = length(qtlNamesListNow))[!whichNull])

          dev.off()
          par(xpd = FALSE, mar = c(5.1, 4.1, 4.1, 2.1))

        }

        end <- Sys.time()
        compTime <- end - start

        rm(resGWAS)
        gc(reset = TRUE); gc(reset = TRUE)

        print(paste0("Computation time: ", round(compTime, 3),
                     " ", attr(compTime, "units"), " for this GWAS model"))
      }
    }

    endGb <- Sys.time()
    compTimeGb <- endGb - startGb
    print(paste0("Computation time: ", round(compTimeGb, 3),
                 " ", attr(compTimeGb, "units"), " for this genetic background"))
    cat("---------------------------------------------------------\n\n")
  }

  endTrait <- Sys.time()
  compTimeTrait <- endTrait - startTrait
  print(paste0("Computation time: ", round(compTimeTrait, 3),
               " ", attr(compTimeTrait, "units"), " for this trait"))
  cat("=========================================================\n\n\n")
}
