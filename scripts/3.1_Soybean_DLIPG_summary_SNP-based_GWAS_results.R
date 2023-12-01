#############################################################################################
######  Title: 3.1_Soybean_DLIPG_summary_SNP-based_GWAS_results                  ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                               ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo      ######
######  Date: 2022/07/11 (Created), 2023/02/05 (Last Updated)                          ######
#############################################################################################



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

scriptIDGWAS <- "2.1"
scriptIDGWASHB <- "2.2"
scriptID <- "3.1"





##### 1.2. Setting some parameters #####
for (nTraitsHere in seq(from = 20, to = 100, by = 20)) {
  # nTraitsHere <- 5

  # nCores <- 30
  # nCores <- 55
  scenarioNo <- 4
  trialNo <- 1

  modelMethodNow <- c("Rio", "GROUP", "PCA", "DAPC", "NI", "HB", "HBxGB")
  # modelMethodNow <- c("Rio", "GROUP", "PCA", "DAPC", "NI", "HB")
  # modelMethodNow <- c("Rio")
  # modelMethodNow <- c("GROUP", "PCA", "DAPC", "NI", "HB", "HBxGB")
  methodParamList <- list(
    list(method = "Rio", param = "true"),
    list(method = "GROUP", param = "true"),
    list(method = "GROUP", param = NULL),
    list(method = "PCA", param = 1),
    list(method = "PCA", param = 2),
    list(method = "DAPC", param = "true"),
    list(method = "DAPC", param = NULL),
    # list(method = "NI", param = 0),
    list(method = "NI", param = "true"),
    list(method = "NI", param = NULL)
    ,
    list(method = "HB", param = NULL),
    list(method = "HBxGB", param = NULL)
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



                      dirMidDLIPGBase <- "midstream/"
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


                      dirSummary <- paste0(dirTrial, scriptID, "_",
                                           "summarized_results/")
                      if (!dir.exists(dirSummary)) {
                        dir.create(dirSummary)
                      }

                      summaryParams <- list(
                        useGivenBlock = TRUE,
                        useHaploBlock = FALSE,
                        sigLevel = c(0.05, 0.01),
                        methodThres = "BH",
                        inflatorPlus = 2,
                        lenLD = 200000,
                        corThres = 0.35,
                        windowSize = 0,
                        nCores = 5,
                        returnAucRes = TRUE,
                        returnBlockRes = TRUE,
                        useImputedBlock = TRUE,
                        dropQtl = FALSE
                      )


                      fileParamsSummary <- paste0(dirSummary, scriptID,
                                                  "_", project, "_summarization_parameters.RData")
                      saveRDS(object = summaryParams,
                              file = fileParamsSummary)






                      ##### 1.3. Import packages #####
                      require(data.table)
                      require(MASS)
                      require(BGLR)
                      require(RAINBOWR)
                      require(ggplot2)




                      dirScript <- "scripts/"
                      # dirScript <- paste0(dirScriptBase, cropName, "/Project/", project, "/")
                      # source(paste0(dirScript, "1.0_Soybean_DLIPG_function_to_draw_original_pair_plots_for_GWAS_results.R"))
                      source(paste0(dirScript, "1.1_Soybean_DLIPG_function_to_perform_GWAS_by_Rio_model_with_MM4LMM.R"))
                      source(paste0(dirScript, "1.2_Soybean_DLIPG_function_to_perform_GWAS_by_SNP_interaction_model.R"))
                      source(paste0(dirScript, "1.3_Soybean_DLIPG_function_to_perform_GWAS_by_HB_interaction_model.R"))
                      source(paste0(dirScript, "1.3.1_Soybean_DLIPG_function_to_define_haplotype_block_for_each_marker.R"))
                      source(paste0(dirScript, "1.4_Soybean_DLIPG_function_to_summary_GWAS_results_haplotype_based.R"))
                      source(paste0(dirScript, "1.4.1_Soybean_DLIPG_function_to_help_plot_GWAS_results_haplotype_based.R"))
                      source(paste0(dirScript, "1.4.2_Soybean_DLIPG_function_to_drop_QTLs_from_GWAS_results_haplotype_based.R"))
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
                      # fileNameAddGrmAdj <- paste0("data/extra/", scriptIDPca, "_",
                      #                             "adjusted_additive_genomic_relationship_",
                      #                             paste(genBackCombAll[[length(genBackCombAll)]],
                      #                                   collapse = "_"), ".rds")
                      # addGrmAdj <- readRDS(file = fileNameAddGrmAdj)

                      if (summaryParams$dropQtl & any(modelMethodNow %in% c("HB", "HBxGB"))) {
                        fileNameAddGrm <- paste0("data/extra/", scriptIDPca, "_",
                                                 "additive_genomic_relationship_",
                                                 paste(genBackCombAll[[length(genBackCombAll)]],
                                                       collapse = "_"), ".rds")
                        addGrm <- readRDS(file = fileNameAddGrm)
                      }


                      ##### 2.5. Read genotype matrix including ancestry information into R #####


                      ##### 2.6. Read simulated phenotypic data into R #####
                      fileNameSimPhenoMatData <- paste0("data/phenotype/",
                                                        scriptIDSimData, "_simulated_phenotype_matrix_Scenario_",
                                                        scenarioNo, "_Trial_", trialNo, ".csv")
                      phenoMat <- data.frame(fread(input = fileNameSimPhenoMatData), row.names = 1)
                      nTraits <- ncol(phenoMat)
                      traitNames <- colnames(phenoMat)



                      ##### 2.7. Modify marker genotype and phenotype data for RAINBOWR #####
                      if (summaryParams$dropQtl & any(modelMethodNow %in% c("HB", "HBxGB"))) {
                        modifyDataRes <- RAINBOWR::modify.data(pheno.mat = phenoMat,
                                                               geno.mat = genoScore,
                                                               map = genoMap,
                                                               return.ZETA = FALSE,
                                                               return.GWAS.format = TRUE)
                        phenoGWAS <- modifyDataRes$pheno.GWAS
                        genoGWAS <- modifyDataRes$geno.GWAS
                      }



                      ##### 2.8. Read `idsAroundMrcList` into R #####
                      if (!summaryParams$useHaploBlock) {
                        fileNameIdsAroundMrkListSnp <- paste0("data/extra/", scriptIDHaploBlock,
                                                              "_ids_around_markers_list_for_long_LD_block_imputed.rds")
                      } else {
                        fileNameIdsAroundMrkListSnp <- paste0("data/extra/", scriptIDHaploBlock,
                                                              "_ids_around_markers_list_for_haplotype_block_imputed.rds")
                      }
                      if (summaryParams$useGivenBlock & file.exists(fileNameIdsAroundMrkListSnp)) {
                        idsAroundMrcListSnp <- readRDS(file = fileNameIdsAroundMrkListSnp)
                      } else {
                        idsAroundMrcListSnp <- NULL
                      }


                      ##### 2.9. Read `idsAroundMrcList` into R #####



                      ##### 2.10. Read haplotype block list for RAINBOWR into R #####
                      if (!summaryParams$useImputedBlock) {
                        fileNameIdsAroundMrkListHaplo <- paste0("data/extra/", scriptIDHaploBlock,
                                                                "_haplotype_block_list_in_LD_block.rds")


                        haploBlockListDf <- data.frame(data.table::fread(input = paste0(
                          "data/extra/", scriptIDHaploBlock,
                          "_haplotype_block_data_for_RAINBOWR.csv"
                        )))
                        rownames(haploBlockListDf) <- haploBlockListDf$marker

                        fileNameHaploBlockMapWithCumPos <- paste0("data/extra/", scriptIDGWASHB,
                                                                  "_genotype_map_with_cumulative_position_",
                                                                  "for_haplotype_block_data.csv")

                      } else {
                        fileNameIdsAroundMrkListHaplo <- paste0("data/extra/", scriptIDHaploBlock,
                                                                "_imputed_haplotype_block_list_in_LD_block.rds")


                        haploBlockListDf <- data.frame(data.table::fread(input = paste0(
                          "data/extra/", scriptIDHaploBlock,
                          "_haplotype_block_data_imputed_for_RAINBOWR.csv"
                        )))
                        rownames(haploBlockListDf) <- haploBlockListDf$marker

                        fileNameHaploBlockMapWithCumPos <- paste0("data/extra/", scriptIDGWASHB,
                                                                  "_genotype_map_with_cumulative_position_",
                                                                  "for_haplotype_block_data_imputed.csv")
                      }

                      if (summaryParams$useGivenBlock & file.exists(fileNameIdsAroundMrkListHaplo)) {
                        if (!summaryParams$useHaploBlock) {
                          idsAroundMrcListHaplo <- readRDS(file = fileNameIdsAroundMrkListHaplo)
                        } else {
                          if (!summaryParams$useImputedBlock) {
                            fileNameHaploVsHaploImp <- paste0("data/extra/", scriptIDHaploBlock,
                                                              "_correspondence_between_haplotype_and_imputed_haplotype_blocks.csv")
                            haploVsHaploImp <- data.frame(data.table::fread(input = fileNameHaploVsHaploImp))

                            idsAroundMrcListHaplo <- split(x = 1:nrow(haploVsHaploImp),
                                                           f = haploVsHaploImp$HBI)[unique(haploVsHaploImp$HBI)]
                          } else {
                            idsAroundMrcListHaplo <- split(x = 1:length(unique(haploBlockListDf$block)),
                                                           f = unique(haploBlockListDf$block))[unique(haploVsHaploImp$HBI)]
                          }
                        }
                      } else {
                        idsAroundMrcListHaplo <- NULL
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



                      ###### 3. Read GWAS results into R and summary GWAS results  ######
                      ##### 3.1. Summary GWAS results for each trait and prepare directory #####
                      if (nTraitsHere >= 10) {
                        fileSaveNumbers <- seq(from = 10, to = (nTraitsHere %/% 10) * 10, by = 10)
                      } else {
                        fileSaveNumbers <- nTraitsHere
                      }
                      fileNamesSummaryStatResListAll <- paste0(dirSummary, scriptID,
                                                               "_all_summary_statistics_results_",
                                                               fileSaveNumbers)
                      fileNamesNullChrResListAll <- paste0(dirSummary, scriptID,
                                                           "_chromosome_under_null_hypothesis_results_",
                                                           fileSaveNumbers)
                      if (summaryParams$useGivenBlock) {
                        if (!summaryParams$useHaploBlock) {
                          fileNamesSummaryStatResListAll <- paste0(fileNamesSummaryStatResListAll,
                                                                   "_given_LD_block")
                          fileNamesNullChrResListAll <- paste0(fileNamesNullChrResListAll,
                                                               "_given_LD_block")
                        } else {
                          fileNamesSummaryStatResListAll <- paste0(fileNamesSummaryStatResListAll,
                                                                   "_given_haplotype_block")
                          fileNamesNullChrResListAll <- paste0(fileNamesNullChrResListAll,
                                                               "_given_haplotype_block")
                        }
                      }
                      if (summaryParams$dropQtl) {
                        fileNamesSummaryStatResListAll <- paste0(fileNamesSummaryStatResListAll,
                                                                 "_drop_true_QTLs")
                        fileNamesNullChrResListAll <- paste0(fileNamesNullChrResListAll,
                                                             "_drop_true_QTLs")
                      }
                      fileNamesSummaryStatResListAll <- paste0(fileNamesSummaryStatResListAll, ".rds")
                      fileNamesNullChrResListAll <- paste0(fileNamesNullChrResListAll, ".rds")


                      if (any(file.exists(fileNamesSummaryStatResListAll))) {
                        maxSaveNo <- max(which(file.exists(fileNamesSummaryStatResListAll)))
                        fileSaveNumberNow <- fileSaveNumbers[maxSaveNo]
                        fileNameSummaryStatResListAll <- fileNamesSummaryStatResListAll[maxSaveNo]
                        summaryStatResListAll <- readRDS(file = fileNameSummaryStatResListAll)
                        fileNameNullChrResListAll <- fileNamesNullChrResListAll[maxSaveNo]
                        nullChrResListAll <- readRDS(file = fileNameNullChrResListAll)
                      } else {
                        fileSaveNumberNow <- 0
                        summaryStatResListAll <- list()
                        nullChrResListAll <- list()
                      }

                      if (fileSaveNumberNow + 1 <= nTraitsHere) {
                        gbInfoUsedList <- rep(list(NULL), length(genBackCombAll))
                        names(gbInfoUsedList) <- names(genBackCombAll)
                        for (traitNo in (fileSaveNumberNow + 1):nTraitsHere) {
                          # traitNo <- 1
                          # traitNo <- 3
                          startTrait <- Sys.time()
                          traitName <- traitNames[traitNo]

                          dirIterSave <- paste0(dirTrial, traitName, "/")
                          if (!dir.exists(dirIterSave)) {
                            dir.create(dirIterSave)
                          }
                          dirIterSummary <- paste0(dirIterSave, scriptID, "_summarized_results/")
                          if (!dir.exists(dirIterSummary)) {
                            dir.create(dirIterSummary)
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

                          qtlNamesAll <- unlist(qtlNamesListAll, recursive = TRUE)
                          qtlNamesList <- unlist(qtlNamesListAll, recursive = FALSE)
                          qtlNamesList <- qtlNamesList[!sapply(X = qtlNamesList, FUN = is.null)]

                          # qtlCandsListQtls <- lapply(X = qtlNamesListQtls,
                          #                            FUN = function(x) {
                          #                              if (!is.null(x)) {
                          #                                qtlCands <- match(x, genoMap$marker)
                          #                              } else {
                          #                                qtlCands <- NULL
                          #                              }
                          #
                          #                              return(qtlCands)
                          #                            })
                          # qtlCandsListHaploQtls <- lapply(X = qtlNamesListHaploQtls,
                          #                                 FUN = function(x) {
                          #                                   if (!is.null(x)) {
                          #                                     qtlCands <- match(x, genoMap$marker)
                          #                                   } else {
                          #                                     qtlCands <- NULL
                          #                                   }
                          #
                          #                                   return(qtlCands)
                          #                                 })
                          #
                          # qtlCandsList <- lapply(X = qtlNamesList,
                          #                        FUN = function(x) {
                          #                          if (!is.null(x)) {
                          #                            qtlCands <- match(x, genoMap$marker)
                          #                          } else {
                          #                            qtlCands <- NULL
                          #                          }
                          #
                          #                          return(qtlCands)
                          #                        })

                          nQtls <- sapply(X = qtlNamesListQtls, FUN = length, simplify = TRUE)
                          nHaploQtls <- sapply(X = qtlNamesListHaploQtls, FUN = length, simplify = TRUE)

                          nQtlTotal <- sum(nQtls) + sum(nHaploQtls)
                          summaryStatResListIter <- list()
                          nullChrResListIter <- list()


                          ##### 3.3. Summary GWAS results with each genetic background and prepape for some parameters  #####
                          for (genBackEach in genBackCombAll) {
                            startGb <- Sys.time()
                            # genBackEach <- genBackCombAll[[1]]
                            # genBackEach <- genBackCombAll[[4]]
                            # genBackEach <- genBackCombAll[[7]]
                            genBackNameEach <- paste(genBackEach, collapse = "_")


                            if (length(genBackEach) == 1) {
                              methodParamListNow <- methodParamList[sapply(methodParamList,
                                                                           function(x) x$method) %in% c("Rio", "PCA", "HB", "HBxGB")]
                            } else {
                              methodParamListNow <- lapply(X = methodParamList,
                                                           FUN = function(x) {
                                                             if ((!(x$method %in% c("HB", "HBxGB"))) &
                                                                 (is.null(x$param))) {
                                                               x$param <- length(genBackEach)
                                                             }

                                                             return(x)
                                                           })
                            }

                            methodParamListPerform <- unlist(lapply(X = methodParamListNow,
                                                                    FUN = function(x) {
                                                                      x$method %in% modelMethodNow
                                                                    }))
                            methodParamListNow <- methodParamListNow[methodParamListPerform]

                            if (summaryParams$dropQtl &
                                any(sapply(X = methodParamList, FUN = function(x) x$method %in% c("HB", "HBxGB")))) {
                              addGrmSubGr <- readRDS(file = paste0("data/extra/", scriptIDPca,
                                                                   "_additive_genomic_relationship_",
                                                                   genBackNameEach, ".rds"))

                              if ((length(genBackEach) >= 2) & is.null(gbInfoUsedList[[genBackNameEach]])) {
                                fcRes <- adegenet::find.clusters(x = addGrmSubGr,
                                                                 n.pca = nPcsDapc,
                                                                 n.clust = length(genBackEach))
                                gbInfoUsedList[[genBackNameEach]] <- gbInfoUsed <- fcRes$grp
                              } else {
                                gbInfoUsed <- gbInfoUsedList[[genBackNameEach]]
                              }

                              includeEpiMethods <- sapply(X = methodParamListNow,
                                                          FUN = function(x) (x$method == "HBxGB") & is.null(x$param))
                              # includeEpiMethods <- rep(FALSE, length(methodParamListNow))

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
                            }


                            summaryStatResListGb <- list()
                            nullChrResListGb <- list()



                            ##### 3.4. Summary GWAS results #####
                            for (methodParams in methodParamListNow) {
                              # methodParams <- methodParamListNow[[1]]
                              # methodParams <- methodParamListNow[[11]]
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
                              methodHaplo <- methodParams$method %in% c("HB", "HBxGB")

                              print(paste0("Trait: ", traitName,
                                           "   Genetic background: ", genBackNameEach,
                                           "   Model: ", modelName))

                              if (!methodHaplo) {
                                fileNameResGWAS <- paste0(dirIterSave, scriptIDGWAS, "_",
                                                          "GWAS_results_for_", traitName,
                                                          "_with_", genBackNameEach,
                                                          "_by_", modelName, "_model.rds")
                              } else {
                                if (!summaryParams$useImputedBlock) {
                                  fileNameResGWAS <- paste0(dirIterSave, scriptIDGWASHB, "_",
                                                            "GWAS_results_for_", traitName,
                                                            "_with_", genBackNameEach,
                                                            "_by_", modelName, "_model.rds")
                                } else {
                                  fileNameResGWAS <- paste0(dirIterSave, scriptIDGWASHB, "_",
                                                            "GWAS_results_for_", traitName,
                                                            "_with_", genBackNameEach,
                                                            "_by_", modelName, "_model_imputed.rds")
                                }

                                if (summaryParams$dropQtl) {
                                  includeEpistasis <- (methodParams$method == "HBxGB") & is.null(methodParams$param)
                                  # includeEpistasis <- FALSE

                                  if (includeEpistasis) {
                                    addGrmSubGrNow <- addGrmSubGrEpi
                                  } else {
                                    addGrmSubGrNow <- addGrmSubGrNonEpi
                                  }
                                }
                              }
                              summaryStatResListStr <- list()
                              nullChrResListStr <- list()

                              if (file.exists(fileNameResGWAS)) {
                                resGWAS <- readRDS(file = fileNameResGWAS)
                                nTargets <- ncol(resGWAS) - 3
                                if (nTargets >= 2) {
                                  targetNames <- colnames(resGWAS)[-c(1:3)]
                                } else {
                                  targetNames <- genBackNameEach
                                }


                                for (targetNo in 1:nTargets) {
                                  # targetNo <- 1
                                  targetName <- targetNames[targetNo]
                                  dirRoc <- paste0(dirIterSummary, scriptID, "_ROC_curve/")
                                  if (!dir.exists(dirRoc)) {
                                    dir.create(dirRoc)
                                  }

                                  saveNameRoc0  <- paste0(dirRoc, scriptID, "_", traitName,
                                                          "_with_", genBackNameEach,
                                                          "_by_", modelName, "_model_",
                                                          targetName, "-targeted_")
                                  fileNameSummaryRes0 <- paste0(dirIterSummary, scriptID, "_",
                                                                "summarized_GWAS_results_for_", traitName,
                                                                "_with_", genBackNameEach,
                                                                "_by_", modelName, "_model_",
                                                                targetName,"-targeted")
                                  if (summaryParams$useGivenBlock) {
                                    if (!summaryParams$useHaploBlock) {
                                      saveNameRoc0 <- paste0(saveNameRoc0, "given_LD_block_")
                                      fileNameSummaryRes0 <- paste0(fileNameSummaryRes0, "_given_LD_block")
                                    } else {
                                      saveNameRoc0 <- paste0(saveNameRoc0, "given_haplotype_block_")
                                      fileNameSummaryRes0 <- paste0(fileNameSummaryRes0, "_given_haplotype_block")
                                    }
                                  }
                                  if (methodHaplo & summaryParams$useImputedBlock) {
                                    saveNameRoc0 <- paste0(saveNameRoc0, "imputed_")
                                    fileNameSummaryRes0 <- paste0(fileNameSummaryRes0, "_imputed")
                                  }
                                  if (summaryParams$dropQtl) {
                                    saveNameRoc0 <- paste0(saveNameRoc0, "drop_QTL_")
                                    fileNameSummaryRes0 <- paste0(fileNameSummaryRes0, "_drop_QTL")
                                  }
                                  saveNameRoc <- saveNameRoc0
                                  fileNameSummaryRes <- paste0(fileNameSummaryRes0, ".rds")


                                  if (!file.exists(fileNameSummaryRes)) {
                                    if (!methodHaplo) {
                                      geneSet <- NULL
                                      idsAroundMrcList <- idsAroundMrcListSnp
                                    } else {
                                      geneSet <- haploBlockListDf
                                      idsAroundMrcList <- idsAroundMrcListHaplo
                                    }

                                    if (summaryParams$dropQtl) {
                                      if (!methodHaplo) {
                                        dropQtlRes <- dropQtlFromData(resGWAS = resGWAS[, c(1:3, targetNo + 3)],
                                                                      genoScore = genoScore[genBackVec %in% genBackEach, , drop = FALSE],
                                                                      genoMap = genoMap,
                                                                      phenoSimInfo = phenoSimInfo,
                                                                      idsAroundMrcList = idsAroundMrcList,
                                                                      selectQtlMethod = "LD",
                                                                      lenLD = summaryParams$lenLD,
                                                                      corThres = summaryParams$corThres,
                                                                      windowSize = summaryParams$windowSize)
                                      } else {
                                        qtlHaploNamesList <- lapply(X = qtlNamesList,
                                                                    FUN = function(qtlNames) {
                                                                      sapply(X = qtlNames, FUN = defineHaploBlock,
                                                                             haploBlockListDf = geneSet,
                                                                             genoMap = genoMap,
                                                                             resGWAS = resGWAS)
                                                                    })
                                        qtlHaploNamesAll <- unlist(qtlHaploNamesList)

                                        resQtlHaplo <- performGWASHBInt(genoGWAS = genoGWAS[!(genoGWAS[, 1] %in% qtlNamesAll), ],
                                                                        phenoGWASNow = phenoGWAS[, c(1, traitNo + 1), drop = FALSE],
                                                                        haploBlockListDf = haploBlockListDf[haploBlockListDf[, 1] %in% qtlHaploNamesAll, ],
                                                                        mapHaploBlock = haploBlockMapWithCumPos[haploBlockMapWithCumPos[, 1] %in% qtlHaploNamesAll, ],
                                                                        interactionKernelOrigin = addGrmSubGr,
                                                                        addGRM = addGrm,
                                                                        addGRMSubGr = addGrmSubGrNow,
                                                                        gbInfo = genBackVec,
                                                                        gbInfoUsed = gbInfoUsed,
                                                                        gbNow = genBackEach,
                                                                        interactionMatMethod = methodParams$method,
                                                                        paramForIntGWAS = methodParams$param,
                                                                        nPcsDapc = nPcsDapc,
                                                                        nCores = summaryParams$nCores,
                                                                        adjustAddGRM = FALSE)
                                        resGWASNow <- resGWAS[, c(1:3, targetNo + 3)]
                                        resGWASNow[resGWASNow[, 1] %in% qtlHaploNamesAll, 4] <- resQtlHaplo[, 4]
                                        dropQtlRes <- dropQtlFromData(resGWAS = resGWASNow,
                                                                      genoScore = genoScore[genBackVec %in% genBackEach, , drop = FALSE],
                                                                      genoMap = genoMap,
                                                                      phenoSimInfo = phenoSimInfo,
                                                                      idsAroundMrcList = idsAroundMrcList,
                                                                      geneSet = geneSet,
                                                                      selectQtlMethod = "LD",
                                                                      lenLD = summaryParams$lenLD,
                                                                      corThres = summaryParams$corThres,
                                                                      windowSize = summaryParams$windowSize)
                                      }
                                      summaryRes <- summaryGwasResults(
                                        resGWAS = dropQtlRes$resGWAS,
                                        genoScore = dropQtlRes$genoScore,
                                        genoMap = dropQtlRes$genoMap,
                                        phenoSimInfo = dropQtlRes$phenoSimInfo,
                                        geneSet = geneSet,
                                        nTopBlockFalsePositive = nQtlTotal,
                                        sigLevel = summaryParams$sigLevel,
                                        methodThres = summaryParams$methodThres,
                                        inflatorPlus = summaryParams$inflatorPlus,
                                        lenLD = summaryParams$lenLD,
                                        corThres = summaryParams$corThres,
                                        windowSize = summaryParams$windowSize,
                                        saveName = saveNameRoc,
                                        plotROC = TRUE,
                                        returnAucRes = summaryParams$returnAucRes,
                                        returnBlockRes = summaryParams$returnBlockRes,
                                        idsAroundMrcList = dropQtlRes$idsAroundMrcList,
                                        nCores = summaryParams$nCores
                                      )
                                    } else {
                                      summaryRes <- summaryGwasResults(
                                        resGWAS = resGWAS[, c(1:3, targetNo + 3)],
                                        genoScore = genoScore[genBackVec %in% genBackEach, , drop = FALSE],
                                        genoMap = genoMap,
                                        phenoSimInfo = phenoSimInfo,
                                        geneSet = geneSet,
                                        nTopBlockFalsePositive = nQtlTotal,
                                        sigLevel = summaryParams$sigLevel,
                                        methodThres = summaryParams$methodThres,
                                        inflatorPlus = summaryParams$inflatorPlus,
                                        lenLD = summaryParams$lenLD,
                                        corThres = summaryParams$corThres,
                                        windowSize = summaryParams$windowSize,
                                        saveName = saveNameRoc,
                                        plotROC = TRUE,
                                        returnAucRes = summaryParams$returnAucRes,
                                        returnBlockRes = summaryParams$returnBlockRes,
                                        idsAroundMrcList = idsAroundMrcList,
                                        nCores = summaryParams$nCores
                                      )
                                    }
                                    saveRDS(object = summaryRes,
                                            file = fileNameSummaryRes)
                                  } else {
                                    summaryRes <- readRDS(file = fileNameSummaryRes)
                                  }


                                  summaryStatResList <- summaryRes$summaryStatResList
                                  summaryStatResList <- lapply(X = summaryStatResList,
                                                               FUN = function(x) {
                                                                 x[-length(x)]
                                                               })
                                  summaryStatResListStr[[targetName]] <- summaryStatResList

                                  nullChrResList <- summaryRes$nullChrRes
                                  nullChrResListStr[[targetName]] <- nullChrResList
                                }
                              }

                              summaryStatResListGb[[modelName]] <- summaryStatResListStr
                              nullChrResListGb[[modelName]] <- nullChrResListStr
                            }

                            summaryStatResListIter[[genBackNameEach]] <- summaryStatResListGb
                            nullChrResListIter[[genBackNameEach]] <- nullChrResListGb

                            endGb <- Sys.time()
                            compTimeGb <- endGb - startGb
                            print(paste0("Computation time: ", round(compTimeGb, 3),
                                         " ", attr(compTimeGb, "units"), " for this genetic background"))
                            cat("---------------------------------------------------------\n\n")
                          }

                          summaryStatResListAll[[traitName]] <- summaryStatResListIter
                          nullChrResListAll[[traitName]] <- nullChrResListIter

                          if (traitNo %in% fileSaveNumbers) {
                            saveRDS(object = summaryStatResListAll,
                                    file = fileNamesSummaryStatResListAll[fileSaveNumbers == traitNo])
                            saveRDS(object = nullChrResListAll,
                                    file = fileNamesNullChrResListAll[fileSaveNumbers == traitNo])
                          }


                          endTrait <- Sys.time()
                          compTimeTrait <- endTrait - startTrait
                          print(paste0("Computation time: ", round(compTimeTrait, 3),
                                       " ", attr(compTimeTrait, "units"), " for this trait"))
                          cat("=========================================================\n\n\n")
                        }
                      } else {
                        dirIterSave <- paste0(dirTrial, traitNames[fileSaveNumberNow], "/")
                        if (!dir.exists(dirIterSave)) {
                          dir.create(dirIterSave)
                        }
                        fileNameIterAllRes <- paste0(dirIterSave,
                                                     "simulated_phenotype_all_results_",
                                                     traitNames[fileSaveNumberNow], ".rds")

                        phenoSimInfo <- readRDS(file = fileNameIterAllRes)

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
                        qtlNamesAll <- unlist(qtlNamesListAll, recursive = TRUE)
                        qtlNamesList <- unlist(qtlNamesListAll, recursive = FALSE)
                        qtlNamesList <- qtlNamesList[!sapply(X = qtlNamesList, FUN = is.null)]

                        # qtlCandsListQtls <- lapply(X = qtlNamesListQtls,
                        #                            FUN = function(x) {
                        #                              if (!is.null(x)) {
                        #                                qtlCands <- match(x, genoMap$marker)
                        #                              } else {
                        #                                qtlCands <- NULL
                        #                              }
                        #
                        #                              return(qtlCands)
                        #                            })
                        # qtlCandsListHaploQtls <- lapply(X = qtlNamesListHaploQtls,
                        #                                 FUN = function(x) {
                        #                                   if (!is.null(x)) {
                        #                                     qtlCands <- match(x, genoMap$marker)
                        #                                   } else {
                        #                                     qtlCands <- NULL
                        #                                   }
                        #
                        #                                   return(qtlCands)
                        #                                 })
                        #
                        # qtlCandsList <- lapply(X = qtlNamesList,
                        #                        FUN = function(x) {
                        #                          if (!is.null(x)) {
                        #                            qtlCands <- match(x, genoMap$marker)
                        #                          } else {
                        #                            qtlCands <- NULL
                        #                          }
                        #
                        #                          return(qtlCands)
                        #                        })


                        nQtls <- sapply(X = qtlNamesListQtls, FUN = length, simplify = TRUE)
                        nHaploQtls <- sapply(X = qtlNamesListHaploQtls, FUN = length, simplify = TRUE)

                        nQtlTotal <- sum(nQtls) + sum(nHaploQtls)
                      }






                      ###### 4. Plot the summarized GWAS results ######
                      ##### 4.1. Preparation of directory and candidates of plot types #####
                      dirBarplot <- paste0(dirSummary, scriptID,
                                           "_barplots_for_summary_statistics")
                      dirBoxplot <- paste0(dirSummary, scriptID,
                                           "_boxplots_for_-log10(p)_values")

                      if (summaryParams$useGivenBlock) {
                        if (!summaryParams$useHaploBlock) {
                          dirBarplot <- paste0(dirBarplot, "_given_LD_block")
                          dirBoxplot <- paste0(dirBoxplot, "_given_LD_block")
                        } else {
                          dirBarplot <- paste0(dirBarplot, "_given_haplotype_block")
                          dirBoxplot <- paste0(dirBoxplot, "_given_haplotype_block")
                        }
                      }
                      if (summaryParams$dropQtl) {
                        dirBarplot <- paste0(dirBarplot, "_drop_QTL")
                        dirBoxplot <- paste0(dirBoxplot, "_drop_QTL")
                      }

                      dirBarplot <- paste0(dirBarplot, "/")
                      dirBoxplot <- paste0(dirBoxplot, "/")

                      if (!dir.exists(dirBarplot)) {
                        dir.create(dirBarplot)
                      }
                      if (!dir.exists(dirBoxplot)) {
                        dir.create(dirBoxplot)
                      }
                      qtlTypeNames <- c("Total", names(qtlNamesList))
                      genBackNamesList <- genBackNames
                      if (summaryParams$returnAucRes & summaryParams$returnBlockRes) {
                        summaryStatNamesAll <- c("RPF-0.05", "RPF-0.01", "AUC")
                      } else if (summaryParams$returnAucRes) {
                        summaryStatNamesAll <- "AUC"
                      } else if (summaryParams$returnBlockRes) {
                        summaryStatNamesAll <- c("RPF-0.05", "RPF-0.01")
                      } else {
                        summaryStatNamesAll <- NULL
                      }

                      summaryStatNamesLogpAll <- c("unadjusted", "adjusted")
                      plotLogpTypes <- c("SNP-based", "LD-based")



                      for (qtlTypeName in qtlTypeNames) {
                        # qtlTypeName <- qtlTypeNames[2]

                        for (genBackNames in genBackNamesList) {
                          ##### 4.2. Plot barplot for summary statistics #####
                          #### 4.2.1. Take an average for summary statistics in each iteration ####
                          # genBackNames <- genBackNamesList[[1]]
                          # genBackNames <- genBackNamesList[[7]]
                          # summaryStatNames <- names(summaryStatResList[[1]])
                          if (!is.null(summaryStatNamesAll)) {
                            if (summaryParams$returnAucRes & summaryParams$returnBlockRes) {
                              summaryStatNames <- c("recall", "precision", "harmonicMean",
                                                    "auc", "aucLdBase", "aucCausal")
                            } else if (summaryParams$returnAucRes) {
                              summaryStatNames <- c("auc", "aucLdBase", "aucCausal")
                            } else if (summaryParams$returnBlockRes) {
                              summaryStatNames <- c("recall", "precision", "harmonicMean")
                            }
                            # summaryStatName <- "auc"
                            summaryStatMeanList <- lapply(X = summaryStatNames,
                                                          FUN = function(summaryStatName) {
                                                            summaryStatArray <- do.call(
                                                              what = abind::abind,
                                                              args = lapply(X = summaryStatResListAll,
                                                                            FUN = function(summaryStatResListIter) {
                                                                              # summaryStatResListIter <- summaryStatResListAll[[1]]
                                                                              summaryStatResListIterUnlist <- unlist(
                                                                                x = unlist(x = summaryStatResListIter[genBackNames], recursive = FALSE),
                                                                                recursive = FALSE
                                                                              )
                                                                              names(summaryStatResListIterUnlist) <-
                                                                                strategyRename(names(summaryStatResListIterUnlist))

                                                                              summaryStatMatNow <- do.call(what = rbind,
                                                                                                           args = lapply(X = summaryStatResListIterUnlist,
                                                                                                                         FUN = function(summaryStatResListEach) {
                                                                                                                           (summaryStatResListEach[[qtlTypeName]])[[summaryStatName]]
                                                                                                                         }))
                                                                              summaryStatOneArrayNow <- array(
                                                                                data = summaryStatMatNow,
                                                                                dim = c(dim(summaryStatMatNow), 1),
                                                                                dimnames = c(dimnames(summaryStatMatNow), list(NULL))
                                                                              )

                                                                              return(summaryStatOneArrayNow)
                                                                            })
                                                            )
                                                            summaryStatMean <- apply(X = summaryStatArray,
                                                                                     MARGIN = c(1, 2),
                                                                                     FUN = mean, na.rm = TRUE)

                                                            if (ncol(summaryStatMean) == 1) {
                                                              colnames(summaryStatMean) <- summaryStatName
                                                            } else {
                                                              colnames(summaryStatMean) <- paste0(
                                                                summaryStatName,
                                                                stringr::str_remove(
                                                                  string = colnames(summaryStatMean),
                                                                  pattern = "BH_0.0"
                                                                ))
                                                            }

                                                            return(summaryStatMean)
                                                          })
                            names(summaryStatMeanList) <- summaryStatNames


                            summaryStatMeanMat <- do.call(what = cbind,
                                                          args = summaryStatMeanList)

                            summaryStatMeanMat[is.na(summaryStatMeanMat)] <- 0

                            #### 4.2.2. Some preparations for drawing barplots with ggplot2 ####
                            if ("China_Japan_Korea" %in% genBackNames) {
                              summaryStatMeanMat <- selectRioModel(summaryStatMeanMat)
                            }

                            summaryStatMeanDf <- data.frame(
                              Strategy = rep(rownames(summaryStatMeanMat), ncol(summaryStatMeanMat)),
                              Statistic = rep(colnames(summaryStatMeanMat), each = nrow(summaryStatMeanMat)),
                              Value = c(summaryStatMeanMat)
                            )

                            if ("China_Japan" %in% genBackNames) {
                              summaryStatMeanDf$Strategy <- stringr::str_remove(
                                string = summaryStatMeanDf$Strategy,
                                pattern = "CJ."
                              )

                              textXSize <- 22
                            } else if ("Japan_Korea" %in% genBackNames) {
                              summaryStatMeanDf$Strategy <- stringr::str_remove(
                                string = summaryStatMeanDf$Strategy,
                                pattern = "JK."
                              )

                              textXSize <- 22
                            } else if ("China_Korea" %in% genBackNames) {
                              summaryStatMeanDf$Strategy <- stringr::str_remove(
                                string = summaryStatMeanDf$Strategy,
                                pattern = "CK."
                              )

                              textXSize <- 22
                            } else if ("China_Japan_Korea" %in% genBackNames) {
                              summaryStatMeanDf$Strategy <- stringr::str_remove(
                                string = summaryStatMeanDf$Strategy,
                                pattern = "CJK."
                              )

                              summaryStatMeanDf$Strategy <- stringr::str_remove(
                                string = summaryStatMeanDf$Strategy,
                                pattern = "M3."
                              )

                              summaryStatMeanDf$Strategy <- stringr::str_replace(
                                string = summaryStatMeanDf$Strategy,
                                pattern = "All",
                                replacement = "M2.All"
                              )

                              textXSize <- 16
                            } else {
                              textXSize <- 22
                            }
                            summaryStatMeanDf$Strategy <- factor(x = summaryStatMeanDf$Strategy,
                                                                 levels = unique(summaryStatMeanDf$Strategy))


                            for (summaryStatNamesNow in summaryStatNamesAll) {
                              # summaryStatNamesNow <- "RPF-0.01"
                              dirSummaryStat <- paste0(dirBarplot, summaryStatNamesNow, "/")
                              if (!dir.exists(dirSummaryStat)) {
                                dir.create(dirSummaryStat)
                              }

                              dirSummaryStatGb <- paste0(dirSummaryStat, paste(genBackNames, collapse = "&"), "/")
                              if (!dir.exists(dirSummaryStatGb)) {
                                dir.create(dirSummaryStatGb)
                              }


                              if (summaryStatNamesNow == "RPF-0.05") {
                                summaryStatNamesRPF5 <- colnames(summaryStatMeanMat)[c(1, 3, 5)]
                                summaryStatMeanDfNow <- summaryStatMeanDf[summaryStatMeanDf$Statistic %in% summaryStatNamesRPF5, ]
                                summaryStatMeanDfNow$Statistic <- factor(x = summaryStatMeanDfNow$Statistic)
                                levels(summaryStatMeanDfNow$Statistic) <- c("F-value", "Precision", "Recall")
                                colsNow <- c("red", "green", "blue")
                              } else if (summaryStatNamesNow == "RPF-0.01") {
                                summaryStatNamesRPF1 <- colnames(summaryStatMeanMat)[c(2, 4, 6)]
                                summaryStatMeanDfNow <- summaryStatMeanDf[summaryStatMeanDf$Statistic %in% summaryStatNamesRPF1, ]
                                summaryStatMeanDfNow$Statistic <- factor(x = summaryStatMeanDfNow$Statistic)
                                levels(summaryStatMeanDfNow$Statistic) <- c("F-value", "Precision", "Recall")

                                colsNow <- c("red", "green", "blue")
                              } else if (summaryStatNamesNow == "AUC") {
                                summaryStatNamesAuc <- colnames(summaryStatMeanMat)[7:9]
                                summaryStatMeanDfNow <- summaryStatMeanDf[summaryStatMeanDf$Statistic %in% summaryStatNamesAuc, ]
                                summaryStatMeanDfNow$Statistic <- factor(x = summaryStatMeanDfNow$Statistic)
                                levels(summaryStatMeanDfNow$Statistic) <- c("AUC", "AUC (block)", "AUC (causal)")

                                colsNow <- c("red", "green", "blue")
                              }


                              #### 4.2.3. Draw barplots for 6 summary statistics and save the results with ggplot2 ####
                              plotSummaryStats <- ggplot(data = summaryStatMeanDfNow,
                                                         aes(x = Strategy, y = Value, fill = Statistic)) +
                                geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
                                ylim(0, 1) +
                                scale_color_manual(values = colsNow) +
                                ggtitle(paste0(qtlTypeName, "   ",
                                               paste(genBackNames, collapse = "&"),
                                               "   ", summaryStatNamesNow)) +
                                theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
                                      axis.title = element_text(size = 22),
                                      axis.text.x = element_text(size = textXSize, face = "bold"),
                                      axis.text.y = element_text(size = 26),
                                      legend.text = element_text(size = 23, face = "bold"),
                                      legend.title = element_text(size = 18),
                                      legend.key.height = unit(31, units = "pt"))


                              fileNamePlotSummaryStats <- paste0(dirSummaryStatGb, scriptID, "_",
                                                                 "barplots_of_", summaryStatNamesNow,
                                                                 "_for_", qtlTypeName, "_QTL_effects_in_",
                                                                 paste(genBackNames, collapse = "&"))
                              if (summaryParams$useGivenBlock) {
                                if (!summaryParams$useHaploBlock) {
                                  fileNamePlotSummaryStats <- paste0(fileNamePlotSummaryStats,
                                                                     "_given_LD_block")
                                } else {
                                  fileNamePlotSummaryStats <- paste0(fileNamePlotSummaryStats,
                                                                     "_given_haplotype_block")
                                }
                              }
                              if (summaryParams$dropQtl) {
                                fileNamePlotSummaryStats <- paste0(fileNamePlotSummaryStats,
                                                                   "_drop_QTL")
                              }


                              fileNamePlotSummaryStats <- paste0(fileNamePlotSummaryStats, ".png")
                              png(fileNamePlotSummaryStats, height = 700, width = 110 * nrow(summaryStatMeanMat))
                              print(plotSummaryStats)
                              dev.off()

                            }
                          }




                          ##### 4.3. Plot boxplot for -log10(p) values #####
                          #### 4.3.1. Extract -log10(p) values information from all the results ####
                          summaryStatNamesLogp <- c("qtlLogp", "maxLogpInLd",
                                                    "qtlAdjustedLogp", "maxAdjustedLogpInLd")
                          # summaryStatNamesLogp <- c("qtlLogp", "maxLogpInLd",
                          #                           "qtlAdjustedLogp")

                          # summaryStatName <- summaryStatNamesLogp[1]
                          rangeLogpBoxplotList <- lapply(X = summaryStatNamesLogp,
                                                         FUN = function(summaryStatName) {
                                                           summaryStatArray <- do.call(
                                                             what = abind::abind,
                                                             args = lapply(X = summaryStatResListAll,
                                                                           FUN = function(summaryStatResListIter) {
                                                                             summaryStatResListIterUnlist <- unlist(
                                                                               x = unlist(x = summaryStatResListIter, recursive = FALSE),
                                                                               recursive = FALSE
                                                                             )
                                                                             names(summaryStatResListIterUnlist) <-
                                                                               strategyRename(names(summaryStatResListIterUnlist))

                                                                             summaryStatMatNow <- do.call(what = rbind,
                                                                                                          args = lapply(X = summaryStatResListIterUnlist,
                                                                                                                        FUN = function(summaryStatResListEach) {
                                                                                                                          (summaryStatResListEach[[qtlTypeName]])[[summaryStatName]]
                                                                                                                        }))
                                                                             summaryStatOneArrayNow <- array(
                                                                               data = summaryStatMatNow,
                                                                               dim = c(dim(summaryStatMatNow), 1),
                                                                               dimnames = c(dimnames(summaryStatMatNow), list(NULL))
                                                                             )

                                                                             return(summaryStatOneArrayNow)
                                                                           })
                                                           )
                                                           summaryStatMean <- apply(X = summaryStatArray,
                                                                                    MARGIN = c(1, 3),
                                                                                    FUN = mean, na.rm = TRUE)
                                                           summaryStatMean[summaryStatMean <= -30] <- NA

                                                           return(range(summaryStatMean, na.rm = TRUE))
                                                         })
                          summaryStatLogpList <- lapply(X = summaryStatNamesLogp,
                                                        FUN = function(summaryStatName) {
                                                          summaryStatArray <- do.call(
                                                            what = abind::abind,
                                                            args = lapply(X = summaryStatResListAll,
                                                                          FUN = function(summaryStatResListIter) {
                                                                            summaryStatResListIterUnlist <- unlist(
                                                                              x = unlist(x = summaryStatResListIter[genBackNames], recursive = FALSE),
                                                                              recursive = FALSE
                                                                            )
                                                                            names(summaryStatResListIterUnlist) <-
                                                                              strategyRename(names(summaryStatResListIterUnlist))

                                                                            summaryStatMatNow <- do.call(what = rbind,
                                                                                                         args = lapply(X = summaryStatResListIterUnlist,
                                                                                                                       FUN = function(summaryStatResListEach) {
                                                                                                                         (summaryStatResListEach[[qtlTypeName]])[[summaryStatName]]
                                                                                                                       }))
                                                                            summaryStatOneArrayNow <- array(
                                                                              data = summaryStatMatNow,
                                                                              dim = c(dim(summaryStatMatNow), 1),
                                                                              dimnames = c(dimnames(summaryStatMatNow), list(NULL))
                                                                            )

                                                                            return(summaryStatOneArrayNow)
                                                                          })
                                                          )
                                                          summaryStatMean <- apply(X = summaryStatArray,
                                                                                   MARGIN = c(1, 3),
                                                                                   FUN = mean, na.rm = TRUE)
                                                          summaryStatMean[summaryStatMean <= -30] <- NA

                                                          return(summaryStatMean)
                                                        })
                          names(summaryStatLogpList) <- summaryStatNamesLogp



                          #### 4.3.2. Some preparation for drawing boxplots with ggplot2 ####
                          for (summaryStatNamesLogpNow in summaryStatNamesLogpAll) {
                            # summaryStatNamesLogpNow <- "unadjusted"
                            dirSummaryStatLogp <- paste0(dirBoxplot, summaryStatNamesLogpNow, "/")
                            if (!dir.exists(dirSummaryStatLogp)) {
                              dir.create(dirSummaryStatLogp)
                            }

                            dirSummaryStatLogpGb <- paste0(dirSummaryStatLogp, paste(genBackNames, collapse = "&"), "/")
                            if (!dir.exists(dirSummaryStatLogpGb)) {
                              dir.create(dirSummaryStatLogpGb)
                            }


                            if (summaryStatNamesLogpNow == "unadjusted") {
                              summaryStatLogpListNow <- summaryStatLogpList[c(1, 2)]
                              rangeLogpBoxplot <- range(unlist(rangeLogpBoxplotList[c(1, 2)]))
                            } else if (summaryStatNamesLogpNow == "adjusted") {
                              summaryStatLogpListNow <- summaryStatLogpList[c(3, 4)]
                              rangeLogpBoxplot <- range(unlist(rangeLogpBoxplotList[c(3, 4)]))
                            }

                            if ("Dent_Flint_Admixed" %in% genBackNames) {
                              summaryStatLogpListNow <- lapply(X = summaryStatLogpListNow,
                                                               FUN = selectRioModel)
                            }

                            summaryStatLogpDf <- data.frame(
                              Strategy = rep(rownames(summaryStatLogpListNow[[1]]),
                                             ncol(summaryStatLogpListNow[[1]]) * length(summaryStatLogpListNow)),
                              Iteration = rep(rep(colnames(summaryStatLogpListNow[[1]]),
                                                  each = nrow(summaryStatLogpListNow[[1]])),
                                              length(summaryStatLogpListNow)),
                              Type = rep(names(summaryStatLogpListNow), each = prod(dim(summaryStatLogpListNow[[1]]))),
                              Value = unlist(lapply(X = summaryStatLogpListNow, FUN = c), use.names = FALSE)
                            )
                            strategyNamesRaw <- rownames(summaryStatLogpListNow[[1]])



                            if ("China_Japan" %in% genBackNames) {
                              summaryStatLogpDf$Strategy <- stringr::str_remove(
                                string = summaryStatLogpDf$Strategy,
                                pattern = "CJ."
                              )
                            } else if ("Japan_Korea" %in% genBackNames) {
                              summaryStatLogpDf$Strategy <- stringr::str_remove(
                                string = summaryStatLogpDf$Strategy,
                                pattern = "JK."
                              )
                            } else if ("China_Korea" %in% genBackNames) {
                              summaryStatLogpDf$Strategy <- stringr::str_remove(
                                string = summaryStatLogpDf$Strategy,
                                pattern = "CK."
                              )
                            } else if ("China_Japan_Korea" %in% genBackNames) {
                              summaryStatLogpDf$Strategy <- stringr::str_remove(
                                string = summaryStatLogpDf$Strategy,
                                pattern = "CJK."
                              )

                              summaryStatLogpDf$Strategy <- stringr::str_remove(
                                string = summaryStatLogpDf$Strategy,
                                pattern = "M3."
                              )

                              summaryStatLogpDf$Strategy <- stringr::str_replace(
                                string = summaryStatLogpDf$Strategy,
                                pattern = "All",
                                replacement = "M2.All"
                              )

                              textXSize <- 18
                            } else {
                              textXSize <- 23
                            }


                            summaryStatLogpDf$Strategy <- factor(x = summaryStatLogpDf$Strategy,
                                                                 levels = unique(summaryStatLogpDf$Strategy))

                            summaryStatLogpDf$Type <- factor(x = summaryStatLogpDf$Type)
                            levels(summaryStatLogpDf$Type) <- c("LD-based", "SNP-based")

                            colsLogp <- rep(NA, length(strategyNamesRaw))
                            names(colsLogp) <- strategyNamesRaw

                            colsLogp[stringr::str_detect(string = strategyNamesRaw,
                                                         pattern = "DA")] <- "blue"
                              colsLogp[stringr::str_detect(string = strategyNamesRaw,
                                                           pattern = "M")] <- "red"
                                colsLogp[stringr::str_detect(string = strategyNamesRaw,
                                                             pattern = "GR")] <- "green"
                                  colsLogp[stringr::str_detect(string = strategyNamesRaw,
                                                               pattern = "PC")] <- "purple"
                                    colsLogp[stringr::str_detect(string = strategyNamesRaw,
                                                                 pattern = "NI")] <- "lightgreen"
                                      colsLogp[stringr::str_detect(string = strategyNamesRaw,
                                                                   pattern = "HB")] <- "pink"
                                        colsLogp[stringr::str_detect(string = strategyNamesRaw,
                                                                     pattern = "HBxGB")] <- "brown"



                                        #### 4.3.3. Draw boxplots for -log10(p) values with ggplot2 ####
                                        for (plotLogpTypeNow in plotLogpTypes) {
                                          # plotLogpTypeNow <- "SNP-based"
                                          plotLogp <- ggplot(data = summaryStatLogpDf[summaryStatLogpDf$Type %in% plotLogpTypeNow, ],
                                                             aes(x = Strategy, y = Value)) +
                                            ylab(expression(-log[10](italic(p)))) +
                                            ylim(rangeLogpBoxplot[1], rangeLogpBoxplot[2]) +
                                            geom_boxplot(fill = colsLogp) +
                                            ggtitle(paste0(qtlTypeName, "   ",
                                                           paste(genBackNames, collapse = "&"),
                                                           "   ", summaryStatNamesLogpNow, "  ", plotLogpTypeNow)) +
                                            theme(plot.title = element_text(size = 25, hjust = 0.5, family = "Verdana",
                                                                            face = "bold"),
                                                  axis.title.x = element_blank(),
                                                  axis.title.y = element_text(size = 30),
                                                  axis.text.x = element_text(size = textXSize, face = "bold"),
                                                  axis.text.y = element_text(size = 29))

                                          fileNamePlotLogp <- paste0(dirSummaryStatLogpGb, scriptID, "_",
                                                                     "boxplots_of_", summaryStatNamesLogpNow,
                                                                     "_", plotLogpTypeNow, "_-log10(p)_values_for_",
                                                                     qtlTypeName, "_QTL_effects_in_",
                                                                     paste(genBackNames, collapse = "&"))
                                          if (summaryParams$useGivenBlock) {
                                            if (!summaryParams$useHaploBlock) {
                                              fileNamePlotLogp <- paste0(fileNamePlotLogp,
                                                                         "_given_LD_block")
                                            } else {
                                              fileNamePlotLogp <- paste0(fileNamePlotLogp,
                                                                         "_given_haplotype_block")
                                            }
                                          }
                                          if (summaryParams$dropQtl) {
                                            fileNamePlotLogp <- paste0(fileNamePlotLogp,
                                                                       "_drop_QTL")
                                          }


                                          fileNamePlotLogp <- paste0(fileNamePlotLogp, ".png")

                                          png(fileNamePlotLogp, height = 1000,
                                              width = 120 * nrow(summaryStatLogpListNow[[1]]))
                                          print(plotLogp)
                                          dev.off()
                                        }
                          }







                          ##### 4.4. Boxplots for recall without using block information #####
                          recallArray <- do.call(
                            what = abind::abind,
                            args = lapply(X = summaryStatResListAll,
                                          FUN = function(summaryStatResListIter) {
                                            summaryStatResListIterUnlist <- unlist(
                                              x = unlist(x = summaryStatResListIter[genBackNames], recursive = FALSE),
                                              recursive = FALSE
                                            )
                                            names(summaryStatResListIterUnlist) <-
                                              strategyRename(names(summaryStatResListIterUnlist))

                                            summaryStatMatNow <- do.call(what = rbind,
                                                                         args = lapply(X = summaryStatResListIterUnlist,
                                                                                       FUN = function(summaryStatResListEach) {
                                                                                         (summaryStatResListEach[[qtlTypeName]])[["recallMaxLogpInLd"]]
                                                                                       }))
                                            summaryStatOneArrayNow <- array(
                                              data = summaryStatMatNow,
                                              dim = c(dim(summaryStatMatNow), 1),
                                              dimnames = c(dimnames(summaryStatMatNow), list(NULL))
                                            )

                                            return(summaryStatOneArrayNow)
                                          })
                          )
                          recallMeanMat <- apply(X = recallArray,
                                                 MARGIN = c(1, 2),
                                                 FUN = mean, na.rm = TRUE)


                          if ("China_Japan_Korea" %in% genBackNames) {
                            recallMeanMat <- selectRioModel(recallMeanMat)
                          }


                          recallMeanDf <- data.frame(
                            Strategy = rep(rownames(recallMeanMat), ncol(recallMeanMat)),
                            Threshold = rep(colnames(recallMeanMat), each = nrow(recallMeanMat)),
                            Value = c(recallMeanMat)
                          )

                          if ("China_Japan" %in% genBackNames) {
                            recallMeanDf$Strategy <- stringr::str_remove(
                              string = recallMeanDf$Strategy,
                              pattern = "CJ."
                            )

                            textXSize <- 22
                          } else if ("Japan_Korea" %in% genBackNames) {
                            recallMeanDf$Strategy <- stringr::str_remove(
                              string = recallMeanDf$Strategy,
                              pattern = "JK."
                            )

                            textXSize <- 22
                          } else if ("China_Korea" %in% genBackNames) {
                            recallMeanDf$Strategy <- stringr::str_remove(
                              string = recallMeanDf$Strategy,
                              pattern = "CK."
                            )

                            textXSize <- 22
                          } else if ("China_Japan_Korea" %in% genBackNames) {
                            recallMeanDf$Strategy <- stringr::str_remove(
                              string = recallMeanDf$Strategy,
                              pattern = "CJK."
                            )

                            recallMeanDf$Strategy <- stringr::str_remove(
                              string = recallMeanDf$Strategy,
                              pattern = "M3."
                            )

                            recallMeanDf$Strategy <- stringr::str_replace(
                              string = recallMeanDf$Strategy,
                              pattern = "All",
                              replacement = "M2.All"
                            )

                            textXSize <- 16
                          } else {
                            textXSize <- 22
                          }

                          recallMeanDf$Strategy <- factor(x = recallMeanDf$Strategy,
                                                          levels = unique(recallMeanDf$Strategy))
                          recallMeanDf$Threshold <- factor(x = recallMeanDf$Threshold,
                                                           levels = unique(recallMeanDf$Threshold))

                          colsRecall <- c("red", "blue")


                          dirRecall <- paste0(dirBarplot, "recall/")
                          if (!dir.exists(dirRecall)) {
                            dir.create(dirRecall)
                          }

                          dirRecallGb <- paste0(dirRecall, paste(genBackNames, collapse = "&"), "/")
                          if (!dir.exists(dirRecallGb)) {
                            dir.create(dirRecallGb)
                          }

                          plotRecalls <- ggplot(data = recallMeanDf,
                                                aes(x = Strategy, y = Value, fill = Threshold)) +
                            geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
                            ylim(0, 1) +
                            scale_color_manual(values = colsRecall) +
                            ggtitle(paste0(qtlTypeName, "   ",
                                           paste(genBackNames, collapse = "&"),
                                           "   Recall")) +
                            theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
                                  axis.title = element_text(size = 22),
                                  axis.text.x = element_text(size = textXSize, face = "bold"),
                                  axis.text.y = element_text(size = 26),
                                  legend.text = element_text(size = 23, face = "bold"),
                                  legend.title = element_text(size = 18),
                                  legend.key.height = unit(31, units = "pt"))



                          fileNamePlotRecalls <- paste0(dirRecallGb, scriptID, "_",
                                                        "barplots_of_recall",
                                                        "_for_", qtlTypeName, "_QTL_effects_in_",
                                                        paste(genBackNames, collapse = "&"))
                          if (summaryParams$useGivenBlock) {
                            if (!summaryParams$useHaploBlock) {
                              fileNamePlotRecalls <- paste0(fileNamePlotRecalls,
                                                            "_given_LD_block")
                            } else {
                              fileNamePlotRecalls <- paste0(fileNamePlotRecalls,
                                                            "_given_haplotype_block")
                            }
                          }
                          if (summaryParams$dropQtl) {
                            fileNamePlotRecalls <- paste0(fileNamePlotRecalls,
                                                          "_drop_QTL")
                          }

                          fileNamePlotRecalls <- paste0(fileNamePlotRecalls, ".png")
                          png(fileNamePlotRecalls, height = 700, width = 110 * nrow(recallMeanMat))
                          print(plotRecalls)
                          dev.off()
                        }
                      }


                      dirQq <- paste0(dirSummary, scriptID,
                                      "_qq_plots_for_chromosome_under_null_hypothesis/")
                      if (!dir.exists(dirQq)) {
                        dir.create(dirQq)
                      }
                      plotThresholds <- c(0.05, 0.01)

                      nullChrNos <- unlist(lapply(X = nullChrResListAll,
                                                  FUN = function(x) (((x[[1]])[[1]])[[1]])$nullChrNo))

                      if (any(nullChrNos != 0)) {

                        nullChrResListAll <- nullChrResListAll[nullChrNos != 0]
                        summaryStatResListAll <- summaryStatResListAll[nullChrNos != 0]

                        fdrNullChrRange <- range(apply(X = do.call(
                          what = abind::abind,
                          args = lapply(X = nullChrResListAll,
                                        FUN = function(nullChrResListIter) {
                                          nullChrResListIterUnlist <- unlist(
                                            x = unlist(x = nullChrResListIter, recursive = FALSE),
                                            recursive = FALSE
                                          )
                                          names(nullChrResListIterUnlist) <-
                                            strategyRename(names(nullChrResListIterUnlist))

                                          nullChrMatNow <- do.call(what = rbind,
                                                                   args = lapply(X = nullChrResListIterUnlist,
                                                                                 FUN = function(nullChrResListEach) {
                                                                                   nullChrResListEach[["fdrNullChrs"]]
                                                                                 }))
                                          nullChrOneArrayNow <- array(
                                            data = nullChrMatNow,
                                            dim = c(dim(nullChrMatNow), 1),
                                            dimnames = c(dimnames(nullChrMatNow), list(NULL))
                                          )

                                          return(nullChrOneArrayNow)
                                        })
                        ), MARGIN = c(1, 2), mean, na.rm = TRUE))



                        rangeProp <- range(do.call(
                          what = abind::abind,
                          args = lapply(X = nullChrResListAll,
                                        FUN = function(nullChrResListIter) {
                                          nullChrResListIterUnlist <- unlist(
                                            x = unlist(x = nullChrResListIter, recursive = FALSE),
                                            recursive = FALSE
                                          )
                                          names(nullChrResListIterUnlist) <-
                                            strategyRename(names(nullChrResListIterUnlist))

                                          nullChrMatNow <- do.call(what = rbind,
                                                                   args = lapply(X = nullChrResListIterUnlist,
                                                                                 FUN = function(nullChrResListEach) {
                                                                                   nullChrResListEach[["propMarkersUnderSigLevels"]]
                                                                                 }))
                                          nullChrOneArrayNow <- array(
                                            data = nullChrMatNow,
                                            dim = c(dim(nullChrMatNow), 1),
                                            dimnames = c(dimnames(nullChrMatNow), list(NULL))
                                          )

                                          return(nullChrOneArrayNow)
                                        })
                        ), na.rm = TRUE)


                        for (genBackNames in genBackNamesList) {
                          fdrNullChrArray <- do.call(
                            what = abind::abind,
                            args = lapply(X = nullChrResListAll,
                                          FUN = function(nullChrResListIter) {
                                            nullChrResListIterUnlist <- unlist(
                                              x = unlist(x = nullChrResListIter[genBackNames], recursive = FALSE),
                                              recursive = FALSE
                                            )
                                            names(nullChrResListIterUnlist) <-
                                              strategyRename(names(nullChrResListIterUnlist))

                                            nullChrMatNow <- do.call(what = rbind,
                                                                     args = lapply(X = nullChrResListIterUnlist,
                                                                                   FUN = function(nullChrResListEach) {
                                                                                     nullChrResListEach[["fdrNullChrs"]]
                                                                                   }))
                                            nullChrOneArrayNow <- array(
                                              data = nullChrMatNow,
                                              dim = c(dim(nullChrMatNow), 1),
                                              dimnames = c(dimnames(nullChrMatNow), list(NULL))
                                            )

                                            return(nullChrOneArrayNow)
                                          })
                          )
                          fdrNullChrMeanMat <- apply(X = fdrNullChrArray,
                                                     MARGIN = c(1, 2),
                                                     FUN = mean, na.rm = TRUE)


                          fdrNullChrMeanDf <- data.frame(
                            Strategy = rep(rownames(fdrNullChrMeanMat), ncol(fdrNullChrMeanMat)),
                            Threshold = rep(colnames(fdrNullChrMeanMat), each = nrow(fdrNullChrMeanMat)),
                            Value = c(fdrNullChrMeanMat)
                          )

                          if ("China_Japan" %in% genBackNames) {
                            fdrNullChrMeanDf$Strategy <- stringr::str_remove(
                              string = fdrNullChrMeanDf$Strategy,
                              pattern = "CJ."
                            )

                            textXSize <- 22
                          } else if ("Japan_Korea" %in% genBackNames) {
                            fdrNullChrMeanDf$Strategy <- stringr::str_remove(
                              string = fdrNullChrMeanDf$Strategy,
                              pattern = "JK."
                            )

                            textXSize <- 22
                          } else if ("China_Korea" %in% genBackNames) {
                            fdrNullChrMeanDf$Strategy <- stringr::str_remove(
                              string = fdrNullChrMeanDf$Strategy,
                              pattern = "CK."
                            )

                            textXSize <- 22
                          } else if ("China_Japan_Korea" %in% genBackNames) {
                            fdrNullChrMeanDf$Strategy <- stringr::str_remove(
                              string = fdrNullChrMeanDf$Strategy,
                              pattern = "CJK."
                            )

                            fdrNullChrMeanDf$Strategy <- stringr::str_remove(
                              string = fdrNullChrMeanDf$Strategy,
                              pattern = "M3."
                            )

                            fdrNullChrMeanDf$Strategy <- stringr::str_replace(
                              string = fdrNullChrMeanDf$Strategy,
                              pattern = "All",
                              replacement = "M2.All"
                            )

                            textXSize <- 16
                          } else {
                            textXSize <- 22
                          }
                          fdrNullChrMeanDf$Strategy <- factor(x = fdrNullChrMeanDf$Strategy,
                                                              levels = unique(fdrNullChrMeanDf$Strategy))
                          fdrNullChrMeanDf$Threshold <- factor(x = fdrNullChrMeanDf$Threshold,
                                                               levels = unique(fdrNullChrMeanDf$Threshold))

                          colsFdr <- c("red", "blue")


                          dirFdr <- paste0(dirBarplot, "fdr_for_chromosome_under_null_hypothesis/")
                          if (!dir.exists(dirFdr)) {
                            dir.create(dirFdr)
                          }

                          dirFdrGb <- paste0(dirFdr, paste(genBackNames, collapse = "&"), "/")
                          if (!dir.exists(dirFdrGb)) {
                            dir.create(dirFdrGb)
                          }

                          plotFdrs <- ggplot(data = fdrNullChrMeanDf,
                                             aes(x = Strategy, y = Value, fill = Threshold)) +
                            geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
                            ylim(0, fdrNullChrRange[2] + 0.05) +
                            scale_color_manual(values = colsFdr) +
                            ggtitle(paste0(paste(genBackNames, collapse = "&"),
                                           "   FDR")) +
                            theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
                                  axis.title = element_text(size = 22),
                                  axis.text.x = element_text(size = textXSize, face = "bold"),
                                  axis.text.y = element_text(size = 26),
                                  legend.text = element_text(size = 23, face = "bold"),
                                  legend.title = element_text(size = 18),
                                  legend.key.height = unit(31, units = "pt")) +
                            geom_hline(yintercept = c(0.05, 0.01),
                                       col = colsFdr,
                                       lwd = 0.5,
                                       lty = 2)



                          fileNamePlotFdrs <- paste0(dirFdrGb, scriptID, "_",
                                                     "barplots_of_fdr_for_chromosome_under_null_hypothesis",
                                                     "_in_", paste(genBackNames, collapse = "&"))

                          if (summaryParams$useGivenBlock) {
                            if (!summaryParams$useHaploBlock) {
                              fileNamePlotFdrs <- paste0(fileNamePlotFdrs,
                                                         "_given_LD_block")
                            } else {
                              fileNamePlotFdrs <- paste0(fileNamePlotFdrs,
                                                         "_given_haplotype_block")
                            }
                          }
                          if (summaryParams$dropQtl) {
                            fileNamePlotFdrs <- paste0(fileNamePlotFdrs,
                                                       "_drop_QTL")
                          }

                          fileNamePlotFdrs <- paste0(fileNamePlotFdrs, ".png")
                          png(fileNamePlotFdrs, height = 700, width = 110 * nrow(fdrNullChrMeanMat))
                          print(plotFdrs)
                          dev.off()





                          propMrkUnderThresArray <- do.call(
                            what = abind::abind,
                            args = lapply(X = nullChrResListAll,
                                          FUN = function(nullChrResListIter) {
                                            nullChrResListIterUnlist <- unlist(
                                              x = unlist(x = nullChrResListIter[genBackNames], recursive = FALSE),
                                              recursive = FALSE
                                            )
                                            names(nullChrResListIterUnlist) <-
                                              strategyRename(names(nullChrResListIterUnlist))

                                            nullChrMatNow <- do.call(what = rbind,
                                                                     args = lapply(X = nullChrResListIterUnlist,
                                                                                   FUN = function(nullChrResListEach) {
                                                                                     nullChrResListEach[["propMarkersUnderSigLevels"]]
                                                                                   }))
                                            nullChrOneArrayNow <- array(
                                              data = nullChrMatNow,
                                              dim = c(dim(nullChrMatNow), 1),
                                              dimnames = c(dimnames(nullChrMatNow), list(NULL))
                                            )

                                            return(nullChrOneArrayNow)
                                          })
                          )
                          # propMrkUnderThresArray[is.nan(propMrkUnderThresArray)] <- NA


                          if (dim(propMrkUnderThresArray)[3] >= 2) {
                            propMrkUnderThresArray <- do.call(
                              what = abind::abind,
                              args = sapply(X = 1:dim(propMrkUnderThresArray)[3],
                                            FUN = function(x){
                                              selectedMat <- propMrkUnderThresArray[, , x]

                                              selectedArray <- array(
                                                data = selectedMat,
                                                dim = c(dim(selectedMat), 1),
                                                dimnames = c(dimnames(selectedMat),
                                                             list(dimnames(propMrkUnderThresArray)[[3]][x]))
                                              )
                                            }, simplify = FALSE)
                            )
                          }
                          propMrkUnderThresArray <- propMrkUnderThresArray[
                            !apply(X = propMrkUnderThresArray,
                                   MARGIN = 1,
                                   FUN = function(x) all(is.na(x))), , , drop = FALSE
                          ]


                          propMrkUnderThresDf <- data.frame(
                            Strategy = rep(rownames(propMrkUnderThresArray),
                                           prod(dim(propMrkUnderThresArray)[2:3])),
                            Threshold = rep(rep(colnames(propMrkUnderThresArray),
                                                each = nrow(propMrkUnderThresArray)),
                                            dim(propMrkUnderThresArray)[3]),
                            Iteration = rep(dimnames(propMrkUnderThresArray)[[3]],
                                            each = prod(dim(propMrkUnderThresArray)[1:2])),
                            Value = c(propMrkUnderThresArray)
                          )






                          if ("China_Japan" %in% genBackNames) {
                            propMrkUnderThresDf$Strategy <- stringr::str_remove(
                              string = propMrkUnderThresDf$Strategy,
                              pattern = "CJ."
                            )
                          } else if ("Japan_Korea" %in% genBackNames) {
                            propMrkUnderThresDf$Strategy <- stringr::str_remove(
                              string = propMrkUnderThresDf$Strategy,
                              pattern = "JK."
                            )
                          } else if ("China_Korea" %in% genBackNames) {
                            propMrkUnderThresDf$Strategy <- stringr::str_remove(
                              string = propMrkUnderThresDf$Strategy,
                              pattern = "CK."
                            )
                          } else if ("China_Japan_Korea" %in% genBackNames) {
                            propMrkUnderThresDf$Strategy <- stringr::str_remove(
                              string = propMrkUnderThresDf$Strategy,
                              pattern = "CJK."
                            )

                            propMrkUnderThresDf$Strategy <- stringr::str_remove(
                              string = propMrkUnderThresDf$Strategy,
                              pattern = "M3."
                            )

                            propMrkUnderThresDf$Strategy <- stringr::str_replace(
                              string = propMrkUnderThresDf$Strategy,
                              pattern = "All",
                              replacement = "M2.All"
                            )

                            textXSize <- 18
                          } else {
                            textXSize <- 23
                          }

                          propMrkUnderThresDf$Strategy <- factor(x = propMrkUnderThresDf$Strategy,
                                                                 levels = unique(propMrkUnderThresDf$Strategy))
                          propMrkUnderThresDf$Threshold <- factor(x = propMrkUnderThresDf$Threshold,
                                                                  levels = unique(propMrkUnderThresDf$Threshold))
                          levels(propMrkUnderThresDf$Threshold) <- plotThresholds
                          strategyNamesRaw <- dimnames(propMrkUnderThresArray)[[1]]
                          colsProp <- rep(NA, length(strategyNamesRaw))
                          names(colsProp) <- strategyNamesRaw

                          colsProp[stringr::str_detect(string = strategyNamesRaw,
                                                       pattern = "DA")] <- "blue"
                            colsProp[stringr::str_detect(string = strategyNamesRaw,
                                                         pattern = "M")] <- "red"
                              colsProp[stringr::str_detect(string = strategyNamesRaw,
                                                           pattern = "GR")] <- "green"
                                colsProp[stringr::str_detect(string = strategyNamesRaw,
                                                             pattern = "PC")] <- "purple"
                                  colsProp[stringr::str_detect(string = strategyNamesRaw,
                                                               pattern = "NI")] <- "lightgreen"
                                    colsProp[stringr::str_detect(string = strategyNamesRaw,
                                                                 pattern = "HB")] <- "pink"
                                      colsProp[stringr::str_detect(string = strategyNamesRaw,
                                                                   pattern = "HBxGB")] <- "brown"



                                      #### 4.3.3. Draw boxplots for -log10(p) values with ggplot2 ####
                                      for (plotThresNow in plotThresholds) {
                                        # plotThresNow <- 0.05
                                        plotProp <- ggplot(data = propMrkUnderThresDf[propMrkUnderThresDf$Threshold %in% plotThresNow, ],
                                                           aes(x = Strategy, y = Value)) +
                                          ylab(paste0("Proportion of markers with p <= ", plotThresNow)) +
                                          ylim(0, rangeProp[2]) +
                                          geom_boxplot(fill = colsProp) +
                                          ggtitle(paste0(paste(genBackNames, collapse = "&"),
                                                         "   Proportion of markers with p <= ", plotThresNow)) +
                                          theme(plot.title = element_text(size = 25, hjust = 0.5, family = "Verdana",
                                                                          face = "bold"),
                                                axis.title.x = element_blank(),
                                                axis.title.y = element_text(size = 30),
                                                axis.text.x = element_text(size = textXSize, face = "bold"),
                                                axis.text.y = element_text(size = 29)) +
                                          geom_hline(yintercept = plotThresNow,
                                                     col = colsFdr[plotThresholds %in% plotThresNow],
                                                     lty = 2, lwd = 0.5)

                                        dirProp <- paste0(dirBoxplot, "proportion_of_markers_with_p_values_under_thresholds/")
                                        if (!dir.exists(dirProp)) {
                                          dir.create(dirProp)
                                        }

                                        dirPropGb <- paste0(dirProp, paste(genBackNames, collapse = "&"), "/")
                                        if (!dir.exists(dirPropGb)) {
                                          dir.create(dirPropGb)
                                        }


                                        fileNamePlotProp <- paste0(dirPropGb, scriptID, "_",
                                                                   "boxplots_of_proportion_of_markers_with_p_values_under_",
                                                                   plotThresNow, "_in_",
                                                                   paste(genBackNames, collapse = "&"))

                                        if (summaryParams$useGivenBlock) {
                                          if (!summaryParams$useHaploBlock) {
                                            fileNamePlotProp <- paste0(fileNamePlotProp,
                                                                       "_given_LD_block")
                                          } else {
                                            fileNamePlotProp <- paste0(fileNamePlotProp,
                                                                       "_given_haplotype_block")
                                          }
                                        }
                                        if (summaryParams$dropQtl) {
                                          fileNamePlotProp <- paste0(fileNamePlotProp,
                                                                     "_drop_QTL")
                                        }

                                        fileNamePlotProp <- paste0(fileNamePlotProp, ".png")
                                        png(fileNamePlotProp, height = 1000,
                                            width = 120 * nrow(propMrkUnderThresArray))
                                        print(plotProp)
                                        dev.off()
                                      }


                                      scoresNullChrList <- lapply(
                                        X = nullChrResListAll,
                                        FUN = function(nullChrResListIter) {
                                          nullChrResListIterUnlist <- unlist(
                                            x = unlist(x = nullChrResListIter[genBackNames], recursive = FALSE),
                                            recursive = FALSE
                                          )
                                          names(nullChrResListIterUnlist) <-
                                            strategyRename(names(nullChrResListIterUnlist))

                                          nullChrMatNow <- do.call(what = rbind,
                                                                   args = lapply(X = nullChrResListIterUnlist,
                                                                                 FUN = function(nullChrResListEach) {
                                                                                   nullChrResListEach[["scoresNullChr"]]
                                                                                 }))
                                          # nullChrOneArrayNow <- array(
                                          #   data = nullChrMatNow,
                                          #   dim = c(dim(nullChrMatNow), 1),
                                          #   dimnames = c(dimnames(nullChrMatNow), list(NULL))
                                          # )

                                          return(nullChrMatNow)
                                        }
                                      )

                                      plotQqRange <- c(0, max(unlist(lapply(scoresNullChrList, max, na.rm = TRUE)),
                                                              na.rm = TRUE))
                                      for (strategyNo in 1:nrow(scoresNullChrList[[1]])) {
                                        strategyNameNow <- rownames(scoresNullChrList[[1]])[strategyNo]
                                        scoresNullChrStrategyNow <- lapply(X = scoresNullChrList,
                                                                           FUN = function(x) {
                                                                             x[strategyNo, ]
                                                                           })

                                        if ("China_Japan" %in% genBackNames) {
                                          strategyNameNow <- stringr::str_remove(
                                            string = strategyNameNow,
                                            pattern = "CJ."
                                          )
                                        } else if ("Japan_Korea" %in% genBackNames) {
                                          strategyNameNow <- stringr::str_remove(
                                            string = strategyNameNow,
                                            pattern = "JK."
                                          )
                                        } else if ("China_Korea" %in% genBackNames) {
                                          strategyNameNow <- stringr::str_remove(
                                            string = strategyNameNow,
                                            pattern = "CK."
                                          )
                                        } else if ("China_Japan_Korea" %in% genBackNames) {
                                          strategyNameNow <- stringr::str_remove(
                                            string = strategyNameNow,
                                            pattern = "CJK."
                                          )

                                          strategyNameNow <- stringr::str_remove(
                                            string = strategyNameNow,
                                            pattern = "M3."
                                          )

                                          strategyNameNow <- stringr::str_replace(
                                            string = strategyNameNow,
                                            pattern = "All",
                                            replacement = "M2.All"
                                          )
                                        }

                                        dirQqGb <- paste0(dirQq, paste(genBackNames, collapse = "&"), "/")
                                        if (!dir.exists(dirQqGb)) {
                                          dir.create(dirQqGb)
                                        }

                                        fileNameQqPlotStrategyNow <- paste0(dirQqGb, scriptID, "_",
                                                                            "qq_plots_for_chromosome_under_null_hypothesis_with_",
                                                                            strategyNameNow, "_model_in_",
                                                                            paste(genBackNames, collapse = "&"),
                                                                            ".png")

                                        png(filename = fileNameQqPlotStrategyNow, width = 700, height = 700)
                                        for (iterNo in 1:length(scoresNullChrStrategyNow)) {
                                          scores <- scoresNullChrStrategyNow[[iterNo]]
                                          remove <- which(scores == 0)
                                          if (length(remove) > 0) {
                                            x <- sort(scores[-remove], decreasing = TRUE)
                                          } else {
                                            x <- sort(scores, decreasing = TRUE)
                                          }
                                          n <- length(x)
                                          unif.p <- -log10(ppoints(n))
                                          if (iterNo == 1) {
                                            plot(unif.p,
                                                 x,
                                                 xlim = plotQqRange,
                                                 ylim = plotQqRange,
                                                 pch = 1,
                                                 xlab = "Expected -log(p)",
                                                 ylab = "Observed -log(p)",
                                                 col = "grey70",
                                                 lwd = 1.1,
                                                 type = "l")
                                          } else {
                                            points(unif.p,
                                                   x,
                                                   pch = 1,
                                                   col = "grey70",
                                                   lwd = 1.1,
                                                   type = "l")
                                          }
                                        }
                                        scoresAll <- do.call(what = c, args = scoresNullChrStrategyNow)
                                        removeAll <- which(scoresAll == 0)
                                        if (length(removeAll) > 0) {
                                          xAll <- sort(scoresAll[-removeAll], decreasing = TRUE)
                                        } else {
                                          xAll <- sort(scoresAll, decreasing = TRUE)
                                        }
                                        nAll <- length(xAll)
                                        unif.pAll <- -log10(ppoints(nAll))
                                        points(unif.pAll,
                                               xAll,
                                               pch = 1,
                                               col = "grey30",
                                               lwd = 2,
                                               type = "l")
                                        lines(c(0, plotQqRange[2]), c(0, plotQqRange[2]),
                                              lty = 2, lwd = 2, col = 2)
                                        title(main = paste0(paste0(genBackNames, collapse = "&"),
                                                            "   ", strategyNameNow),
                                              cex = 2)
                                        dev.off()
                                      }
                        }
                      }


}
