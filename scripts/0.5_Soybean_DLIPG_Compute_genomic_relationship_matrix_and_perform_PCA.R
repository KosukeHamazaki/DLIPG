##############################################################################################
######  Title: 0.5_Soybean_DLIPG_Compute_genomic_relationship_matrix_and_perform_PCA    ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                                ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo       ######
######  Date: 2022/06/25 (Created), 2022/07/13 (Last Updated)                           ######
##############################################################################################



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
# scriptIDHaploBlock <- "0.4"
scriptID <- "0.5"





##### 1.2. Setting some parameters #####
mafThres <- 0.05
ldThres <- 0.95
colsGenBack <- c("red", "blue", "green")
names(colsGenBack) <- c("China", "Japan", "Korea")

pchPlt <- 19
alphaPlt <- 0.5
cexPlt <- 1.2
cexPltAxis <- 1.5
cexPltLab <- 2
cexPltTitle <- 2.5
cexLegend <- 1.8


dirMidDLIGBBase <- "midstream/"

dirPca <- paste0(dirMidDLIGBBase,
                 scriptID, "_PCA/")
if (!dir.exists(dirPca)) {
  dir.create(dirPca)
}







##### 1.3. Import packages #####
require(data.table)
require(MASS)
require(BGLR)
require(RAINBOWR)
require(ggplot2)
require(gaston)




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



##### 2.2. Read genotype map for all genetic backgrounds into R #####
fileNameGenoMap <- "data/extra/L900_Soybean_a1_chrnum_maf0.025_ms0.05_ChKoJa_300x3_maf0.025_imputed_maf0.025_SNP_geno_map.csv"
genoMap <- data.frame(fread(input = fileNameGenoMap))
colnames(genoMap)[1] <- "marker"



##### 2.3. Read marker genotype into R #####
fileNameGenoScore <- paste0("data/genotype/L900_Soybean_a1_chrnum_maf0.025_ms0.05_ChKoJa_300x3_maf0.025_imputed_maf0.025_SNP_geno.tsv")

genoScoreDf <- data.frame(fread(input = fileNameGenoScore,
                                check.names = FALSE),
                          check.names = FALSE)
genoScore <- t(genoScoreDf[, -c(1:2)])
colnames(genoScore) <- genoMap$marker

genBackVec <- genBackVecRaw[rownames(genoScore)]





###### 3. Compute GRM, perform PCA, and extract common markers ######
mrkNamesList <- list()
for (genBackComb in genBackCombAll) {
  ##### 3.1. Process marker genotype in each subpopulation #####
  # genBackComb <- genBackCombAll[[6]]
  genBackName <- paste(genBackComb,
                       collapse = "_")
  genBackVecNow <- genBackVec[genBackVec %in% genBackComb]
  genoScoreNow <- genoScore[genBackVec %in% genBackComb, , drop = FALSE]
  mafCutRes <- MAF.cut(x.0 = genoScoreNow,
                       map.0 = genoMap,
                       min.MAF = mafThres)
  genoScoreMaf <- mafCutRes$x
  genoMapMaf <- mafCutRes$map
  mrkNamesMaf <- genoMapMaf$marker

  gastonDataMaf <- gaston::as.bed.matrix(x = genoScoreMaf + 1)
  gastonDataMaf@snps[, c(2, 1, 4)] <- genoMapMaf

  gastonDataLd <- gaston::LD.thin(x = gastonDataMaf,
                                  threshold = ldThres)
  genoScoreLd <- as.matrix(gastonDataLd) - 1
  genoMapLd <- gastonDataLd@snps[, c(2, 1, 4)]
  colnames(genoMapLd) <- colnames(genoMapMaf)
  mrkNamesNow <- genoMapLd$marker


  mrkNamesList[[genBackName]] <- mrkNamesNow



  ##### 3.2. Compute GRM and save the image #####
  grmNow <- calcGRM(genoMat = genoScoreLd)
  grmAdjNow <- calcGRM(genoMat = genoScoreLd,
                       subpop = genBackVecNow)

  fileNameAddGrmEach <- paste0("data/extra/", scriptID, "_",
                               "additive_genomic_relationship_",
                               genBackName, ".rds")
  fileNameAddGrmEachPca <- paste0(dirPca, scriptID, "_",
                                  "additive_genomic_relationship_",
                                  genBackName, ".rds")

  saveRDS(object = grmNow,
          file = fileNameAddGrmEach)
  saveRDS(object = grmNow,
          file = fileNameAddGrmEachPca)

  fileNameAddGrmAdjEach <- paste0("data/extra/", scriptID, "_",
                                  "adjusted_additive_genomic_relationship_",
                                  genBackName, ".rds")
  fileNameAddGrmAdjEachPca <- paste0(dirPca, scriptID, "_",
                                     "adjusted_additive_genomic_relationship_",
                                     genBackName, ".rds")

  saveRDS(object = grmNow,
          file = fileNameAddGrmEach)
  saveRDS(object = grmNow,
          file = fileNameAddGrmEachPca)
  saveRDS(object = grmAdjNow,
          file = fileNameAddGrmAdjEach)
  saveRDS(object = grmAdjNow,
          file = fileNameAddGrmAdjEachPca)


  fileNameAddGrmImg <- paste0(dirPca, scriptID, "_",
                              "additive_genomic_relationship_",
                              genBackName, "_image.png")
  png(filename = fileNameAddGrmImg,
      width = 600 * length(genBackComb), height = 600 * length(genBackComb))
  image(grmNow)
  title(main = paste0("Additive GRM of ", genBackName),
        cex.main = 1.5 * length(genBackComb))
  dev.off()


  fileNameAddGrmAdjImg <- paste0(dirPca, scriptID, "_",
                                 "adjusted_additive_genomic_relationship_",
                                 genBackName, "_image.png")
  png(filename = fileNameAddGrmAdjImg,
      width = 600 * length(genBackComb), height = 600 * length(genBackComb))
  image(grmAdjNow)
  title(main = paste0("Adjusted additive GRM of ", genBackName),
        cex.main = 1.5 * length(genBackComb))
  dev.off()



  ##### 3.3. Perform PCA and save the image #####
  #### 3.3.1. For original GRM ####
  pcaRes <- prcomp(x = grmNow)

  fileNamePcaRes <- paste0(dirPca, scriptID, "_",
                           "principal_component_scores_",
                           genBackName, ".png")


  pcaVar <- pcaRes$sdev ^ 2
  pcaContrib <- pcaVar / sum(pcaVar)


  png(filename = fileNamePcaRes,
      width = 1700, height = 850)
  par(mfrow = c(1, 2),
      mar = c(6.1, 5.1, 5.1, 3.1))
  plot(pcaRes$x[, 1:2],
       col = scales::alpha(colsGenBack[genBackVecNow],
                           alphaPlt),
       xlab = paste0("PC1 (",
                     round(pcaContrib[1] * 100, 2),
                     "%)"),
       ylab = paste0("PC2 (",
                     round(pcaContrib[2] * 100, 2),
                     "%)"),
       main = paste0("PC1 vs PC2 (",
                     genBackName, ")"),
       cex = cexPlt,
       cex.axis = cexPltAxis,
       cex.lab = cexPltLab,
       cex.main = cexPltTitle,
       pch = pchPlt)
  legend("bottomleft",
         legend = genBackComb,
         col = scales::alpha(colsGenBack[genBackComb],
                             alphaPlt),
         cex = cexLegend,
         pch = pchPlt)

  plot(pcaRes$x[, 3:4],
       col = scales::alpha(colsGenBack[genBackVecNow],
                           alphaPlt),
       xlab = paste0("PC3 (",
                     round(pcaContrib[3] * 100, 2),
                     "%)"),
       ylab = paste0("PC4 (",
                     round(pcaContrib[4] * 100, 2),
                     "%)"),
       main = paste0("PC3 vs PC4 (",
                     genBackName, ")"),
       cex = cexPlt,
       cex.axis = cexPltAxis,
       cex.lab = cexPltLab,
       cex.main = cexPltTitle,
       pch = pchPlt)
  legend("bottomleft",
         legend = genBackComb,
         col = scales::alpha(colsGenBack[genBackComb],
                             alphaPlt),
         cex = cexLegend,
         pch = pchPlt)
  par(mfrow = c(1, 1),
      mar = c(5.1, 4.1, 4.1, 2.1))
  dev.off()


  #### 3.3.2. For adjusted GRM ####
  pcaResAdj <- prcomp(x = grmAdjNow)

  fileNamePcaResAdj <- paste0(dirPca, scriptID, "_",
                              "adjusted_principal_component_scores_",
                              genBackName, ".png")


  pcaVarAdj <- pcaResAdj$sdev ^ 2
  pcaContribAdj <- pcaVarAdj / sum(pcaVarAdj)


  png(filename = fileNamePcaResAdj,
      width = 1700, height = 850)
  par(mfrow = c(1, 2),
      mar = c(6.1, 5.1, 5.1, 3.1))
  plot(pcaResAdj$x[, 1:2],
       col = scales::alpha(colsGenBack[genBackVecNow],
                           alphaPlt),
       xlab = paste0("PC1 (",
                     round(pcaContribAdj[1] * 100, 2),
                     "%)"),
       ylab = paste0("PC2 (",
                     round(pcaContribAdj[2] * 100, 2),
                     "%)"),
       main = paste0("PC1 vs PC2 (",
                     genBackName, ")"),
       cex = cexPlt,
       cex.axis = cexPltAxis,
       cex.lab = cexPltLab,
       cex.main = cexPltTitle,
       pch = pchPlt)
  legend("bottomleft",
         legend = genBackComb,
         col = scales::alpha(colsGenBack[genBackComb],
                             alphaPlt),
         cex = cexLegend,
         pch = pchPlt)

  plot(pcaResAdj$x[, 3:4],
       col = scales::alpha(colsGenBack[genBackVecNow],
                           alphaPlt),
       xlab = paste0("PC3 (",
                     round(pcaContribAdj[3] * 100, 2),
                     "%)"),
       ylab = paste0("PC4 (",
                     round(pcaContribAdj[4] * 100, 2),
                     "%)"),
       main = paste0("PC3 vs PC4 (",
                     genBackName, ")"),
       cex = cexPlt,
       cex.axis = cexPltAxis,
       cex.lab = cexPltLab,
       cex.main = cexPltTitle,
       pch = pchPlt)
  legend("bottomleft",
         legend = genBackComb,
         col = scales::alpha(colsGenBack[genBackComb],
                             alphaPlt),
         cex = cexLegend,
         pch = pchPlt)
  par(mfrow = c(1, 1),
      mar = c(5.1, 4.1, 4.1, 2.1))
  dev.off()
}



##### 3.4. Extract common markers to all population #####
mrkNamesCommon <- Reduce(f = intersect,
                         x = mrkNamesList)
fileNameMrkNamesCommon <- paste0("data/extra/", scriptID, "_",
                                 "marker_names_common_to_all_population_combinations.rds")
saveRDS(object = mrkNamesCommon,
        file = fileNameMrkNamesCommon)
