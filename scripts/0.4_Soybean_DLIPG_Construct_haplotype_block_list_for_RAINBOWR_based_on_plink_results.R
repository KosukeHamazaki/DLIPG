###########################################################################################################
######  Title: 0.4_Soybean_DLIPG_Construct_haplotype_block_list_for_RAINBOWR_based_on_plink_results  ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                                             ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo                    ######
######  Date: 2022/06/22 (Created), 2022/06/22 (Last Updated)                                        ######
###########################################################################################################





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
    dirResearchBase <- "/media/hamazaki/HDD2/research/"   ### for Ubunutu 1
    dirScriptBase <- "~/GitHub/research_secret/"
  } else if (stringr::str_detect(string = os, pattern = "Ubuntu 20")) {
    dirResearchBase <- "/media/hamazaki-64core/39F16BAC190DFF39/research/"     ### for Ubuntu 2
    dirScriptBase <- "~/GitHub/research_secret/"
  } else {
    stop("Which type of work-station do you use?")
  }

  setwd(paste0(dirResearchBase, cropName, "/Project/", project))
}

scriptID <- "0.4"



##### 1.2. Setting some parameters #####
dirMidDLIPGBase <- "midstream/"

fileParamsDLIPG <- paste0(dirMidDLIPGBase, scriptID,
                          "_", project, "_all_parameters.RData")
save.image(fileParamsDLIPG)


dirScript <- paste0(dirMidDLIPGBase, scriptID,
                    "_haplotype_block/")
if (!dir.exists(dirScript)) {
  dir.create(dirScript)
}


##### 1.3. Import packages #####
require(data.table)
require(gaston)
require(RAINBOWR)


.charSeq <- function(prefix = "",
                     seq,
                     suffix = "") {
  sprintf(fmt = paste0(prefix,
                       "%0", floor(log10(max(seq))) + 1, "i",
                       suffix),
          seq)
}




##### 1.4. Project options #####
options(stringAsFactors = FALSE)





###### 2. Modification of data ######
##### 2.1. Read original files into R #####
fileNameGenoMap <- "data/extra/L900_Soybean_a1_chrnum_maf0.025_ms0.05_ChKoJa_300x3_maf0.025_imputed_maf0.025_SNP_geno_map.csv"

if (!file.exists(fileNameGenoMap)) {
  gastonDataRaw <- read.vcf("raw_data/genotype/L900_Soybean_a1_chrnum_maf0.025_ms0.05_ChKoJa_300x3_maf0.025_imputed_maf0.025_SNP.vcf.gz")
  genoMap <- gastonDataRaw@snps[, c("id", "chr", "pos")]

  fwrite(x = genoMap,
         file = fileNameGenoMap,
         row.names = FALSE,
         col.names = TRUE)

  rm(gastonDataRaw); gc(reset = TRUE); gc(reset = TRUE)
} else {
  genoMap <- data.frame(fread(input = fileNameGenoMap))
}
mrkNames <- genoMap$id






##### 2.2. Convert haplotype block list for RAINBOWR and impute missing blocks with one SNP #####
#### 2.2.1. Convert haplotype block list for RAINBOWR ####
fileNameHaploBlocksDet <- "raw_data/genotype/L900_Soybean_a1_chrnum_maf0.025_ms0.05_ChKoJa_300x3_maf0.025_imputed_maf0.025_plink_minmaf0.025_haplotype_block_data.blocks.det"

haploBlockDf <- convertBlockList(
  fileNameBlocksDetPlink = fileNameHaploBlocksDet,
  map = genoMap,
  blockNamesHead = "HB_",
  imputeOneSNP = FALSE,
  insertZeros = TRUE,
  n.core = 2
)

fileNameHaploBlocksDf <- paste0("data/extra/", scriptID,
                                "_haplotype_block_data_for_RAINBOWR.csv")
fwrite(x = haploBlockDf,
       file = fileNameHaploBlocksDf)


#### 2.2.2. Impute missing haplotype blocks with one SNP ####
haploBlockImpDf <- convertBlockList(
  fileNameBlocksDetPlink = fileNameHaploBlocksDet,
  map = genoMap,
  blockNamesHead = "HBI_",
  imputeOneSNP = TRUE,
  insertZeros = TRUE,
  n.core = 2
)

idsAroundMrkListHaplo <- split(x = 1:length(mrkNames),
                               f = haploBlockImpDf$block)


fileNameHaploBlocksImpDf <- paste0("data/extra/", scriptID,
                                   "_haplotype_block_data_imputed_for_RAINBOWR.csv")
fileNameIdsAroundMrkListHaplo <- paste0("data/extra/", scriptID,
                                        "_ids_around_markers_list_for_haplotype_block_imputed.rds")
fwrite(x = haploBlockImpDf,
       file = fileNameHaploBlocksImpDf)
saveRDS(object = idsAroundMrkListHaplo,
        file = fileNameIdsAroundMrkListHaplo)




##### 2.3. Convert LD block list for RAINBOWR and impute missing blocks with one SNP #####
#### 2.3.1. Convert LD block list for RAINBOWR ####
fileNameLdBlocksDet <- "raw_data/genotype/L900_Soybean_a1_chrnum_maf0.025_ms0.05_ChKoJa_300x3_maf0.025_imputed_maf0.025_plink_minmaf0_long_LD_block_data.blocks.det"
ldBlockDf <- convertBlockList(
  fileNameBlocksDetPlink = fileNameLdBlocksDet,
  map = genoMap,
  blockNamesHead = "LD_",
  imputeOneSNP = FALSE,
  insertZeros = TRUE,
  n.core = 2
)

fileNameLdBlocksDf <- paste0("data/extra/", scriptID,
                             "_long_LD_block_data_for_RAINBOWR.csv")
fwrite(x = ldBlockDf,
       file = fileNameLdBlocksDf)



#### 2.3.2. Impute missing LD blocks with one SNP ####
ldBlockImpDf <- convertBlockList(
  fileNameBlocksDetPlink = fileNameLdBlocksDet,
  map = genoMap,
  blockNamesHead = "LDI_",
  imputeOneSNP = TRUE,
  insertZeros = TRUE,
  n.core = 2
)

idsAroundMrkListLd <- split(x = 1:length(mrkNames),
                            f = ldBlockImpDf$block)


fileNameLdBlocksImpDf <- paste0("data/extra/", scriptID,
                                "_long_LD_block_data_imputed_for_RAINBOWR.csv")
fileNameIdsAroundMrkListLd <- paste0("data/extra/", scriptID,
                                     "_ids_around_markers_list_for_long_LD_block_imputed.rds")
fwrite(x = ldBlockImpDf,
       file = fileNameLdBlocksImpDf)
saveRDS(object = idsAroundMrkListLd,
        file = fileNameIdsAroundMrkListLd)









##### 2.4. Correspondence between imputed haplotype and LD blocks #####
haploLdBlockImpDf <- data.frame(HBI = haploBlockImpDf$block,
                                LDI = ldBlockImpDf$block,
                                marker = mrkNames)


ldNamesSplitHaplo <- split(
  x = haploLdBlockImpDf$LDI,
  f = haploLdBlockImpDf$HBI
)[unique(haploLdBlockImpDf$HBI)]



ldBlockNamesRepresented <- unlist(pbmcapply::pbmclapply(
  X = ldNamesSplitHaplo,
  FUN = function(ldNamesEach) {
    uniqueLdNamesEach <- unique(ldNamesEach)
    if (length(uniqueLdNamesEach) == 1) {
      ldNameRepresent <- uniqueLdNamesEach
    } else {
      tableLdNamesEach <- table(ldNamesEach)
      if (length(unique(tableLdNamesEach)) == 1) {
        ldNameRepresent <- ldNamesEach[ceiling(length(ldNamesEach) / 2)]
      } else {
        ldNameRepresent <- names(sort(tableLdNamesEach,
                                      decreasing = TRUE)[1])
      }
    }

    return(ldNameRepresent)
  },
  mc.cores = 5
))


haploImpVsLdImp <- data.frame(HBI = unique(haploBlockImpDf$block),
                              LDI = ldBlockNamesRepresented)
rownames(haploImpVsLdImp) <- 1:nrow(haploImpVsLdImp)


haploImpInLdList <- split(x = 1:nrow(haploImpVsLdImp),
                          f = haploImpVsLdImp$LDI)[unique(haploImpVsLdImp$LDI)]



fileNameHaploLdBlockImpDf <- paste0("data/extra/", scriptID,
                                    "_haplotype_and_LD_block_data_imputed.csv")
fileNameHaploImpVsLdImp <- paste0("data/extra/", scriptID,
                                  "_correspondence_between_imputed_haplotype_and_LD_blocks.csv")
fileNameHaploImpInLdList <- paste0("data/extra/", scriptID,
                                   "_imputed_haplotype_block_list_in_LD_block.rds")
fwrite(x = haploLdBlockImpDf,
       file = fileNameHaploLdBlockImpDf)
fwrite(x = haploImpVsLdImp,
       file = fileNameHaploImpVsLdImp)
saveRDS(object = haploImpInLdList,
        file = fileNameHaploImpInLdList)




##### 2.5. Correspondence between haplotype and imputed haplotype blocks #####
haploBlockImpBlocks <- haploBlockImpDf$block
names(haploBlockImpBlocks) <- mrkNames
haploVsHaploImp <- data.frame(HB = unique(haploBlockDf$block),
                              HBI = unique(haploBlockImpBlocks[haploBlockDf$marker]))
ldBlockNamesRepresentedForHaplo <- ldBlockNamesRepresented[haploVsHaploImp$HBI]
haploInLdList <- split(
  x = 1:nrow(haploVsHaploImp),
  f = ldBlockNamesRepresentedForHaplo
)[unique(ldBlockNamesRepresentedForHaplo)]


fileNameHaploVsHaploImp <- paste0("data/extra/", scriptID,
                                  "_correspondence_between_haplotype_and_imputed_haplotype_blocks.csv")
fileNameHaploInLdList <- paste0("data/extra/", scriptID,
                                "_haplotype_block_list_in_LD_block.rds")
fwrite(x = haploVsHaploImp,
       file = fileNameHaploVsHaploImp)
saveRDS(object = haploInLdList,
        file = fileNameHaploInLdList)




##### 2.6. Correspondence between LD and imputed LD blocks #####
ldBlockImpBlocks <- ldBlockImpDf$block
names(ldBlockImpBlocks) <- mrkNames
ldVsLdImp <- data.frame(LD = unique(ldBlockDf$block),
                        LDI = unique(ldBlockImpBlocks[ldBlockDf$marker]))
fileNameLdVsLdImp <- paste0("data/extra/", scriptID,
                            "_correspondence_between_LD_and_imputed_LD_blocks.csv")
fwrite(x = ldVsLdImp,
       file = fileNameLdVsLdImp)
