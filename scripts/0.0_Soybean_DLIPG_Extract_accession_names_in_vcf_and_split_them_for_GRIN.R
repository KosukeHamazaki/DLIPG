###############################################################################################
######  Title: 0.0_Soybean_DLIPG_Extract_accession_names_in_vcf_and_split_them_for_GRIN  ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                                 ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo        ######
######  Date: 2022/06/18 (Created), 2022/06/21 (Last Updated)                            ######
###############################################################################################





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

scriptID <- "0.0"



##### 1.2. Setting some parameters #####
dirMidDLIPGBase <- "midstream/"
nSplit <- 10
nEach <- 20000 / nSplit


fileParamsDLIPG <- paste0(dirMidDLIPGBase, scriptID,
                          "_", project, "_all_parameters.RData")
save.image(fileParamsDLIPG)



##### 1.3. Import packages #####
require(RAINBOWR)
require(gaston)
require(stringr)
require(readr)
require(data.table)




##### 1.4. Project options #####
options(stringAsFactors = FALSE)






###### 2. Modification of data ######
##### 2.1. Read vcf file into R #####
gastonDataRaw <- read.vcf(file = "raw_data/genotype/soysnp50k_wm82.a1_41317.vcf.gz",
                          convert.chr = FALSE)
See(gastonDataRaw)
See(gastonDataRaw@ped)
See(gastonDataRaw@snps)



##### 2.2. Extract accession names with PI #####
accessionNamesRaw <- gastonDataRaw@ped$id
accessionNames <- accessionNamesRaw[str_detect(string = accessionNamesRaw, pattern = "PI")]

accessionNamesNo <- parse_number(x = accessionNames)



##### 2.3. Change accession names for GRIN #####
accessionNamesFinalChr <- sapply(X = accessionNames,
                                 FUN = function(accessionName) {
                                   str_sub(string = accessionName,
                                           start = str_count(accessionName),
                                           end = str_count(accessionName))
                                 })
accessionNamesFinalChrNum <- as.numeric(accessionNamesFinalChr)

accessionNamesNew <- sapply(
  X = 1:length(accessionNames),
  FUN = function(accessionNameId) {
    if (is.na(accessionNamesFinalChrNum[accessionNameId])) {
      accessionNameNew <- paste0("PI ", accessionNamesNo[accessionNameId],
                                 " ", accessionNamesFinalChr[accessionNameId])
    } else {
      accessionNameNew <- paste0("PI ", accessionNamesNo[accessionNameId])
    }

    return(accessionNameNew)
  }
)

accessionNamesNewMat <- as.matrix(accessionNamesNew)



##### 2.4. Split accession names for GRIN and save the files #####
for (i in 1:nSplit) {
  fileNameAccessionNamesNew <- paste0("raw_data/extra/", scriptID,
                                      "_accession_names_for_search_in_GRIN-",
                                      i, ".csv")
  dataStart <- (i - 1) * nEach + 1
  dataEnd <- i * nEach
  dataEnd <- min(dataEnd, nrow(accessionNamesNewMat))

  fwrite(x = data.frame(accessionNamesNewMat[dataStart:dataEnd, ]),
         file = fileNameAccessionNamesNew,
         row.names = FALSE,
         col.names = FALSE)
}
