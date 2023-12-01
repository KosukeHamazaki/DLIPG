####################################################################################################
######  Title: 0.1_Soybean_DLIPG_Select_accessions_based_on_origin_information_and_k-medoids  ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                                      ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo             ######
######  Date: 2022/06/20 (Created), 2022/06/21 (Last Updated)                                 ######
####################################################################################################





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

scriptID <- "0.1"



##### 1.2. Setting some parameters #####
nSplit <- 10
nProp <- 0.05

msThres <- 0.01
mafThres <- 0.05

nTopNation <- 3
nSampled <- 300

dirMidDLIPGBase <- "midstream/"

fileParamsDLIPG <- paste0(dirMidDLIPGBase, scriptID,
                          "_", project, "_all_parameters.RData")
save.image(fileParamsDLIPG)


dirScript <- paste0(dirMidDLIPGBase, scriptID,
                    "_selection_of_accessions/")
if (!dir.exists(dirScript)) {
  dir.create(dirScript)
}


##### 1.3. Import packages #####
require(RAINBOWR)
require(gaston)
require(stringr)
require(readr)
require(data.table)
require(cluster)
require(utils)

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
##### 2.1. Read required files into R #####
#### 2.1.1. Read vcf file into R ####
gastonDataRaw <- read.vcf(file = "raw_data/genotype/soysnp50k_wm82.a1_41317.vcf.gz",
                          convert.chr = FALSE)
See(gastonDataRaw)

accessionNamesRaw <- gastonDataRaw@ped$id
accessionNames <- accessionNamesRaw[str_detect(string = accessionNamesRaw, pattern = "PI")]


#### 2.1.2. Read accession information files into R ####
accessionInfoRaw <- NULL
for (i in 1:nSplit) {
  fileNameAccessionInfoNow <- paste0("raw_data/extra/Search Accessions GRIN-Global-",
                                     i, ".csv")
  accessionInfoNow <- data.frame(fread(input = fileNameAccessionInfoNow))
  accessionInfoRaw <- rbind(accessionInfoRaw, accessionInfoNow)
}
accessionInfoRaw <- accessionInfoRaw[!duplicated(accessionInfoRaw$ACCESSION), ]

rownames(accessionInfoRaw) <- unlist(lapply(
  X = str_split(string = accessionInfoRaw$ACCESSION,
                pattern = " "),
  FUN = paste, collapse = ""
))
See(accessionInfoRaw)

accessionInfo <- accessionInfoRaw[accessionNames, ]
See(accessionInfo)




##### 2.2. Extract origin data and convert it into country of origin #####
accessionOrigins <- accessionInfo$ORIGIN
table(accessionOrigins)

accessionNations <- unlist(lapply(
  X = str_split(string = accessionOrigins,
                pattern = ", "),
  FUN = function(x) {
    ret <- x[length(x)]
    if (ret %in% c("South", "North")) {
      ret <- x[(length(x) - 1)]
    }

    return(ret)
  }
))

accessionNationsTabSorted <- sort(table(accessionNations), decreasing = TRUE)
accessionNationSorted <- names(accessionNationsTabSorted)




##### 2.3. Filter some SNPs in gaston data and convert it to scored marker genotype matrix #####
gastonDataPI <- select.inds(x = gastonDataRaw,
                            condition = id %in% accessionNames)
gastonDataChr <- select.snps(x = gastonDataPI,
                             condition = chr %in% .charSeq("Gm", 1:20))
gastonDataMiss <- select.snps(x = gastonDataChr,
                              condition = NAs < (nrow(gastonDataChr@ped) * msThres))
gastonDataMaf <- select.snps(x = gastonDataMiss,
                             condition = maf > mafThres)
genoData <- as.matrix(gastonDataMaf)
genoDataSp <- as(object = genoData, Class = "sparseMatrix")

rm(gastonDataPI, gastonDataChr,
   gastonDataMiss, gastonDataMaf,
   genoData)
gc(reset = TRUE); gc(reset = TRUE)






###### 3. Select accessions using origin-based k-medoids ######
##### 3.1. Select accessions with the equal proportion for each origin from the original population #####
#### 3.1.1. Select accessions using origin-based k-medoids ####
fileNameAccessionNamesSel <- paste0("raw_data/extra/", scriptID,
                                    "_accession_names_selected_by_origin_based_kmedoids_list.rds")
accessionNamesSel <- list()
pb <- txtProgressBar(min = 1,
                     max = length(accessionNationSorted),
                     style = 3)
for (accessionNationNoEach in 1:length(accessionNationSorted)) {
  setTxtProgressBar(pb = pb,
                    value = accessionNationNoEach)
  accessionNationEach <- rev(accessionNationSorted)[accessionNationNoEach]
  whichAccession <- which(accessionNations %in% accessionNationEach)

  if (length(whichAccession) >= (1 / nProp)) {
    fileNameGenoDistEach <- paste0(dirScript, scriptID,
                                   "_distance_matrix_of_accessions_from_",
                                   accessionNationEach, ".rds")

    if (!file.exists(fileNameGenoDistEach)) {
      genoDataSpEach <- genoDataSp[whichAccession, , drop = FALSE]
      genoDistEach <- dist(x = genoDataSpEach)

      saveRDS(object = genoDistEach,
              file = fileNameGenoDistEach)
      rm(genoDataSpEach)
    } else {
      genoDistEach <- readRDS(file = fileNameGenoDistEach)
    }

    clusterResEach <- pam(x = genoDistEach, k = floor(nrow(genoDataSpEach) * nProp))
    accessionNamesSel[[accessionNationEach]] <- clusterResEach$medoids

    if (length(whichAccession) >= 500) {
      rm(genoDistEach)
      gc(reset = TRUE); gc(reset = TRUE)
      saveRDS(object = accessionNamesSel,
              file = fileNameAccessionNamesSel)
    }
  } else {
    accessionNamesSel[[accessionNationEach]] <- NULL
  }
}



#### 3.1.2. Save the results ####
accessionNamesSelRev <- accessionNamesSel[accessionNationSorted]
accessionNamesSelRev <- accessionNamesSelRev[!sapply(accessionNamesSelRev, is.null)]
accessionNamesSelVec <- unlist(accessionNamesSelRev)

accessionNationsRep <- rep(accessionNationSorted,
                           floor(accessionNationsTabSorted * nProp))

accessionNamesSelDf <- data.frame(accessionNamesSelVec, accessionNamesSelVec)
accessionInfoOriginSelDf <- data.frame(accessionNamesSelVec, accessionNationsRep)
colnames(accessionInfoOriginSelDf) <- c("Accession", "Country")
accessionInfoAllSelDf <- accessionInfo[accessionNamesSelVec, ]


fileNameAccessionNamesSelDf <- paste0("data/extra/", scriptID,
                                      "_accession_names_selected_by_origin_based_kmedoids.csv")
fileNameAccessionInfoOriginSelDf <- paste0("data/extra/", scriptID,
                                           "_accession_names_with_country_of_origin_selected_by_origin_based_kmedoids.csv")
fileNameAccessionInfoAllSelDf <- paste0("data/extra/", scriptID,
                                        "_accession_information_selected_by_origin_based_kmedoids.csv")

fwrite(x = accessionNamesSelDf,
       file = fileNameAccessionNamesSelDf,
       row.names = FALSE,
       col.names = FALSE)
fwrite(x = accessionInfoOriginSelDf,
       file = fileNameAccessionInfoOriginSelDf,
       row.names = FALSE,
       col.names = TRUE)
fwrite(x = accessionInfoAllSelDf,
       file = fileNameAccessionInfoAllSelDf,
       row.names = FALSE,
       col.names = TRUE)



#### 3.1.3. Select accessions from China, Korea, and Japan, then save the results ####
accessionNamesSelRevTopNation <- accessionNamesSelRev[1:nTopNation]
accessionNamesSelRevTopNation <- accessionNamesSelRevTopNation[
  !sapply(accessionNamesSelRevTopNation, is.null)
]
accessionNamesSelTopNationVec <- unlist(accessionNamesSelRevTopNation)

accessionNationsRepTopNation <- rep(accessionNationSorted[1:nTopNation],
                                    floor(accessionNationsTabSorted[1:nTopNation] * nProp))

accessionNamesSelTopNationDf <- data.frame(accessionNamesSelTopNationVec, accessionNamesSelTopNationVec)
accessionInfoOriginSelTopNationDf <- data.frame(accessionNamesSelTopNationVec, accessionNationsRepTopNation)
colnames(accessionInfoOriginSelTopNationDf) <- c("Accession", "Country")
accessionInfoAllSelTopNationDf <- accessionInfo[accessionNamesSelTopNationVec, ]


fileNameAccessionNamesSelTopNationDf <- paste0("data/extra/", scriptID,
                                               "_accession_names_selected_by_origin_based_kmedoids_top_countires_of_origin.csv")
fileNameAccessionInfoOriginSelTopNationDf <- paste0("data/extra/", scriptID,
                                                    "_accession_names_with_country_of_origin_selected_by_origin_based_kmedoids_top_countires_of_origin.csv")
fileNameAccessionInfoAllSelTopNationDf <- paste0("data/extra/", scriptID,
                                                 "_accession_information_selected_by_origin_based_kmedoids_top_countires_of_origin.csv")

fwrite(x = accessionNamesSelTopNationDf,
       file = fileNameAccessionNamesSelTopNationDf,
       row.names = FALSE,
       col.names = FALSE)
fwrite(x = accessionInfoOriginSelTopNationDf,
       file = fileNameAccessionInfoOriginSelTopNationDf,
       row.names = FALSE,
       col.names = TRUE)
fwrite(x = accessionInfoAllSelTopNationDf,
       file = fileNameAccessionInfoAllSelTopNationDf,
       row.names = FALSE,
       col.names = TRUE)








##### 3.2. Select accessions with the equal number for each origin from China, Korea, and Japan #####
#### 3.2.1. Select accessions using origin-based k-medoids ####
accessionNationSortedSel <- accessionNationSorted[1:nTopNation]
fileNameAccessionNamesSelEqSampled <- paste0("raw_data/extra/", scriptID,
                                             "_accession_names_selected_by_origin_based_kmedoids_equal_sampling_list.rds")
accessionNamesSelEqSampled <- list()
for (accessionNationNoEach in 1:length(accessionNationSortedSel)) {
  accessionNationEach <- accessionNationSortedSel[accessionNationNoEach]
  print(accessionNationEach)


  whichAccession <- which(accessionNations %in% accessionNationEach)

  if (length(whichAccession) >= nSampled) {
    fileNameGenoDistEach <- paste0(dirScript, scriptID,
                                   "_distance_matrix_of_accessions_from_",
                                   accessionNationEach, ".rds")

    if (!file.exists(fileNameGenoDistEach)) {
      genoDataSpEach <- genoDataSp[whichAccession, , drop = FALSE]
      genoDistEach <- dist(x = genoDataSpEach)

      saveRDS(object = genoDistEach,
              file = fileNameGenoDistEach)
      rm(genoDataSpEach)
    } else {
      genoDistEach <- readRDS(file = fileNameGenoDistEach)
    }
    clusterResEachEqSampled <- pam(x = genoDistEach,
                                   k = nSampled)
    accessionNamesSelEqSampled[[accessionNationEach]] <- clusterResEachEqSampled$medoids

    rm(genoDistEach)
    gc(reset = TRUE); gc(reset = TRUE)

    saveRDS(object = accessionNamesSelEqSampled,
            file = fileNameAccessionNamesSelEqSampled)
  } else {
    accessionNamesSelEqSampled[[accessionNationEach]] <- NULL
  }
}



#### 3.2.2. Save the results ####
accessionNamesSelEqSampledVec <- unlist(accessionNamesSelEqSampled)

accessionNationsEqSampledRep <- rep(accessionNationSorted[1:nTopNation],
                                    each = nSampled)

accessionNamesSelEqSampledDf <- data.frame(accessionNamesSelEqSampledVec,
                                           accessionNamesSelEqSampledVec)
accessionInfoOriginSelEqSampledDf <- data.frame(accessionNamesSelEqSampledVec,
                                                accessionNationsEqSampledRep)
colnames(accessionInfoOriginSelEqSampledDf) <- c("Accession", "Country")
accessionInfoAllSelEqSampledDf <- accessionInfo[accessionNamesSelEqSampledVec, ]


fileNameAccessionNamesSelEqSampledDfTxt <- paste0("data/extra/", scriptID,
                                                  "_accession_names_selected_by_origin_based_kmedoids_equal_sampling.txt")
fileNameAccessionNamesSelEqSampledDf <- paste0("data/extra/", scriptID,
                                               "_accession_names_selected_by_origin_based_kmedoids_equal_sampling.csv")
fileNameAccessionInfoOriginSelEqSampledDf <- paste0("data/extra/", scriptID,
                                                    "_accession_names_with_country_of_origin_selected_by_origin_based_kmedoids_equal_sampling.csv")
fileNameAccessionInfoAllSelEqSampledDf <- paste0("data/extra/", scriptID,
                                                 "_accession_information_selected_by_origin_based_kmedoids_equal_sampling.csv")


fwrite(x = accessionNamesSelEqSampledDf,
       file = fileNameAccessionNamesSelEqSampledDfTxt,
       sep = "\t",
       row.names = FALSE,
       col.names = FALSE)
fwrite(x = accessionNamesSelEqSampledDf,
       file = fileNameAccessionNamesSelEqSampledDf,
       row.names = FALSE,
       col.names = FALSE)
fwrite(x = accessionInfoOriginSelEqSampledDf,
       file = fileNameAccessionInfoOriginSelEqSampledDf,
       row.names = FALSE,
       col.names = TRUE)
fwrite(x = accessionInfoAllSelEqSampledDf,
       file = fileNameAccessionInfoAllSelEqSampledDf,
       row.names = FALSE,
       col.names = TRUE)
