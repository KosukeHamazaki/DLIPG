#############################################################################################
######  Title: 0.2_Soybean_DLIPG_Drop_scaffolds_and_convert_chromosome_ids_to_numeric  ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                               ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo      ######
######  Date: 2022/06/21 (Created), 2022/06/21 (Last Updated)                          ######
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

scriptID <- "0.2"



##### 1.2. Setting some parameters #####
dirMidDLIPGBase <- "midstream/"

fileParamsDLIPG <- paste0(dirMidDLIPGBase, scriptID,
                          "_", project, "_all_parameters.RData")
save.image(fileParamsDLIPG)


# dirScript <- paste0(dirMidDLIPGBase, scriptID,
#                     "_selection_of_accessions/")
# if (!dir.exists(dirScript)) {
#   dir.create(dirScript)
# }


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
#### 2.1.1. Read vcf files into R ####
gastonDataRaw <- read.vcf(file = "raw_data/genotype/soysnp50k_wm82.a1_41317.vcf.gz",
                          convert.chr = FALSE)
See(gastonDataRaw)


gastonDataChr <- select.snps(x = gastonDataRaw,
                             condition = chr %in% .charSeq("Gm", 1:20))
gastonDataChr@snps$chr <- readr::parse_number(gastonDataChr@snps$chr)


write.bed.matrix(x = gastonDataChr,
                 basename = "raw_data/genotype/soysnp50k_wm82.a1_41317_chrnum")
