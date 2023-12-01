#############################################################################################
######  Title: 1.4.1_Soybean_DLIPG_function_to_help_plot_GWAS_results_haplotype_based  ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                               ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo      ######
######  Date: 2022/01/25 (Created), 2022/01/25 (Last Updated)                          ######
#############################################################################################

require(stringr)


strategyRename <- function(summaryStatStrategyNames) {
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "China_Japan_Korea",
                                                   replacement = "CJK")
  summaryStatStrategyNames <- stringr::str_remove(string = summaryStatStrategyNames,
                                                  pattern = ".China_Japan_Korea")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "China_Japan",
                                                   replacement = "CJ")
  summaryStatStrategyNames <- stringr::str_remove(string = summaryStatStrategyNames,
                                                  pattern = ".China_Japan")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "Japan_Korea",
                                                   replacement = "JK")
  summaryStatStrategyNames <- stringr::str_remove(string = summaryStatStrategyNames,
                                                  pattern = ".Japan_Korea")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "China_Korea",
                                                   replacement = "CK")
  summaryStatStrategyNames <- stringr::str_remove(string = summaryStatStrategyNames,
                                                  pattern = ".China_Korea")
  summaryStatStrategyNames <- stringr::str_remove(string = summaryStatStrategyNames,
                                                  pattern = "Rio-")

  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "China",
                                                   replacement = "C")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "Japan",
                                                   replacement = "J")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "Korea",
                                                   replacement = "K")
  summaryStatStrategyNames <- stringr::str_remove(string = summaryStatStrategyNames,
                                                  pattern = ".China")
  summaryStatStrategyNames <- stringr::str_remove(string = summaryStatStrategyNames,
                                                  pattern = ".Japan")
  summaryStatStrategyNames <- stringr::str_remove(string = summaryStatStrategyNames,
                                                  pattern = ".Korea")

  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "GROUP-0",
                                                   replacement = "NI0")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "GROUP-",
                                                   replacement = "KM")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "PCA-",
                                                   replacement = "PC")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "DAPC-",
                                                   replacement = "DA")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "NI-",
                                                   replacement = "NI")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "true",
                                                   replacement = "T")

  return(summaryStatStrategyNames)
}

strategyRename2 <- function(summaryStatStrategyNames) {
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "China_Japan_Korea",
                                                   replacement = "CJK")
  summaryStatStrategyNames <- stringr::str_remove(string = summaryStatStrategyNames,
                                                  pattern = ".China_Japan_Korea")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "China_Japan",
                                                   replacement = "CJ")
  summaryStatStrategyNames <- stringr::str_remove(string = summaryStatStrategyNames,
                                                  pattern = ".China_Japan")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "Japan_Korea",
                                                   replacement = "JK")
  summaryStatStrategyNames <- stringr::str_remove(string = summaryStatStrategyNames,
                                                  pattern = ".Japan_Korea")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "China_Korea",
                                                   replacement = "CK")
  summaryStatStrategyNames <- stringr::str_remove(string = summaryStatStrategyNames,
                                                  pattern = ".China_Korea")
  summaryStatStrategyNames <- stringr::str_remove(string = summaryStatStrategyNames,
                                                  pattern = "Rio-")

  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "China",
                                                   replacement = "C")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "Japan",
                                                   replacement = "J")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "Korea",
                                                   replacement = "K")
  summaryStatStrategyNames <- stringr::str_remove(string = summaryStatStrategyNames,
                                                  pattern = ".China")
  summaryStatStrategyNames <- stringr::str_remove(string = summaryStatStrategyNames,
                                                  pattern = ".Japan")
  summaryStatStrategyNames <- stringr::str_remove(string = summaryStatStrategyNames,
                                                  pattern = ".Korea")

  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "GROUP-0",
                                                   replacement = "SNP0")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "GROUP-",
                                                   replacement = "KM")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "PCA-",
                                                   replacement = "PC")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "DAPC-",
                                                   replacement = "DA")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "NI-",
                                                   replacement = "SNP")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "true",
                                                   replacement = "T")

  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "KM3",
                                                   replacement = "KM")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "PC2",
                                                   replacement = "PC")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "DA3",
                                                   replacement = "DA")
  summaryStatStrategyNames <- stringr::str_replace(string = summaryStatStrategyNames,
                                                   pattern = "SNP3",
                                                   replacement = "SNP")

  return(summaryStatStrategyNames)
}



selectRioModel <- function(mat) {
  strategyNames <- rownames(mat)
  rioModel <- stringr::str_detect(string = strategyNames, pattern = "CJK.M3")
  notTotal <- !stringr::str_detect(string = strategyNames, pattern = "All")
  if (qtlTypeName == "Total") {
    whichRemove <- which(rioModel & notTotal)
  } else if (str_detect(string = qtlTypeName, pattern = "Common")) {
    divergent <- stringr::str_detect(string = strategyNames, pattern = "-")
    whichRemove <- which(rioModel & divergent)
  } else if (str_detect(string = qtlTypeName, pattern = "China")) {
    china <- stringr::str_count(string = strategyNames, pattern = "C") >= 2
    whichRemove <- which(rioModel & (!(china | (!notTotal))))
  } else if (str_detect(string = qtlTypeName, pattern = "Japan")) {
    japan <- stringr::str_count(string = strategyNames, pattern = "J") >= 2
    whichRemove <- which(rioModel & (!(japan | (!notTotal))))
  } else if (str_detect(string = qtlTypeName, pattern = "Korea")) {
    korea <- stringr::str_count(string = strategyNames, pattern = "K") >= 2
    whichRemove <- which(rioModel & (!(korea | (!notTotal))))
  } else if (str_detect(string = qtlTypeName, pattern = "C-J")) {
    divergent <- stringr::str_detect(string = strategyNames, pattern = "C-J")
    whichRemove <- which(rioModel & (!divergent) & notTotal)
  } else if (str_detect(string = qtlTypeName, pattern = "J-K")) {
    divergent <- stringr::str_detect(string = strategyNames, pattern = "J-K")
    whichRemove <- which(rioModel & (!divergent) & notTotal)
  } else if (str_detect(string = qtlTypeName, pattern = "K-C")) {
    divergent <- stringr::str_detect(string = strategyNames, pattern = "C-K")
    whichRemove <- which(rioModel & (!divergent) & notTotal)
  } else if (any(str_detect(string = qtlTypeName, pattern = c("PC", "Polygenes")))) {
    common <- stringr::str_detect(string = strategyNames, pattern = "\\+")
    divergent <- stringr::str_detect(string = strategyNames, pattern = "-")
    notTotal <- !stringr::str_detect(string = strategyNames, pattern = "All")
    whichRemove <- which(rioModel & !(common | divergent) & notTotal)
  }

  mat <- mat[-whichRemove, , drop = FALSE]

  return(mat)
}


selectRioModel2 <- function(mat) {
  strategyNames <- rownames(mat)

  # rioModel <- stringr::str_detect(string = strategyNames, pattern = "CJK.M3")
  # notTotal <- !stringr::str_detect(string = strategyNames, pattern = "All")
  # notChina <- strategyNames != "CJK.M3.C"
  # notCJK <- strategyNames != "CJK.M3.C+J+K"
  # whichRemove <- which(rioModel & notTotal & notChina & notCJK)
  #
  # mat <- mat[-whichRemove, , drop = FALSE]

  whichRemain <- strategyNames %in% c("CJK.M3.All", "CJK.M3.C",
                                      "CJK.M3.C+J+K", "CJK.KM",
                                      "CJK.PC", "CJK.DA", "CJK.SNP",
                                      "CJK.HB", "CJK.HBxGB")
  strategyNames <- strategyNames[whichRemain]

  mat <- mat[whichRemain, , drop = FALSE]
  return(mat)
}


selectRioModelStrategyNames2 <- function(strategyNames) {
  # rioModel <- stringr::str_detect(string = strategyNames, pattern = "CJK.M3")
  # notTotal <- !stringr::str_detect(string = strategyNames, pattern = "All")
  # notChina <- strategyNames != "CJK.M3.C"
  # notCJK <- strategyNames != "CJK.M3.C+J+K"
  # notSNPT <-
  # whichRemove <- which(rioModel & notTotal & notChina & notCJK)

  whichRemain <- strategyNames %in% c("CJK.M3.All", "CJK.M3.C",
                                      "CJK.M3.C+J+K", "CJK.KM",
                                      "CJK.PC", "CJK.DA", "CJK.SNP",
                                      "CJK.HB", "CJK.HBxGB")
  strategyNames <- strategyNames[whichRemain]

  return(strategyNames)
}
