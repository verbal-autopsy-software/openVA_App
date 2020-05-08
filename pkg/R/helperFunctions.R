#' Get the CODs assigned to each individual in the data
#'
#' @param x An fitted object returned from openVA::codeVA()
#' @param top An integer specifying the number of top causes included
#' in the results.
#'
#' @return A matrix with the death ID, top causes, and associated
#' probability for the top causes.
#'
#' @export
#'
indivCOD <- function (x, top = 3) {

  probs <- openVA::getIndivProb(x)
  cods  <- colnames(probs)

  out <- matrix(NA, nrow = nrow(probs), ncol = top*2)

  for (i in 1:nrow(probs)) {
    probsOrdered <- order(probs[i,], decreasing = TRUE)
    newTop <- top
    if(length(unique(probsOrdered)) <= top) newTop = (top - 1)
    if(newTop < 1){
      cat("Error: not enough unique causes \n")
      next
    }
    for (j in 1:newTop) {
      k <- seq(1, top*2, by = 2)[j]
      out[i, k  ] <- cods[probsOrdered[j]]
      out[i, k+1] <- round(probs[i, probsOrdered[j]], 4)
    }
  }
  
  
  out <- cbind(rownames(probs), out)
  colnames(out) <- c("ID", "Most Likely Cause", "Probability",
                     "Second Most Likely Cause", "Probability",
                     "Third Most Likely Cause", "Probability")
  out <- data.frame(out)

  # reorder output
  if(!is.null(x$ID_orig)){
    new <- data.frame(ID = x$ID_orig)
    out <- merge(new, out, all.x = TRUE, by = "ID")
    out <- out[match(x$ID_orig, out$ID), ]
    if(sum(is.na(out[,2])) > 0){
        # In the future, we may want to add something here to replace NA
    }
  }else{
    warning("No initial order specified. The output may not be in the correct order!!")
  }

  out
}

#' Format InterVA5 & InSilicoVA summary objects to write as CSV file
#'
#' @param x A summary of object returned from openVA::codeVA()
#'
#' @export
#'
csmfSummaryCSV <- function (x) {

  alg <- class(x)
  if (grepl("interVA", alg)) {
    capturedHead <- utils::capture.output(print(x))[1:5]
    capturedTail <- utils::tail(utils::capture.output(print(x)), n = 9)
    #startTail <- which(grepl('Top', capturedTail)) - 1
    top <- x$top
  } else {
    capturedHead <- utils::capture.output(print(x))[1:8]
    capturedTail <- NULL
    top <- x$showTop
  }
  numRows <- top + length(capturedHead) + 1 + length(capturedTail)
  numRows <- top + length(capturedHead) + 1 + 9
  numCols <- ifelse(grepl("interVA", alg), 2, 6)
  outMat <- matrix("", nrow = numRows, ncol = numCols)
  outMat[1:length(capturedHead), 1] <- capturedHead
  if (grepl("interVA", alg)) {
    csmfCols <- 2
    } else {
      csmfCols <- 1:5
  }
  csmf.out.ordered <- x$csmf.ordered[1:top, ]
  csmf.out.ordered[, csmfCols] <- round(csmf.out.ordered[, csmfCols], 4)
  rowStart <- length(capturedHead) + 1
  if (grepl("interVA", alg)) {
    outMat[rowStart,] <- names(csmf.out.ordered)
    rowStart <- rowStart + 1
    rowStop <- rowStart + nrow(csmf.out.ordered) - 1
    outMat[rowStart:rowStop, ] <- as.matrix(csmf.out.ordered)
  } else {
    outMat[rowStart,] <- c("cause", colnames(csmf.out.ordered))
    rowStart <- rowStart + 1
    rowStop <- rowStart + nrow(csmf.out.ordered) - 1
    outMat[rowStart:rowStop, ] <- cbind(rownames(as.matrix(csmf.out.ordered)),
                                        as.matrix(csmf.out.ordered))
  }
  if (!is.null(capturedTail)) {
    rowStart <- rowStop + 1
    rowStop <- rowStart + 2
    outMat[rowStart:rowStop, ] <- matrix(c(
                                           "", "Top 6 Circumstances of Mortality Category:",
                                           "cause", "", "", "likelihood"),
                                         nrow = 3)
    rowStart <- rowStop + 1
    rowStop <- numRows
    csmf.out.ordered <- x$comcat[1:6, ]
    csmf.out.ordered[, 2] <- round(csmf.out.ordered[, 2], 4)
    outMat[rowStart:rowStop, ] <- as.matrix(csmf.out.ordered)
  }
  return (outMat)
}

#' Separate InterVA5 results by age and sex.
#'
#' @param x A fitted object returned from openVA::codeVA()
#'
#' @export
#'
sepVAResults <- function (x) {
  
  index <- which(colnames(x$checkedData) == 'i019a')
  idMale <- x$checkedData[x$checkedData[, index] == 1, 1]
  index <- which(colnames(x$checkedData) == 'i019b')
  idFemale <- x$checkedData[x$checkedData[, index] == 1, 1]
  index <- which(colnames(x$checkedData) == 'i022g')
  idNeonate <- x$checkedData[x$checkedData[, index] == 1 & 
                             !is.na(x$checkedData[, index]), 1]
  index1 <- which(colnames(x$checkedData) == 'i022f')
  index2 <- which(colnames(x$checkedData) == 'i022e')
  index3 <- which(colnames(x$checkedData) == 'i022d')
  idChild <- x$checkedData[(x$checkedData[, index1] == 1 &
                              !is.na(x$checkedData[, index1])) |
                            (x$checkedData[, index2] == 1 &
                              !is.na(x$checkedData[, index2])) |
                            (x$checkedData[, index3] == 1 &
                              !is.na(x$checkedData[, index3])), 1]
  index1 <- which(colnames(x$checkedData) == 'i022a')
  index2 <- which(colnames(x$checkedData) == 'i022b')
  index3 <- which(colnames(x$checkedData) == 'i022c')
  idAdult <- x$checkedData[(x$checkedData[, index1] == 1 &
                            !is.na(x$checkedData[, index1])) |
                           (x$checkedData[, index2] == 1 &
                            !is.na(x$checkedData[, index2])) |
                           (x$checkedData[, index3] == 1 &
                            !is.na(x$checkedData[, index3])), 1]
  idMNeonate <- idNeonate[idNeonate %in% idMale]
  idMChild <- idChild[idChild %in% idMale]
  idMAdult <- idAdult[idAdult %in% idMale]
  idFNeonate <- idNeonate[idNeonate %in% idFemale]
  idFChild <- idChild[idChild %in% idFemale]
  idFAdult <- idAdult[idAdult %in% idFemale]
  
  results <- list(all = x)
  keep <- x$ID %in% idMale
  if (sum(keep) == 0) {
    results$male <- NULL
  } else {
    results$male <- list(ID = x$ID[keep], VA5 = x$VA5[keep])
    class(results$male) <- class(x)
  }
  keep <- x$ID %in% idFemale
  if (sum(keep) == 0) {
    results$female <- NULL
  } else {
    results$female <- list(ID = x$ID[keep], VA5 = x$VA5[keep])
    class(results$female) <- class(x)
  }
  keep <- x$ID %in% idNeonate
  if (sum(keep) == 0) {
    results$neonate <- NULL
  } else {
    results$neonate <- list(ID = x$ID[keep], VA5 = x$VA5[keep])
    class(results$neonate) <- class(x)
  }
  keep <- x$ID %in% idChild
  if (sum(keep) == 0) {
    results$child <- NULL
  } else {
    results$child <- list(ID = x$ID[keep], VA5 = x$VA5[keep])
    class(results$child) <- class(x)
  }
  keep <- x$ID %in% idAdult
  if (sum(keep) == 0) {
    results$adult <- NULL
  } else {
    results$adult <- list(ID = x$ID[keep], VA5 = x$VA5[keep])
    class(results$adult) <- class(x)
  }
  keep <- x$ID %in% idMNeonate
  if (sum(keep) == 0) {
    results$mNeonate <- NULL
  } else {
    results$mNeonate <- list(ID = x$ID[keep], VA5 = x$VA5[keep])
    class(results$mNeonate) <- class(x)
  }
  keep <- x$ID %in% idMChild
  if (sum(keep) == 0) {
    results$mChild <- NULL
  } else {
    results$mChild <- list(ID = x$ID[keep], VA5 = x$VA5[keep])
    class(results$mChild) <- class(x)
  }
  keep <- x$ID %in% idMAdult
  if (sum(keep) == 0) {
    results$mAdult <- NULL
  } else {
    results$mAdult <- list(ID = x$ID[keep], VA5 = x$VA5[keep])
    class(results$mAdult) <- class(x)
  }
  keep <- x$ID %in% idFNeonate
  if (sum(keep) == 0) {
    results$fNeonate <- NULL
  } else {
    results$fNeonate <- list(ID = x$ID[keep], VA5 = x$VA5[keep])
    class(results$fNeonate) <- class(x)
  }
  keep <- x$ID %in% idFChild
  if (sum(keep) == 0) {
    results$fChild <- NULL
  } else {
    results$fChild <- list(ID = x$ID[keep], VA5 = x$VA5[keep])
    class(results$fChild) <- class(x)
  }
  keep <- x$ID %in% idFAdult
  if (sum(keep) == 0) {
    results$fAdult <- NULL
  } else {
    results$fAdult <- list(ID = x$ID[keep], VA5 = x$VA5[keep])
    class(results$fAdult) <- class(x)
  }
  return(results)
}

#' Separate InterVA5 log messages by age and sex.
#'
#' @param x A fitted object returned from openVA::codeVA()
#' @param names_runs Character vector with names of demographic
#'   groups to include in the stratification of results.
#'   Possible values include: all, male, female, neonate, child,
#'   adult, mNeonate, mChild, mAdult, fNeonate, fChild, and fAdult.
#' @param log_file String with path to log file.
#'
#' @export
#'
sepVALog <- function (x, names_runs, log_file) {
  
  idAll <- x$checkedData[, 1]
  index <- which(colnames(x$checkedData) == "i019a")
  idMale <- x$checkedData[x$checkedData[, index] == 1, 1]
  index <- which(colnames(x$checkedData) == "i019b")
  idFemale <- x$checkedData[x$checkedData[, index] == 1, 1]
  index <- which(colnames(x$checkedData) == "i022g")
  idNeonate <- x$checkedData[x$checkedData[, index] == 1 & 
                               !is.na(x$checkedData[, index]), 1]
  index1 <- which(colnames(x$checkedData) == "i022f")
  index2 <- which(colnames(x$checkedData) == "i022e")
  index3 <- which(colnames(x$checkedData) == "i022d")
  idChild <- x$checkedData[(x$checkedData[, index1] == 1 &
                              !is.na(x$checkedData[, index1])) |
                             (x$checkedData[, index2] == 1 &
                                !is.na(x$checkedData[, index2])) |
                             (x$checkedData[, index3] == 1 &
                                !is.na(x$checkedData[, index3])), 1]
  index1 <- which(colnames(x$checkedData) == "i022a")
  index2 <- which(colnames(x$checkedData) == "i022b")
  index3 <- which(colnames(x$checkedData) == "i022c")
  idAdult <- x$checkedData[(x$checkedData[, index1] == 1 &
                              !is.na(x$checkedData[, index1])) |
                             (x$checkedData[, index2] == 1 &
                                !is.na(x$checkedData[, index2])) |
                             (x$checkedData[, index3] == 1 &
                                !is.na(x$checkedData[, index3])), 1]
  idMNeonate <- idNeonate[idNeonate %in% idMale]
  idMChild <- idChild[idChild %in% idMale]
  idMAdult <- idAdult[idAdult %in% idMale]
  idFNeonate <- idNeonate[idNeonate %in% idFemale]
  idFChild <- idChild[idChild %in% idFemale]
  idFAdult <- idAdult[idAdult %in% idFemale]
  
  outFileName <- paste0("__tmp__", names_runs, "/errorlogV5.txt")
  groupName <- gsub("^(.)", "id\\U\\1", names_runs, perl = TRUE)
  for (i in 1:length(names_runs)) {
    ids <- get(groupName[i])
    rmIDs <- idAll[!(idAll %in% ids)]
    sedIDs <- paste0(rmIDs, collapse = " |")
    sedCMD1 <- "sed -E"
    sedCMD2 <- paste0("'/^(", sedIDs, ") /d'")
    sedCMD3 <- paste0(log_file, " > ", outFileName[i])
    system(paste(sedCMD1, sedCMD2, sedCMD3))
  }
}


#' Return list of countries and abbreviations for Tariff2 input
#'
#'
#' @return A list of country names and abbreviations.
#'
#' @noRd
smartVA_countries <- function () {
  c(
    "Unknown" = "Unknown",
    "Afghanistan" = "AFG",
    "Albania" = "ALB",
    "Algeria" = "DZA",
    "Andorra" = "AND",
    "Angola" = "AGO",
    "Antigua and Barbuda" = "ATG",
    "Argentina" = "ARG",
    "Armenia" = "ARM",
    "Australia" = "AUS",
    "Austria" = "AUT",
    "Azerbaijan" = "AZE",
    "Bahrain" = "BHR",
    "Bangladesh" = "BGD",
    "Barbados" = "BRB",
    "Belarus" = "BLR",
    "Belgium" = "BEL",
    "Belize" = "BLZ",
    "Benin" = "BEN",
    "Bhutan" = "BTN",
    "Bolivia" = "BOL",
    "Bosnia and Herzegovina" = "BIH",
    "Botswana" = "BWA",
    "Brazil" = "BRA",
    "Brunei" = "BRN",
    "Bulgaria" = "BGR",
    "Burkina Faso" = "BFA",
    "Burundi" = "BDI",
    "Cambodia" = "KHM",
    "Cameroon" = "CMR",
    "Canada" = "CAN",
    "Cape Verde" = "CPV",
    "Central African Republic" = "CAF",
    "Chad" = "TCD",
    "Chile" = "CHL",
    "China" = "CHN",
    "Colombia" = "COL",
    "Comoros" = "COM",
    "Congo" = "COG",
    "Costa Rica" = "CRI",
    "Cote d'Ivoire" = "CIV",
    "Croatia" = "HRV",
    "Cuba" = "CUB",
    "Cyprus" = "CYP",
    "Czech Republic" = "CZE",
    "Democratic Republic of the Congo" = "COD",
    "Denmark" = "DNK",
    "Djibouti" = "DJI",
    "Dominica" = "DMA",
    "Dominican Republic" = "DOM",
    "Ecuador" = "ECU",
    "Egypt" = "EGY",
    "El Salvador" = "SLV",
    "Equatorial Guinea" = "GNQ",
    "Eritrea" = "ERI",
    "Estonia" = "EST",
    "Ethiopia" = "ETH",
    "Federated States of Micronesia" = "FSM",
    "Fiji" = "FJI",
    "Finland" = "FIN",
    "France" = "FRA",
    "Gabon" = "GAB",
    "Georgia" = "GEO",
    "Germany" = "DEU",
    "Ghana" = "GHA",
    "Greece" = "GRC",
    "Grenada" = "GRD",
    "Guatemala" = "GTM",
    "Guinea" = "GIN",
    "Guinea-Bissau" = "GNB",
    "Guyana" = "GUY",
    "Haiti" = "HTI",
    "Honduras" = "HND",
    "Hungary" = "HUN",
    "Iceland" = "ISL",
    "India" = "IND",
    "Indonesia" = "IDN",
    "Iran" = "IRN",
    "Iraq" = "IRQ",
    "Ireland" = "IRL",
    "Israel" = "ISR",
    "Italy" = "ITA",
    "Jamaica" = "JAM",
    "Japan" = "JPN",
    "Jordan" = "JOR",
    "Kazakhstan" = "KAZ",
    "Kenya" = "KEN",
    "Kiribati" = "KIR",
    "Kuwait" = "KWT",
    "Kyrgyzstan" = "KGZ",
    "Laos" = "LAO",
    "Latvia" = "LVA",
    "Lebanon" = "LBN",
    "Lesotho" = "LSO",
    "Liberia" = "LBR",
    "Libya" = "LBY",
    "Lithuania" = "LTU",
    "Luxembourg" = "LUX",
    "Macedonia" = "MKD",
    "Madagascar" = "MDG",
    "Malawi" = "MWI",
    "Malaysia" = "MYS",
    "Maldives" = "MDV",
    "Mali" = "MLI",
    "Malta" = "MLT",
    "Marshall Islands" = "MHL",
    "Mauritania" = "MRT",
    "Mauritius" = "MUS",
    "Mexico" = "MEX",
    "Moldova" = "MDA",
    "Mongolia" = "MNG",
    "Montenegro" = "MNE",
    "Morocco" = "MAR",
    "Mozambique" = "MOZ",
    "Myanmar" = "MMR",
    "Namibia" = "NAM",
    "Nepal" = "NPL",
    "Netherlands" = "NLD",
    "New Zealand" = "NZL",
    "Nicaragua" = "NIC",
    "Niger" = "NER",
    "Nigeria" = "NGA",
    "North Korea" = "PRK",
    "Norway" = "NOR",
    "Oman" = "OMN",
    "Pakistan" = "PAK",
    "Palestine" = "PSE",
    "Panama" = "PAN",
    "Papua New Guinea" = "PNG",
    "Paraguay" = "PRY",
    "Peru" = "PER",
    "Philippines" = "PHL",
    "Poland" = "POL",
    "Portugal" = "PRT",
    "Qatar" = "QAT",
    "Romania" = "ROU",
    "Russia" = "RUS",
    "Rwanda" = "RWA",
    "Saint Lucia" = "LCA",
    "Saint Vincent and the Grenadines" = "VCT",
    "Samoa" = "WSM",
    "Sao Tome and Principe" = "STP",
    "Saudi Arabia" = "SAU",
    "Senegal" = "SEN",
    "Serbia" = "SRB",
    "Seychelles" = "SYC",
    "Sierra Leone" = "SLE",
    "Singapore" = "SGP",
    "Slovakia" = "SVK",
    "Slovenia" = "SVN",
    "Solomon Islands" = "SLB",
    "Somalia" = "SOM",
    "South Africa" = "ZAF",
    "South Korea" = "KOR",
    "Spain" = "ESP",
    "Sri Lanka" = "LKA",
    "Sudan" = "SDN",
    "Suriname" = "SUR",
    "Swaziland" = "SWZ",
    "Sweden" = "SWE",
    "Switzerland" = "CHE",
    "Syria" = "SYR",
    "Taiwan" = "TWN",
    "Tajikistan" = "TJK",
    "Tanzania" = "TZA",
    "Thailand" = "THA",
    "The Bahamas" = "BHS",
    "The Gambia" = "GMB",
    "Timor-Leste" = "TLS",
    "Togo" = "TGO",
    "Tonga" = "TON",
    "Trinidad and Tobago" = "TTO",
    "Tunisia" = "TUN",
    "Turkey" = "TUR",
    "Turkmenistan" = "TKM",
    "Uganda" = "UGA",
    "Ukraine" = "UKR",
    "United Arab Emirates" = "ARE",
    "United Kingdom" = "GBR",
    "United States" = "USA",
    "Uruguay" = "URY",
    "Uzbekistan" = "UZB",
    "Vanuatu" = "VUT",
    "Venezuela" = "VEN",
    "Vietnam" = "VNM",
    "Yemen" = "YEM",
    "Zambia" = "ZMB",
    "Zimbabwe" = "ZWE")
}
