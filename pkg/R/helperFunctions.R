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
  out
}

#' Create CSMF plot with Undetermined for InterVA4 Output
#'
#' @param x A fitted object returned from openVA::codeVA()
#' @param top An integer specifying the number of top causes included
#' in the results.
#'
#' @return 
#' @export
#'
CSMF2 <- function (x, top) {

  csmf  <- InterVA4::CSMF(x, InterVA.rule=TRUE, noplot=TRUE, min.prob=.00001)
  csmf2 <- sort(csmf, decreasing=TRUE)[1:top]
  isUndet <- "Undetermined" %in% names(csmf2)
  if(isUndet){
    idUndet <- which(names(csmf2) == "Undetermined")
    newBars <- c(csmf2[-idUndet], csmf2[idUndet])
    par(las = 2)
    par(mar = c(5, 15, 4, 2))
    barplot(newBars[top:1], horiz=TRUE, names.arg=names(newBars)[top:1],
            cex.names = 0.8, xlab = "Probability",
            col = rev(grey.colors(length(newBars))))
  } else{
    InterVA4::CSMF(x, InterVA.rule=TRUE, min.prob=.00001)
  }
}

#' Return list of countries and abbreviations for Tariff2 input
#'
#'
#' @return A list of country names and abbreviations.
#'
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
    "Cote d’Ivoire" = "CIV",
    "Croatia" = "HRV",
    "Cuba" = "CUB",
    "Cyprus" = "CY
    P",
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
