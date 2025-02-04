#' Data for CYP2C9 diplotypes
#' 
#' Contains CYP2C9 diplotypes and corresponding activity scores and metabolizer phenotypes
#' 
#' @format A data frame with 3 variables and 15 rows
#' \describe{
#'     \item{CYP2C9.Diplotype}{CYP2C9 star allele}
#'     \item{Activity.Score}{CYP2C9 activity score based on the star allele}
#'     \item{Phenotype}{CYP2C9 metabolizer phenotype based on the star allele}
#'     }
#'     
#'@source {Created in house based on the alleles tested}
#'
#'@examples
#'data(CYP2C9_data)       #Lazy loading
"CYP2C9_data"

#' Data for CYP2C19 diplotypes
#' 
#' Contains CYP2C19 diplotypes and corresponding activity scores and metabolizer phenotypes
#' 
#' @format A data frame with 3 variables and 23 rows
#' \describe{
#'     \item{CYP2C19.Diplotype}{CYP2C19 star allele}
#'     \item{Activity.Score}{CYP2C19 activity score based on the star allele}
#'     \item{Phenotype}{CYP2C19 metabolizer phenotype based on the star allele}
#'     }
#'     
#'@source {Created in house based on the alleles tested}
#'
#'@examples
#'data(CYP2C19_data)       #Lazy loading
"CYP2C19_data"

#' Data for CYP2D6 star alleles
#' 
#' Contains CYP2D6 star alleles and corresponding activity scores 
#' 
#' @format A data frame with 3 variables and 125 rows
#' \describe{
#'     \item{CYP2D6.Diplotype}{CYP2D6 star allele}
#'     \item{Activity.Score}{CYP2D6 activity score based on the star allele}
#'     \item{Phenotype}{CYP2D6 metabolizer phenotype based on the star allele}
#'     }
#'     
#'@source {Created in house based on the CYP2D6 tier 1 variant alleles. Also includes: 8, 14, 2A, 36+10, 4.013}
#'
#'@examples
#'data(CYP2D6_data)       #Lazy loading
"CYP2D6_data"

#' Data for CYP2B6 diplotypes
#' 
#' Contains CYP2B6 diplotypes and corresponding activity scores and metabolizer phenotypes
#' 
#' @format A data frame with 3 variables and 16 rows
#' \describe{
#'     \item{CYP2B6.Diplotype}{CYP2B6 star allele}
#'     \item{Activity.Score}{CYP2B6 activity score based on the star allele}
#'     \item{Phenotype}{CYP2B6 metabolizer phenotype based on the star allele}
#'     }
#'     
#'@source {Created in house based on the CYP2B6 alleles tested}
#'
#'@examples
#'data(CYP2B6_data)       #Lazy loading
"CYP2B6_data"

#' Data for CYP3A5 diplotypes
#' 
#' Contains CYP3A5 diplotypes and corresponding activity scores and metabolizer phenotypes
#' 
#' @format A data frame with 3 variables and 21 rows
#' \describe{
#'     \item{CYP3A5.Diplotype}{CYP3A5 star allele}
#'     \item{Activity.Score}{CYP3A5 activity score based on the star allele}
#'     \item{Phenotype}{CYP3A5 metabolizer phenotype based on the star allele}
#'     }
#'     
#'@source {Downloaded from PharmGKB}
#'
#'@examples
#'data(CYP3A5_data)       #Lazy loading
"CYP3A5_data"

#' Sample Data with CYP2D6 genotypes
#' 
#' Contains sample genotypes for 15 CYP2D6 snps
#' 
#' @format A data frame with 15 variables and 4 rows
#' \describe{
#'     \item{CYP2D6_CNVInt6}{CYP2D6 snp}
#'     \item{CYP2D6_CNVx9}{CYP2D6 snp}
#'     \item{CYP2D6_rs1065852}{CYP2D6 snp}
#'     \item{CYP2D6_rs1080985}{CYP2D6 snp}
#'     \item{CYP2D6_rs1135840}{CYP2D6 snp}
#'     \item{CYP2D6_rs16947}{CYP2D6 snp}
#'     \item{CYP2D6_rs28371706}{CYP2D6 snp}
#'     \item{CYP2D6_rs28371725}{CYP2D6 snp}
#'     \item{CYP2D6_rs35742686}{CYP2D6 snp}
#'     \item{CYP2D6_rs3892097}{CYP2D6 snp}
#'     \item{CYP2D6_rs5030655}{CYP2D6 snp}
#'     \item{CYP2D6_rs5030656}{CYP2D6 snp}
#'     \item{CYP2D6_rs5030865A}{CYP2D6 snp}
#'     \item{CYP2D6_rs5030865T}{CYP2D6 snp}
#'     \item{CYP2D6_rs59421388}{CYP2D6 snp}
#'     }
#'     
#'@source {Created in house to be used as an example}
#'
#'@examples
#'data(sample_CYP2D6_genotypes)       #Lazy loading
"sample_CYP2D6_genotypes"

#' Allele Definition Table for CYP2B6 
#' 
#' Contains CYP2B6 information about which variants define star alleles. Preprocessed to only include star allele and rs ID columns.
#' 
#' @format A data frame with 49 variables and 55 rows
#' \describe{
#'     \item{CYP2B6.Allele}{CYP2B6 star alleles}
#'     \item{rsID}{CYP2B6  rs IDs}
#'     }
#'     
#'@source {Downloaded from Pharmgkb}
#'
#'@examples
#'data(CYP2B6_Allele_def)       #Lazy loading
"CYP2B6_Allele_def"

#' Allele Definition Table for CYP2C9 
#' 
#' Contains CYP2C9 information about which variants define star alleles. Preprocessed to only include star allele and rs ID columns.
#' 
#' @format A data frame with 81 variables and 86 rows
#' \describe{
#'     \item{CYP2C9.Allele}{CYP2C9 star alleles}
#'     \item{rsID}{CYP2C9  rs IDs}
#'     }
#'     
#'@source {Downloaded from Pharmgkb}
#'
#'@examples
#'data(CYP2C9_Allele_def)       #Lazy loading
"CYP2C9_Allele_def"

#' Allele Definition Table for CYP2C19 
#' 
#' Contains CYP2C19 information about which variants define star alleles. Preprocessed to only include star allele and rs ID columns.
#' 
#' @format A data frame with 81 variables and 86 rows
#' \describe{
#'     \item{CYP2C19.Allele}{CYP2C9 star alleles}
#'     \item{rsID}{CYP2C19  rs IDs}
#'     }
#'     
#'@source {Downloaded from Pharmgkb}
#'
#'@examples
#'data(CYP2C19_Allele_def)       #Lazy loading
"CYP2C19_Allele_def"

#' Allele Definition Table for CYP2D6 
#' 
#' Contains CYP2D6 information about which variants define star alleles. Preprocessed to only include star allele and rs ID columns.
#' 
#' @format A data frame with 81 variables and 86 rows
#' \describe{
#'     \item{CYP2D6.Allele}{CYP2C9 star alleles}
#'     \item{rsID}{CYP2D6  rs IDs}
#'     }
#'     
#'@source {Downloaded from Pharmgkb}
#'
#'@examples
#'data(CYP2D6_Allele_def)       #Lazy loading
"CYP2D6_Allele_def"

#' Allele Definition Table for CYP1A2 
#' 
#' Contains CYP1A2 information about which variants define star alleles. Preprocessed to only include star allele and rs ID columns.
#' 
#' @format A data frame with 42 variables and 47 rows
#' \describe{
#'     \item{CYP1A2}{CYP1A2 star alleles}
#'     \item{rsID}{CYP1A2 rs IDs}
#'     }
#'     
#'@source {Downloaded from Pharmgkb}
#'
#'@examples
#'data(CYP1A2_Allele_def)       #Lazy loading
"CYP1A2_Allele_def"

#' Allele Definition Table for CYP3A4 
#' 
#' Contains CYP3A4 information about which variants define star alleles. Preprocessed to only include star allele and rs ID columns.
#' 
#' @format A data frame with 46 variables and 43 rows
#' \describe{
#'     \item{CYP3A4}{CYP3A4 star alleles}
#'     \item{rsID}{CYP3A4 rs IDs}
#'     }
#'     
#'@source {Downloaded from Pharmgkb}
#'
#'@examples
#'data(CYP3A4_Allele_def)       #Lazy loading
"CYP3A4_Allele_def"

#' Allele Definition Table for CYP3A5 
#' 
#' Contains CYP3A5 information about which variants define star alleles. Preprocessed to only include star allele and rs ID columns.
#' 
#' @format A data frame with 6 variables and 7 rows
#' \describe{
#'     \item{CYP3A5}{CYP3A5 star alleles}
#'     \item{rsID}{CYP3A5 rs IDs}
#'     }
#'     
#'@source {Downloaded from Pharmgkb}
#'
#'@examples
#'data(CYP3A5_Allele_def)       #Lazy loading
"CYP3A5_Allele_def"

#' Allele Definition Table V2 for CYP1A2 
#' 
#' Contains CYP1A2 information about which variants define star alleles based on the update from Pharmvar in Dec. 2024. Preprocessed to only include star allele and rs ID columns.
#' 
#' @format A data frame with 42 variables and 47 rows
#' \describe{
#'     \item{CYP1A2}{CYP1A2 star alleles}
#'     \item{rsID}{CYP1A2 rs IDs}
#'     }
#'     
#'@source {Downloaded from Pharmvar}
#'
#'@examples
#'data(CYP1A2_Allele_def_2)       #Lazy loading
"CYP1A2_Allele_def_2"