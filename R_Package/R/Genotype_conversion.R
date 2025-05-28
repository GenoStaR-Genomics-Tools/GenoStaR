
#' @title Find Missing Genotypes 
#' @description This function checks the input genotype dataframe for any missing genotype entries. If genotypes are missing, they will be returned in a list specifying which patient index and column contains the missing data.
#' @param df A data frame containing genotype data with specific columns. Each row represents a patient, and columns represent genotypes for specific SNPs. 
#' Column names should follow the format "gene_snp", for example, "CYP2D6_rs1065852".
#' @return A list of character vectors, where each element specifies the patient index (row) and column name (SNP) of missing data (NA).
#' @examples
#' find_missing_data(data)
#'
#' @export
find_missing_data <- function(df) {
  #check which positions contain empty values
  missing_positions <- which(df == " ", arr.ind = TRUE)
  if (nrow(missing_positions) == 0) {
    return("No Missing Genotypes")
  }
  missing_messages <- apply(missing_positions, 1, function(pos) {
    #get positions of missing data
    row <- pos[1]
    col <- pos[2]
    #label the positions with the specific patient index and snp column
    patient_index <- df$Patient_index[row]
    col_name <- colnames(df)[col]
    paste("Missing Data: patient index", patient_index, ",", col_name)
  })
  return(as.list(missing_messages))
}

#' @title Fill empty cells in dataframe
#' @description This function iterates through a dataframe (allele definition table in this use case) and fills the empty cells with the genotypes from the first row (wildtype)
#' @param df A data frame containing genotype data with specific columns. The column names should be snps, for example rs1065852.
#' @return A data frame with all the empty cells assigned the wildtype for that snp. 
#' @examples
#' data_filled <- fill_empty_cells(data_filtered)
#' @export
fill_empty_cells <- function(df) {

  for (i in seq_along(df)) {
    df[, i][df[, i] == ""] <- df[1, i]
  }
  return(df)
}


#' @title Call CYP2D6 Diplotype
#' @description This function will use index polymorphisms to determine the CYP2D6 diplotype. If no diplotype is found, it will expand the attempt using more of the allele definition table. Used within the find_diplotype() function, not meant to be used individually. 
#' @param genotypes the preprocessed genotypes, ready to be assessed as haplotypes
#' @param row_idx indicates current row of the dataframe 
#' @param diplotypes a character vector to store any found diplotype matches
#' @param cnv_x9 copy number information for exon 9
#' @param cnv_int6 copy number information for intron 6
#' @param cnv_5utr copy information for 5 prime utr 
#' @param score the current score designated for each found star allele 
#' @param total_score the combined scores to give the total score associated with a diplotype
#' @return A matrix containing any unmatched haplotypes, indicating further processing is required
#' @export
process_cyp2d6 <- function(genotypes, row_idx, diplotypes, cnv_x9, cnv_int6, cnv_5utr, score, total_score) {
  #initialize variables 
  allele_matches <- vector("character", 2)
  genotype_match <- FALSE
  seq_unmatched <- matrix(nrow = 0, ncol = 0)
 #if there's no CNV information, explicitly set to NA for handling
  if(length(cnv_x9) <= 0){
    cnv_x9 <- NA
  }
  if(length(cnv_int6) <= 0){
    cnv_int6 <- NA
  }
  if(length(cnv_5utr) <= 0){
    cnv_5utr <- NA
  }

  #Extract and format CNV information 
  #Split the "1 1" and "2 2" strings into individual numeric values (for checking *36)
  cnv_x9_vals <- if (length(cnv_x9) > 0 && (!is.na(cnv_x9))) as.numeric(unlist(strsplit(as.character(cnv_x9), " "))) else NA
  cnv_int6_vals <- if (length(cnv_int6) > 0 && (!is.na(cnv_int6))) as.numeric(unlist(strsplit(as.character(cnv_int6), " "))) else NA
  CNV_flag <- TRUE
  #Find the initial star alleles based on the genotypes, focused on the tier 1 variants 
  #iterate through one haplotype at a time 
  for(i in 1:2){
    points <- 0
    if(i == 2){
      CNV_flag <- FALSE
    }
    star_allele <- "" #initialize as empty 

    if(CNV_flag){ #check cnv conditions only once
      if(length(cnv_x9) > 0 & length(cnv_int6) > 0 & length(cnv_5utr) > 0 & !is.na(cnv_x9) & !is.na(cnv_int6) & !is.na(cnv_5utr)) { # All CNV info is provided
        #check star alleles that are based on cnv information and genotypes at specific snps 
        if (sum(cnv_x9_vals) < sum(cnv_int6_vals) && "CYP2D6_rs1065852" %in% rownames(genotypes) && !"A" %in% genotypes["CYP2D6_rs1065852", i] && 
            "CYP2D6_rs1135840" %in% rownames(genotypes) && !"G" %in% genotypes["CYP2D6_rs1135840", i]){
          star_allele <- "*36"
          genotype_match <- TRUE
          points <- 100
        }else if (all(cnv_x9 == "1 1") & all(cnv_int6 == "1 1") & all(cnv_5utr == "1 1")){ 
          star_allele <- "*5"  
          genotype_match <- TRUE
          points <- 100
        }else if (all(cnv_x9 == "2 2") & (all(cnv_int6 == "1 1") | all(cnv_int6 == "2 2")) & all(cnv_5utr == "1 1")){
          star_allele <- "*13"  
          genotype_match <- TRUE
          points <- 100
        }else if ((all(cnv_x9 == "1 1") | all(cnv_x9 == "0 0")) & (all(cnv_int6 == "1 1") | all(cnv_int6 == "0 0")) & all(cnv_5utr == "2 2")){
          star_allele <- "*68" 
          genotype_match <- TRUE
          points <- 100
        }
        else if(all(cnv_x9 == "2 2") & all(cnv_int6 == "2 2") & all(cnv_5utr == "3 3") & "CYP2D6_rs3892097" %in% rownames(genotypes) & "T" %in% genotypes["CYP2D6_rs3892097", i]){
          star_allele <- "*68+*4"
          genotype_match <- TRUE
          points <- 100
        }
        else if(sum(cnv_x9_vals) < sum(cnv_int6_vals) && "CYP2D6_rs1065852" %in% rownames(genotypes) && "A" %in% genotypes["CYP2D6_rs1065852", i] && 
                "CYP2D6_rs1135840" %in% rownames(genotypes) && "G" %in% genotypes["CYP2D6_rs1135840", i] && 
                "CYP2D6_rs16947" %in% rownames(genotypes) && !"A" %in% genotypes["CYP2D6_rs16947", i]){
          star_allele <- "*36+*10"
          genotype_match <- TRUE
          points <- 100
        }
      }
      else if(length(cnv_x9) > 0 & length(cnv_int6) > 0 & (length(cnv_5utr) <= 0 | is.na(cnv_5utr)) & !is.na(cnv_x9) & !is.na(cnv_int6)){ #if only given cnv info for exon 9 and intron 6 provided 
        if (sum(cnv_x9_vals) < sum(cnv_int6_vals) && "CYP2D6_rs1065852" %in% rownames(genotypes) && !"A" %in% genotypes["CYP2D6_rs1065852", i] && 
            "CYP2D6_rs1135840" %in% rownames(genotypes) && !"G" %in% genotypes["CYP2D6_rs1135840", i]){
          star_allele <- "*36"
          genotype_match <- TRUE
          points <- 100
        }else if (all(cnv_x9 == "1 1") & all(cnv_int6 == "1 1")){ 
          star_allele <- "*5"  
          genotype_match <- TRUE
          points <- 100
        }
        else if (all(cnv_x9 == "2 2") & (all(cnv_int6 == "1 1") | all(cnv_int6 == "2 2"))){
          star_allele <- "*13"  
          genotype_match <- TRUE
          points <- 100
        }else if ((all(cnv_x9 == "1 1") | all(cnv_x9 == "0 0")) & (all(cnv_int6 == "1 1") | all(cnv_int6 == "0 0"))){
          star_allele <- "*68" 
          genotype_match <- TRUE
          points <- 100
        }
        else if(all(cnv_x9 == "2 2") & all(cnv_int6 == "2 2") & "CYP2D6_rs3892097" %in% rownames(genotypes) & "T" %in% genotypes["CYP2D6_rs3892097", i]){
          star_allele <- "*68+*4"
          genotype_match <- TRUE
          points <- 100
        }
        else if(sum(cnv_x9_vals) < sum(cnv_int6_vals) && "CYP2D6_rs1065852" %in% rownames(genotypes) && "A" %in% genotypes["CYP2D6_rs1065852", i] && 
                "CYP2D6_rs1135840" %in% rownames(genotypes) && "G" %in% genotypes["CYP2D6_rs1135840", i] && 
                "CYP2D6_rs16947" %in% rownames(genotypes) && !"A" %in% genotypes["CYP2D6_rs16947", i]){
          star_allele <- "*36+*10"
          genotype_match <- TRUE
          points <- 100
        }
      }
      else if(length(cnv_x9 > 0) & !is.na(cnv_x9) & (length(cnv_int6) <= 0 | is.na(cnv_int6)) & (length(cnv_5utr) <= 0 | is.na(cnv_5utr))){ #when only given x9 copy information
        if (all(cnv_x9 == "1 1")){ 
          star_allele <- "*5"  
          genotype_match <- TRUE
          points <- 100
        }
      }
      else{#if no cnv information given, skip the above conditionals and return warning 
        warning("No CNV columns found in the dataframe. Proceeding without CNVs.")
        CNV_flag <- FALSE
      }
      if (!genotype_match) {
        # Manually set CNV_flag to FALSE and allow the next block to run
        CNV_flag <- FALSE
      }
    }
    if (!CNV_flag || (!genotype_match && i == 1) || i == 2) {
      # Check for *4 star allele (e.g., rs3892097)
      if ("CYP2D6_rs3892097" %in% rownames(genotypes) && "T" %in% genotypes["CYP2D6_rs3892097", i]){ 
        star_allele  <- "*4"
        genotype_match <- TRUE
        points <- 50
      }
      # Check for *3 star allele (e.g., rs35742686)
      else if ("CYP2D6_rs35742686" %in% rownames(genotypes) && "-" %in% genotypes["CYP2D6_rs35742686", i]) {
        star_allele  <- "*3"
        genotype_match <- TRUE
        points <- 50
      }
      # Check for *6 star allele (e.g., rs5030655)
      else if ("CYP2D6_rs5030655" %in% rownames(genotypes) && "-" %in% genotypes["CYP2D6_rs5030655", i]) {
        star_allele  <- "*6"
        genotype_match <- TRUE
        points <- 50
      }
      # Check for *9 star allele (e.g., rs5030656)
      else if ("CYP2D6_rs5030656" %in% rownames(genotypes) && "-" %in% genotypes["CYP2D6_rs5030656", i]) {
        star_allele  <- "*9"
        genotype_match <- TRUE
        points <- 50
      }
      # Check for *8 star allele (e.g., rs5030865)
      else if ("CYP2D6_rs5030865A" %in% rownames(genotypes) && "A" %in% genotypes["CYP2D6_rs5030865A", i] | 
               "CYP2D6_rs5030865T" %in% rownames(genotypes) && "A" %in% genotypes["CYP2D6_rs5030865T", i]){
        star_allele  <- "*14"
        genotype_match <- TRUE
        points <- 60
      }
      else if("CYP2D6_rs5030865A" %in% rownames(genotypes) && "T" %in% genotypes["CYP2D6_rs5030865A", i] |
              "CYP2D6_rs5030865T" %in% rownames(genotypes) && "T" %in% genotypes["CYP2D6_rs5030865T", i]){
        star_allele  <- "*8"
        genotype_match <- TRUE
        points <- 60
      }
      else if("CYP2D6_rs5030867" %in% rownames(genotypes) && "G" %in% genotypes["CYP2D6_rs5030867", i]){
        star_allele <- "*7"
        genotype_match <- TRUE
        points <- 50
      }
      else if("CYP2D6_rs201377835" %in% rownames(genotypes) && "G" %in% genotypes["CYP2D6_rs201377835", i]){
        star_allele <- "*11"
        genotype_match <- TRUE
        points <- 50
      }
      else if("CYP2D6_rs28371706" %in% rownames(genotypes) && "A" %in% genotypes["CYP2D6_rs28371706", i] && 
              "CYP2D6_rs1135840" %in% rownames(genotypes) && "G" %in% genotypes["CYP2D6_rs1135840", i] && 
              "CYP2D6_rs16947" %in% rownames(genotypes) && "A" %in% genotypes["CYP2D6_rs16947", i]){
        star_allele <- "*17"
        genotype_match <- TRUE
        points <- 70
      }
      else if("CYP2D6_rs59421388" %in% rownames(genotypes) && "T" %in% genotypes["CYP2D6_rs59421388", i]){ 
        star_allele  <- "*29"
        genotype_match <- TRUE
        points <- 50
      }
      else if("CYP2D6_rs28371725" %in% rownames(genotypes) && "T" %in% genotypes["CYP2D6_rs28371725", i] && 
              "CYP2D6_rs16947" %in% rownames(genotypes) && "A" %in% genotypes["CYP2D6_rs16947", i] &&
              "CYP2D6_rs1135840" %in% rownames(genotypes) && "G" %in% genotypes["CYP2D6_rs1135840", i]){
        star_allele  <- "*41"
        genotype_match <- TRUE
        points <- 70
      }
      
      #Check for *2 star allele
      else if ("CYP2D6_rs16947" %in% rownames(genotypes) && "A" %in% genotypes["CYP2D6_rs16947", i] &&
               "CYP2D6_rs1135840" %in% rownames(genotypes) && "G" %in% genotypes["CYP2D6_rs1135840", i]) {
        star_allele <- "*2"
        genotype_match <- TRUE
        points <- 40 
      }
      #Check for *10 
      else if("CYP2D6_rs1065852" %in% rownames(genotypes) && "A" %in% genotypes["CYP2D6_rs1065852", i] && 
              "CYP2D6_rs1135840" %in% rownames(genotypes) && "G" %in% genotypes["CYP2D6_rs1135840", i] && 
              "CYP2D6_rs16947" %in% rownames(genotypes) && !"A" %in% genotypes["CYP2D6_rs16947", i] && 
              "CYP2D6_rs3892097" %in% rownames(genotypes) && !"T" %in% genotypes["CYP2D6_rs3892097", i]){
        star_allele  <- "*10"
        genotype_match <- TRUE
        points <- 40 
      }
      else if("CYP2D6_rs1065852" %in% rownames(genotypes) && "G" %in% genotypes["CYP2D6_rs1065852", i] &&
              "CYP2D6_rs1135840" %in% rownames(genotypes) && "C" %in% genotypes["CYP2D6_rs1135840", i] &&
              "CYP2D6_rs16947" %in% rownames(genotypes) && "G" %in% genotypes["CYP2D6_rs16947", i] &&
              "CYP2D6_rs28371706" %in% rownames(genotypes) && "G" %in% genotypes["CYP2D6_rs28371706", i] &&
              "CYP2D6_rs28371725" %in% rownames(genotypes) && "C" %in% genotypes["CYP2D6_rs28371725", i] &&
              "CYP2D6_rs35742686" %in% rownames(genotypes) && "T" %in% genotypes["CYP2D6_rs35742686", i] &&
              "CYP2D6_rs3892097" %in% rownames(genotypes) && "C" %in% genotypes["CYP2D6_rs3892097", i] &&
              "CYP2D6_rs5030655" %in% rownames(genotypes) && "A" %in% genotypes["CYP2D6_rs5030655", i] &&
              "CYP2D6_rs5030656" %in% rownames(genotypes) && "TCT" %in% genotypes["CYP2D6_rs5030656", i] &&
              "CYP2D6_rs5030865A" %in% rownames(genotypes) && "C" %in% genotypes["CYP2D6_rs5030865A", i] &&
              "CYP2D6_rs5030865T" %in% rownames(genotypes) && "C" %in% genotypes["CYP2D6_rs5030865T", i] &&
              "CYP2D6_rs59421388" %in% rownames(genotypes) && "C" %in% genotypes["CYP2D6_rs59421388",i ]){
        star_allele <- "*1"
        genotype_match <- TRUE 
        points <- 40
      }
      else if("CYP2D6_rs1065852" %in% rownames(genotypes) && "G" %in% genotypes["CYP2D6_rs1065852", i] &&
              "CYP2D6_rs1135840" %in% rownames(genotypes) && "C" %in% genotypes["CYP2D6_rs1135840", i] &&
              "CYP2D6_rs16947" %in% rownames(genotypes) && "G" %in% genotypes["CYP2D6_rs16947", i] &&
              "CYP2D6_rs28371706" %in% rownames(genotypes) && "G" %in% genotypes["CYP2D6_rs28371706", i] &&
              "CYP2D6_rs28371725" %in% rownames(genotypes) && "C" %in% genotypes["CYP2D6_rs28371725", i] &&
              "CYP2D6_rs35742686" %in% rownames(genotypes) && "T" %in% genotypes["CYP2D6_rs35742686", i] &&
              "CYP2D6_rs3892097" %in% rownames(genotypes) && "C" %in% genotypes["CYP2D6_rs3892097", i] &&
              "CYP2D6_rs5030655" %in% rownames(genotypes) && "A" %in% genotypes["CYP2D6_rs5030655", i] &&
              "CYP2D6_rs5030656" %in% rownames(genotypes) && "TCT" %in% genotypes["CYP2D6_rs5030656", i] &&
              "CYP2D6_rs5030865A" %in% rownames(genotypes) && "C" %in% genotypes["CYP2D6_rs5030865A", i] &&
              "CYP2D6_rs5030865T" %in% rownames(genotypes) && "C" %in% genotypes["CYP2D6_rs5030865T", i]){
        star_allele <- "*1"
        genotype_match <- TRUE
        points <- 40
      }
    } 
    if(genotype_match){
      allele_matches[i] <- star_allele
      score[i] <- points
    }
  }
    #if allele matches variable isn't empty 
    if (all(allele_matches != "")) {
      # Format diplotype for both alleles
      diplotype <- paste(allele_matches[1], allele_matches[2], sep = "/")
      total_score <- score[1] + score[2]
      #check for special cases based on cnv info available
      if(length(cnv_x9) > 0 & !is.na(cnv_x9) & length(cnv_int6) > 0 & !is.na(cnv_int6) & length(cnv_5utr) > 0 & !is.na(cnv_5utr)){ #if all cnv info available
        if(all(cnv_x9 == "2 2") & all(cnv_int6 == "4 4") & all(cnv_5utr == "4 4")){
          diplotype <- "*36+*10/*36+*10"
          total_score <- 100
        }
        if(all(cnv_x9 == "2 2") & all(cnv_int6 == "3 3") & all(cnv_5utr == "3 3") & (allele_matches[1] == "*10" | allele_matches[2] == "*10")){
          diplotype <- "*10/*36+*10"
          total_score <- 100
        }
        if(allele_matches[1] == "*4" & allele_matches[2] == "*4" & all(cnv_x9 == "2 2") & all(cnv_int6 == "2 2") & all(cnv_5utr == "3 3")){
          diplotype <- "*4/*68+*4"
          total_score <- 100
        }
      }
      else if(length(cnv_x9) > 0 & length(cnv_int6) > 0 & length(cnv_5utr) <= 0){
        if(all(cnv_x9 == "2 2") & all(cnv_int6 == "4 4")){
          diplotype <- "*36+*10/*36+*10"
          total_score <- 100
        }
        if(all(cnv_x9 == "2 2") & all(cnv_int6 == "3 3") & (allele_matches[1] == "*10" | allele_matches[2] == "*10")){
          diplotype <- "*10/*36+*10"
          total_score <- 100
        }
        if(allele_matches[1] == "*4" & allele_matches[2] == "*4" & all(cnv_x9 == "2 2") & all(cnv_int6 == "2 2")){
          diplotype <- "*4/*68+*4"
          total_score <- 100
        }
      }
     #save the diplotype, calculate total score and then reset variables 
      diplotypes[row_idx] <- diplotype
      no_match <- 0
      one_match <- NA
      matched_allele <- NA
      total_score <- score[1] + score[2]
    }
    
    #Handle unmatched sequences if any star alleles couldn't be found from index polymorphisms
    else if((allele_matches[1] != "" & allele_matches[2] == "") | (allele_matches[1] == "" & allele_matches[2] != "")){
      no_match <- 1
      one_match <- which(allele_matches != "")
      matched_allele <- allele_matches[one_match]
      score[one_match] <- score[one_match]
      total_score <- score[one_match]
    }
    #the case where no matches were found 
    else if(all(allele_matches == "")){
      no_match <- 2
      one_match <- NA
      matched_allele <- NA
      total_score <- 0
    }
  return(list(no_match, diplotypes, matched_allele, one_match, score, total_score))
}

#' @title Call CYP1A2 Diplotype
#' @description This function will use index polymorphisms to determine the CYP1A2 diplotype. If no diplotype is found, it will expand the attempt using more of the allele definition table. Used within the find_diplotype() function, not meant to be used individually. 
#' @param genotypes the preprocessed genotypes, ready to be assessed as haplotypes
#' @param row_idx indicates current row of the dataframe 
#' @param diplotypes a character vector to store any found diplotype matches
#' @param score the current score designated for each found star allele 
#' @param total_score the combined scores to give the total score associated with a diplotype
#' @param CYP1A2_name indicates which nomenclature version to use, either new or old
#' @return A matrix containing any unmatched haplotypes, indicating further processing is required
#' @export
process_cyp1a2 <- function(genotypes, row_idx, diplotypes, score, total_score, CYP1A2_name){
  #initialize variables 
  allele_matches <- vector("character", 2)
  genotype_match <- FALSE
  seq_unmatched <- matrix(nrow = 0, ncol = 0)
  
  #Find the initial star alleles based on the genotypes
  #iterate through one haplotype at a time 
  #conditions to assign star alleles based on index polymorphisms
  if(CYP1A2_name == "new"){
    for(i in 1:2){
      points <- 0
      if ("CYP1A2_rs56160784" %in% rownames(genotypes) && "G" %in% genotypes["CYP1A2_rs56160784", i]) {
        star_allele <- "*2"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs56276455" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs56276455", i]) {
        star_allele <- "*3"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs72547516" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs72547516", i] &&
                 "CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i]) {
        star_allele <- "*4"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs55889066" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs55889066", i]) {
        star_allele <- "*5"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs28399424" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs28399424", i]) {
        star_allele <- "*6"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs56107638" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs56107638", i]) {
        star_allele <- "*7"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs72547517" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs72547517", i]) {
        star_allele <- "*8"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                 "CYP1A2_rs138652540" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs138652540", i]) {
        star_allele <- "*9"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                 "CYP1A2_rs72547512" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_rs72547512", i]) {
        star_allele <- "*10"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs72547513" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs72547513", i]) {
        star_allele <- "*11"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs758748797" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs758748797", i]) {
        star_allele <- "*12"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                 "CYP1A2_rs35796837" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs35796837", i]) {
        star_allele <- "*13"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                 "CYP1A2_rs45486893" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs45486893", i]) {
        star_allele <- "*14"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                 "CYP1A2_rs72547511" %in% rownames(genotypes) && "G" %in% genotypes["CYP1A2_rs72547511", i]) {
        star_allele <- "*15"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs72547515" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs72547515", i]) {
        star_allele <- "*16"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
               "CYP1A2_rs149928755" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs149928755", i]) {
        star_allele <- "*17"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_X3468AC" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_X3468AC", i]) {
        star_allele <- "*18"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs374094758" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs374094758", i]) {
        star_allele <- "*19"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs144148965" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs144148965", i]) {
        star_allele <- "*20"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                "CYP1A2_rs17861157" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs17861157", i] &&
                "CYP1A2_rs143193369" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs143193369", i]){
        star_allele <- "*21"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                "CYP1A2_rs17861157" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs17861157", i] &&
                "CYP1A2_rs143193369" %in% rownames(genotypes) && !"A" %in% genotypes["CYP1A2_rs143193369", i]){
        star_allele <- "*22"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                  "CYP1A2_rs17861152" %in% rownames(genotypes) && "G" %in% genotypes["CYP1A2_rs17861152", i] &&
                  "CYP1A2_rs17861157" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs17861157", i] &&
                  "CYP1A2_rs143193369" %in% rownames(genotypes) && !"A" %in% genotypes["CYP1A2_rs143193369", i]){
        star_allele <- "*23"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                "CYP1A2_rs568839929" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs568839929", i]){
        star_allele <- "*24"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                  "CYP1A2_rs201537008" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs201537008", i]){
        star_allele <- "*25"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                "CYP1A2_rs45442197" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs45442197", i]){
        star_allele <- "*26"
        genotype_match <- TRUE
        points <- 50
      } else if ("CYP1A2_rs758752876" %in% rownames(genotypes) && "-" %in% genotypes["CYP1A2_rs758752876", i]){
        star_allele <- "*27"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs60086777" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs60086777", i]){
        star_allele <- "*28"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                "CYP1A2_rs571663822" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs571663822", i]){
        star_allele <- "*29"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                "CYP1A2_rs17861157" %in% rownames(genotypes) && !"A" %in% genotypes["CYP1A2_rs17861157", i] &&
                "CYP1A2_rs143193369" %in% rownames(genotypes) && !"A" %in% genotypes["CYP1A2_rs143193369", i] &&
                "CYP1A2_rs17861152" %in% rownames(genotypes) && !"G" %in% genotypes["CYP1A2_rs17861152", i] &&
                "CYP1A2_rs45442197" %in% rownames(genotypes) && !"T" %in% genotypes["CYP1A2_rs45442197", i]  &&
                "CYP1A2_rs201537008" %in% rownames(genotypes) && !"A" %in% genotypes["CYP1A2_rs201537008", i] &&
                "CYP1A2_rs568839929" %in% rownames(genotypes) && !"A" %in% genotypes["CYP1A2_rs568839929", i] &&
                "CYP1A2_rs149928755" %in% rownames(genotypes) && !"T" %in% genotypes["CYP1A2_rs149928755", i] &&
                "CYP1A2_rs72547511" %in% rownames(genotypes) && !"G" %in% genotypes["CYP1A2_rs72547511", i] &&
                "CYP1A2_rs45486893" %in% rownames(genotypes) && !"T" %in% genotypes["CYP1A2_rs45486893", i] &&
                "CYP1A2_rs35796837" %in% rownames(genotypes) && !"A" %in% genotypes["CYP1A2_rs35796837", i] &&
                "CYP1A2_rs72547512" %in% rownames(genotypes) && !"C" %in% genotypes["CYP1A2_rs72547512", i] &&
                "CYP1A2_rs72547512" %in% rownames(genotypes) && !"C" %in% genotypes["CYP1A2_rs72547512", i]&&
                "CYP1A2_rs138652540" %in% rownames(genotypes) && !"T" %in% genotypes["CYP1A2_rs138652540", i] &&
                "CYP1A2_rs72547516" %in% rownames(genotypes) && !"T" %in% genotypes["CYP1A2_rs72547516", i] &&
                "CYP1A2_rs147333000" %in% rownames(genotypes) && !"T" %in% genotypes["CYP1A2_rs147333000", i] &&
                "CYP1A2_rs192799115" %in% rownames(genotypes) && !"A" %in% genotypes["CYP1A2_rs192799115", i] &&
                "CYP1A2_rs148157092" %in% rownames(genotypes) && !"T" %in% genotypes["CYP1A2_rs148157092", i] &&
                "CYP1A2_rs45565238" %in% rownames(genotypes) && !"A" %in% genotypes["CYP1A2_rs45565238", i] &&
                "CYP1A2_rs34067076" %in% rownames(genotypes) && !"A" %in% genotypes["CYP1A2_rs34067076", i] &&
                "CYP1A2_rs34151816" %in% rownames(genotypes) && !"T" %in% genotypes["CYP1A2_rs34151816", i] &&
                "CYP1A2_rs199528490" %in% rownames(genotypes) && !"C" %in% genotypes["CYP1A2_rs199528490", i] &&
                "CYP1A2_rs200075745" %in% rownames(genotypes) && !"T" %in% genotypes["CYP1A2_rs200075745", i] &&
                "CYP1A2_rs45540640" %in% rownames(genotypes) && !"G" %in% genotypes["CYP1A2_rs45540640", i]){
        star_allele <- "*30"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                "CYP1A2_rs147333000" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs147333000", i]){
        star_allele <- "*31"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                "CYP1A2_rs192799115" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs192799115", i]){
        star_allele <- "*32"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs201763966" %in% rownames(genotypes) && "G" %in% genotypes["CYP1A2_rs201763966", i]){
        star_allele <- "*33"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                "CYP1A2_rs148157092" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs148157092", i]){
        star_allele <- "*34"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs45468096" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs45468096", i]){
        star_allele <- "*35"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                "CYP1A2_rs17861157" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs17861157", i] &&
                "CYP1A2_rs45565238" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs45565238", i] &&
                "CYP1A2_rs143193369" %in% rownames(genotypes) && !"A" %in% genotypes["CYP1A2_rs143193369", i]){
        star_allele <- "*36"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
               "CYP1A2_rs34067076" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs34067076", i]){
        star_allele <- "*37"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs150164960" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs150164960", i]){
        star_allele <- "*38"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs773105304" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs773105304", i]){
        star_allele <- "*39"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                "CYP1A2_rs34151816" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs34151816", i]){
        star_allele <- "*40"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                "CYP1A2_rs369511887" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs369511887", i]){
        star_allele <- "*41"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                "CYP1A2_rs202191520" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs202191520", i]){
        star_allele <- "*42"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                "CYP1A2_rs566851431" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs566851431", i]){
        star_allele <- "*43"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                "CYP1A2_rs199528490" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_rs199528490", i]){
        star_allele <- "*44"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                "CYP1A2_rs200075745" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs200075745", i]){
        star_allele <- "*45"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs201551575" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs201551575", i]){
        star_allele <- "*46"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
                "CYP1A2_rs45540640" %in% rownames(genotypes) && "G" %in% genotypes["CYP1A2_rs45540640", i]){
        star_allele <- "*47"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i]){
        star_allele <- "*30"
        genotype_match <- TRUE
        points <- 50
      }else {
        star_allele <- "*1"
        genotype_match <- TRUE
        points <- 50
      }
    if(genotype_match){
      allele_matches[i] <- star_allele
      score[i] <- points
      }
    }
  }
  else if(CYP1A2_name == "old"){
    for(i in 1:2){
      star_allele <- "" #initialize as empty
      if ("CYP1A2_rs2470890" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_rs2470890", i]){ 
        star_allele  <- "*1B"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs2069514" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_rs2069514", i]){ 
        star_allele  <- "*1C"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs35694136" %in% rownames(genotypes) && "-" %in% genotypes["CYP1A2_rs35694136", i]){ 
        star_allele  <- "*1D"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs2069526" %in% rownames(genotypes) && "G" %in% genotypes["CYP1A2_rs2069526", i]){ 
        star_allele  <- "*1E"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
               "CYP1A2_rs2069526" %in% rownames(genotypes) && !"G" %in% genotypes["CYP1A2_rs2069526", i]){ 
        star_allele  <- "*1F"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs2069526" %in% rownames(genotypes) && "G" %in% genotypes["CYP1A2_rs2069526", i] &&
               "CYP1A2_rs2470890" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_rs2470890", i]){ 
        star_allele  <- "*1G"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_X2025AC" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_X2025AC", i] &&
               "CYP1A2_rs2470890" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_rs2470890", i]){ 
        star_allele  <- "*1H"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs2069526" %in% rownames(genotypes) && "G" %in% genotypes["CYP1A2_rs2069526", i] &&
               "CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i]){ 
        star_allele  <- "*1J"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs12720461" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs12720461", i]){ 
        star_allele  <- "*1K"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs2069514" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs2069514", i] &&
               "CYP1A2_rs35694136" %in% rownames(genotypes) && "-" %in% genotypes["CYP1A2_rs35694136", i] &&
               "CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
               "CYP1A2_rs2470890" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_rs2470890", i]){ 
        star_allele  <- "*1L"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
               "CYP1A2_rs2472304" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs2472304", i]){ 
        star_allele  <- "*1M"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_X3954TG" %in% rownames(genotypes) && "G" %in% genotypes["CYP1A2_X3954TG", i] &&
               "CYP1A2_rs35694136" %in% rownames(genotypes) && "-" %in% genotypes["CYP1A2_rs35694136", i] &&
               "CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i] &&
               "CYP1A2_X2321GC" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_X2321GC", i] &&
               "CYP1A2_rs2470890" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_rs2470890", i] &&
               "CYP1A2_X5521AG" %in% rownames(genotypes) && "G" %in% genotypes["CYP1A2_X5521AG", i]){ 
        star_allele  <- "*1N"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_X733GC" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_X733GC", i]){ 
        star_allele  <- "*1P"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_X2808AC" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_X2808AC", i]){ 
        star_allele  <- "*1Q"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_X367CT" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_X367CT", i]){ 
        star_allele  <- "*1R"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_X3053AG" %in% rownames(genotypes) && "G" %in% genotypes["CYP1A2_X3053AG", i]){ 
        star_allele  <- "*1S"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_X2667TG" %in% rownames(genotypes) && "G" %in% genotypes["CYP1A2_X2667TG", i]){ 
        star_allele  <- "*1T"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_X678CT" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_X678CT", i]){ 
        star_allele  <- "*1U"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs35694136" %in% rownames(genotypes) && "-" %in% genotypes["CYP1A2_rs35694136", i] &&
               "CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i]){ 
        star_allele  <- "*1V"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_X3113AG" %in% rownames(genotypes) && "G" %in% genotypes["CYP1A2_X3113AG", i]){ 
        star_allele  <- "*1W"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs56160784" %in% rownames(genotypes) && "G" %in% genotypes["CYP1A2_rs56160784", i]){ 
        star_allele  <- "*2"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs56276455" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs56276455", i]){ 
        star_allele  <- "*3"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs72547516" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs72547516", i]){ 
        star_allele  <- "*4"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs55889066" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs55889066", i]){ 
        star_allele  <- "*5"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs28399424" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs28399424", i]){ 
        star_allele  <- "*6"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs56107638" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs56107638", i]){ 
        star_allele  <- "*7"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs72547517" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs72547517", i]){ 
        star_allele  <- "*8"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs138652540" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs138652540", i]){ 
        star_allele  <- "*9"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs72547512" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_rs72547512", i]){ 
        star_allele  <- "*10"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs72547513" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs72547513", i]){ 
        star_allele  <- "*11"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_X634AT" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_X634AT", i]){ 
        star_allele  <- "*12"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs35796837" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs35796837", i]){ 
        star_allele  <- "*13"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs45486893" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs45486893", i]){ 
        star_allele  <- "*14"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs2470890" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_rs2470890", i] &&
                "CYP1A2_rs72547511" %in% rownames(genotypes) && "G" %in% genotypes["CYP1A2_rs72547511", i]){ 
        star_allele  <- "*15"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs2470890" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_rs2470890", i] &&
                "CYP1A2_rs72547515" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs72547515", i]){ 
        star_allele  <- "*16"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs2470890" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_rs2470890", i] &&
                "CYP1A2_rs149928755" %in% rownames(genotypes) && "T" %in% genotypes["CYP1A2_rs149928755", i] &&
                "CYP1A2_rs2472304" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs2472304", i] &&
                "CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i]){ 
        star_allele  <- "*17"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs2470890" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_rs2470890", i] &&
                "CYP1A2_X3468AC" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_X3468AC", i]){ 
        star_allele  <- "*18"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs2470890" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_rs2470890", i] &&
                "CYP1A2_X5328GA" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_X5328GA", i]){ 
        star_allele  <- "*19"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs2470890" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_rs2470890", i] &&
                "CYP1A2_rs144148965" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs144148965", i]){ 
        star_allele  <- "*20"
        genotype_match <- TRUE
        points <- 50
      }else if ("CYP1A2_rs2470890" %in% rownames(genotypes) && "C" %in% genotypes["CYP1A2_rs2470890", i] &&
                "CYP1A2_rs143193369" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs143193369", i] &&
                "CYP1A2_rs17861157" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs17861157", i] &&
                "CYP1A2_rs762551" %in% rownames(genotypes) && "A" %in% genotypes["CYP1A2_rs762551", i]){ 
        star_allele  <- "*17"
        genotype_match <- TRUE
        points <- 50
      }else {
        star_allele <- "*1"
        genotype_match <- TRUE
        points <- 50
      }
      if(genotype_match){
        allele_matches[i] <- star_allele
        score[i] <- points
      }
    }
  }
  if (all(allele_matches != "")) {
    # Format diplotype for both alleles
    diplotype <- paste(allele_matches[1], allele_matches[2], sep = "/")
    total_score <- score[1] + score[2]
    
    #save the diplotype, calculate total score and then reset variables 
    diplotypes[row_idx] <- diplotype
    no_match <- 0
    one_match <- NA
    matched_allele <- NA
    total_score <- score[1] + score[2]
  }
  #Handle unmatched sequences if any star alleles couldn't be found from index polymorphisms
  else if((allele_matches[1] != "" & allele_matches[2] == "") | (allele_matches[1] == "" & allele_matches[2] != "")){
    no_match <- 1
    one_match <- which(allele_matches != "")
    matched_allele <- allele_matches[one_match]
    score[one_match] <- score[one_match]
    total_score <- score[one_match]
  }
  #the case where no matches were found 
  else if(all(allele_matches == "")){
    no_match <- 2
    one_match <- NA
    matched_allele <- NA
    total_score <- 0
  }
  return(list(no_match, diplotypes, matched_allele, one_match, score, total_score))
}


#' @title Assign Star Allele Genotypes  
#' @description This function assigns diplotypes based on genotypes dynamically using an allele definition table. Using the SNPs in the input data, this assigns the corresponding diplotype.
#' "-" will be treated as an indel. 
#' Genotype data should have the form "A G", "T T", etc. 
#' @param df A data frame containing genotype data with specific columns. Each row represents a patient, and columns represent genotypes for specific SNPs. 
#' Column names should follow the format "gene_snp", for example, "CYP2D6_rs1065852".
#' @param genes The specific CYP genes, as a list of strings
#' @return A data frame with an additional column `gene_diplotype` indicating the assigned star allele genotype and messages indicating missing SNPs and/or genotypes.
#' @examples
#' find_diplotype(data, c("CYP2D6", "CYP2C9"))
#' @export
find_diplotype <- function(df, genes, CYP1A2_name){
  #keep genes consistent as upper case 
  genes <- toupper(genes)
  #initialize total score 
  total_score <- 0
  for (gene in genes) {
      # Load allele definitions based on the gene
    #implement special behaviour for CYP1A2 depending on if user wants new or old nomenclature for star alleles
      data <- switch(gene,
                     "CYP2B6" = CYP2B6_Allele_def,
                     "CYP2C9" = CYP2C9_Allele_def,
                     "CYP2D6" = CYP2D6_Allele_def, 
                     "CYP2C19" = CYP2C19_Allele_def,
                     "CYP1A2" = if (CYP1A2_name == "new") {
                       CYP1A2_Allele_def_2
                     } else if (CYP1A2_name == "old") {
                       CYP1A2_Allele_def
                     } else {
                       stop("Invalid value for CYP1A2_name. Must be 'new' or 'old'.")
                     }, 
                     "CYP3A4" = CYP3A4_Allele_def, 
                     "CYP3A5" = CYP3A5_Allele_Def)
      #check for missing genotypes in the data
      missing_geno <- find_missing_data(df)
      # Set the column names from the first row (row 6 in the original file from pharmkgb)
      columns <- colnames(data)
      data <- data[-1, ]  # Remove the row with SNP names
      # Extract the star alleles (assuming the first column contains star alleles)
      star_alleles <- data[, 1]
      data <- data[, -1]  # Remove the first column with star alleles
      
      # Filter the input dataframe columns to include only those relevant to the current gene
      gene_snps <- grep(paste0("^", gene, "_"), colnames(df), value = TRUE)
    
      # Extract the rs part from the input column names
      snps_clean <- gsub(".*_X?(rs[0-9]+|X[0-9A-Z]+)", "\\1", gene_snps)

      # Extract base SNP names from snps_clean
      base_snps <- ifelse(grepl("^X[0-9]+[A-Z]+$", snps_clean), snps_clean, gsub("[A-Z]$", "", snps_clean))

      # Match base SNPs with columns in the allele definition data
      matching_snps <- intersect(base_snps, colnames(data))

      # Filter out SNPs not found in the allele definition data
      valid_snps_clean <- snps_clean[base_snps %in% matching_snps]
      # Filter df_gene to only contain SNPs in valid_snps_clean with the gene prefix
      #Avoids processing genes that are not present in allele definition data
      valid_gene_snps <- gene_snps[snps_clean %in% valid_snps_clean]
      df_gene <- df[, valid_gene_snps, drop = FALSE]

      # Ensure there are matching SNPs before proceeding
      if (length(valid_snps_clean) > 0) {
        #Create a new data_filtered dataframe
        data_filtered <- data.frame() 
        
        # Loop through valid_snps_clean and add the corresponding column from data
        for (snp in valid_snps_clean) {
          base_snp <- gsub("[A-Z]$", "", snp)  # Extract base SNP name
          
          # Add the corresponding column from data to data_filtered
          if (nrow(data_filtered) == 0) {
            if (length(matching_snps) == 1) { #in the case where there's only one column, explicitly keep as dataframe
              data_filtered <- data.frame(data[[base_snp]], drop = FALSE, stringsAsFactors = FALSE)
            } else {
              data_filtered <- data.frame(data[[base_snp]], stringsAsFactors = FALSE)
            }
            colnames(data_filtered) <- snp  # Name it according to valid_snps_clean
          } else {
            data_filtered[[snp]] <- data[[base_snp]]
          }
        }
  
        
        #Filter data_filtered again to ensure only star alleles with their functional snps are included in data 
        # Create a vector to store the original row indices
         original_row_indices <- seq_len(nrow(data))
        # Filter rows based on conditions
        row_indices <- apply(data, 1, function(row, idx) {
          # Check if there is at least one non-empty cell in the SNP columns of data_filtered
          non_empty_in_snp <- any(row[colnames(data_filtered)] != "" & !is.na(row[colnames(data_filtered)]))
          # Determine if the row should be kept
          keep_row <- non_empty_in_snp
          return(keep_row)
        }, idx = original_row_indices)
       # Ensure the first row is always included
         row_indices[1] <- TRUE
         data_filtered <- data_filtered[row_indices, ]
 
         #Map back to original star alleles
         filtered_star_alleles <- star_alleles[row_indices]
       
         #fill the empty cells with reference genotype
         data_filled <- fill_empty_cells(data_filtered)
   
        # Initialize a vector to store diplotypes
        diplotypes <- vector("character", nrow(df))
        
        for (row_idx in 1:nrow(df_gene)) {
          # Initialize score vector
          score <- c(0, 0)
       
          #splits the alleles into two parts across the current row of df_gene, goes from "A G" to "A" "G"
          genotypes <- t(sapply(df_gene[row_idx, ], function(x) unlist(strsplit(x, " "))))
          allele_matches <- vector("character", 2)
          genotype_match <- FALSE
          # Initialize variables relevant for CYP2D6
          one_match <- NA
          matched_allele <- NA
          no_match <- 0
          #special matching for tier 1 variants for cyp2d6 in another function 
          if (gene == "CYP2D6"){
            cnv_x9_col <- grep("CNVx9", colnames(df), value = TRUE)
            cnv_int6_col <- grep("CNVInt6", colnames(df), value = TRUE)
            cnv_5utr_col <- grep("CNV5UTR", colnames(df), value = TRUE)
            
            # Access the CNV_x9, CNV_5utr and CNV_Int6 columns to get copy numbers 
            cnv_x9 <- as.character(df[row_idx, cnv_x9_col])
            cnv_int6 <- as.character(df[row_idx, cnv_int6_col])
            cnv_5utr <- as.character(df[row_idx, cnv_5utr_col])

            #call helper function to assign cyp2d6 star alleles
            cyp2d6_result <- process_cyp2d6(genotypes, row_idx, diplotypes, cnv_x9, cnv_int6, cnv_5utr, score, total_score) 
            no_match <- cyp2d6_result[[1]]
            diplotypes <- cyp2d6_result[[2]]
            matched_allele <- cyp2d6_result[[3]]
            one_match <- cyp2d6_result[[4]]
            score <- cyp2d6_result[[5]]
            total_score <- cyp2d6_result[[6]]
          }
          
          if(gene == "CYP1A2"){
            cyp1a2_result <- process_cyp1a2(genotypes, row_idx, diplotypes, score, total_score, CYP1A2_name)
            no_match <- cyp1a2_result[[1]]
            diplotypes <- cyp1a2_result[[2]]
            matched_allele <- cyp1a2_result[[3]]
            one_match <- cyp1a2_result[[4]]
            score <- cyp1a2_result[[5]]
            total_score <- cyp1a2_result[[6]]
          }
       
          #case where gene is not cyp2d6 or both haplotypes didn't match in process_cyp2d6 or process_cyp1a2
          if((no_match >= 1 && gene == "CYP2D6") || 
             (no_match >= 1 && gene == "CYP1A2") || 
             gene != "CYP2D6" &&
             gene != "CYP1A2"){ #(no_match >=1 && (gene == "CYP2D6" | gene == "CYP1A2")) | (gene != "CYP2D6" | gene != "CYP1A2")
      
            for (i in 1:2) {  # Loop over the genotypes
              genotype_match <- TRUE  # Assume a match until proven otherwise
              
              for (j in 1:nrow(data_filled)) {  # Iterate through each row of data_filled
                genotype_match <- TRUE  # Reset to TRUE for each row check
              
                for (k in 1:length(valid_snps_clean)) {  # Check each SNP in the input dataframe
                  allele <- strsplit(df_gene[[k]][row_idx], " ")[[1]][i]  # Get the allele for the current SNP and genotype
                  snp_data <- data_filled[j, k]  # SNP data for the current row
                  base_snp <- gsub("[A-Z]$", "", valid_snps_clean[k])# Strip suffix from valid_snps_clean[k] for comparison
                  
                  
                  #special case for cyp2d6 and rs5030865 
                  if (gene == "CYP2D6" && base_snp == "rs5030865") {
                    allele_T <- df_gene[["CYP2D6_rs5030865T"]][row_idx]
                    allele_A <- df_gene[["CYP2D6_rs5030865A"]][row_idx]
                    
                    if (allele_T == "C C" && allele_A == "C C") {
                      # Both are homozygous C C, it doesn't matter which one is used
                      allele <- "C"
                    } else if (allele_T != "C C" && allele_A == "C C" ) {# && allele_A == "C C"
                      allele <- strsplit(allele_T, " ")[[1]][i]  # Use rs5030865T if it's not C C
                    } else if (allele_A != "C C" && allele_T == "C C"){
                      allele <- strsplit(allele_A, " ")[[1]][i]  # Otherwise, use rs5030865A
                    }else if(allele_T != "C C" && allele_A != "C C"){ 
                      #Allele A and T are heterozygous, combine non-C alleles
                      alleles_T <- strsplit(allele_T, " ")[[1]]
                      alleles_A <- strsplit(allele_A, " ")[[1]]
                      #Find non-C alleles
                      non_C_T <- alleles_T[alleles_T != "C"]
                      non_C_A <- alleles_A[alleles_A != "C"]
                      
                      #Combine non-C alleles
                      if (length(non_C_T) > 0 && length(non_C_A) > 0) {
                        allele_combo <- paste(non_C_T[1], non_C_A[1])
                        allele <- strsplit(allele_combo, " ")[[1]][i] 
                      } else {
                        allele <- NA  # Handle cases where no non-C alleles are found
                      }
                    }
                    else{ 
                      allele <- NA
                    }
                  }
                  #checking iupac codes found in allele definition tables for matches in input data 
                  if (!(allele %in% unlist(strsplit(snp_data, " "))) &&  # Direct match check
                      !(snp_data == "R" && allele %in% c("A", "G")) &&  # R matches A or G
                      !(snp_data == "Y" && allele %in% c("C", "T")) &&  # Y matches C or T
                      !(snp_data == "M" && allele %in% c("A", "C")) &&  # M matches A or C
                      !(snp_data == "S" && allele %in% c("G", "C")) &&  # S matches G or C
                      !(snp_data == "K" && allele %in% c("G", "T")) &&  # K matches G or T
                      !(snp_data == "W" && allele %in% c("A", "T")) &&  # W matches A or T
                      !(snp_data == "H" && allele %in% c("A", "T", "C")) &&
                      !(snp_data == "B" && allele %in% c("C", "T", "G")) &&
                      !(snp_data == "D" && allele %in% c("A", "T", "G")) &&
                      !(snp_data == "V" && allele %in% c("C", "A", "G")) &&
                      !(grepl("del", snp_data) && allele == "-")) {  # "del" matches "-"
                    genotype_match <- FALSE
                    break  # Exit the SNP loop early if a mismatch is found
                  }
                }
                
                if (genotype_match) {
                  allele_matches[i] <- filtered_star_alleles[j] 
                  break  # Exit the row loop early if a matching row is found
                }
              }
              #if no match found, save message
              if (!genotype_match) {
                diplotypes[row_idx] <- "No matching star alleles found"
                break
              }
            }
            
            if (genotype_match) {
              #if matched allele found from process_cyp2d6, replace the allele_matches at that index with it
              if(!is.na(matched_allele)){
                allele_matches[one_match] <- matched_allele
                score[one_match] <- 50
              }
              #format star alleles properly (ex.*1/*3)
              diplotype <- paste(allele_matches[1], allele_matches[2], sep = "/") 
              
              #linkage disequilibrium for CYP2B6
              if(gene == "CYP2B6" && (diplotype == "*4/*9" | diplotype == "*9/*4")){  
                diplotype <- "*1/*6"
                total_score <- 50
              }
              
              #Check for specific CYP2D6 star alleles given CNV information 
              if(gene == "CYP2D6"){
                #linkage disequilibrium case for cyp2d6
                if(diplotype == "*4/*10" | diplotype == "*10/*4"){
                  diplotype <- "*1/*4"
                  total_score <- 50
                }
                #cnv checks if all info availble 
                if((length(cnv_x9) > 0 && !is.na(cnv_x9)) && 
                   (length(cnv_int6) > 0 && !is.na(cnv_int6)) && 
                   (length(cnv_5utr) > 0 && !is.na(cnv_5utr))){
                  #whole gene deletion
                  if(all(cnv_x9 == "0 0" & cnv_int6 == "0 0" & cnv_5utr == "0 0")){ 
                    diplotype <- "*5/*5"
                    total_score <- 50
                  }
                  #other special cases for cyp2d6
                  if(all(cnv_x9 == "2 2" & cnv_int6 == "4 4" & cnv_5utr == "4 4")){
                    diplotype <- "*36+*10/*36+*10"
                    total_score <- 50
                  }
                  if(all(cnv_x9 == "2 2" & cnv_int6 == "3 3" & cnv_5utr == "3 3")){
                    diplotype <- "*10/*36+*10"
                    total_score <- 50
                  }
                  if(all(allele_matches[1] == "4" & allele_matches[2] == "4" & cnv_x9 == "2 2" & cnv_int6 == "2 2" & cnv_5utr == "3 3")){
                    diplotype <- "*4/*68+*4"
                    total_score <- 50
                  }
                }
                #if only data from exon 9 and intron 6 available
                else if(length(cnv_x9) > 0 && !is.na(cnv_x9) && 
                        length(cnv_int6) > 0 && !is.na(cnv_int6) && 
                        (length(cnv_5utr) <= 0 || is.na(cnv_5utr))){ 
                  #whole gene deletion
                  if(all(cnv_x9 == "0 0" & cnv_int6 == "0 0")){ #& CNV_flag
                    diplotype <- "*5/*5"
                    total_score <- 50
                  }
                  #other special cases for cyp2d6
                  if(all(cnv_x9 == "2 2" & cnv_int6 == "4 4")){
                    diplotype <- "*36+*10/*36+*10"
                    total_score <- 50
                  }
                  if(all(cnv_x9 == "2 2" & cnv_int6 == "3 3")){
                    diplotype <- "*10/*36+*10"
                    total_score <- 50
                  }
                  if(all(allele_matches[1] == "4" & allele_matches[2] == "4" & cnv_x9 == "2 2" & cnv_int6 == "2 2")){
                    diplotype <- "*4/*68+*4"
                    total_score <- 50
                  }
                }
                #if only cnv_x9 CNV data available
                else if(length(cnv_x9) > 0 && !is.na(cnv_x9) && 
                        (length(cnv_int6) <= 0 || is.na(cnv_int6)) && 
                        (length(cnv_5utr) <= 0 || is.na(cnv_5utr))){
                  #whole gene deletion check
                  if(all(cnv_x9 == "0 0")){ 
                    diplotype <- "*5/*5"
                    total_score <- 50
                  }
                }
              }
              if(gene == "CYP1A2"){ #case of *6 for CYP1A2 based on presence of index polymorphism
                if(diplotype == "*1A/*1A" & (df$CYP1A2_rs28399424[row_idx] == "C T" | df$CYP1A2_rs28399424[row_idx] == "T C")){
                  diplotype <- "*1A/*6"
                }
                if(diplotype == "*1A/*1A" & df$CYP1A2_rs28399424[row_idx] == "T T"){
                  diplotype <- "*6/*6"
                }
              }
              #add the diplotype to the row to save 
              diplotypes[row_idx] <- diplotype
            }
          }
          
          if(total_score == 0){
            for (i in 1:2) {
              if (!is.na(allele_matches[i])) {
                matched_allele <- allele_matches[i]
                
                # Find the row index in data_filled corresponding to the matched allele
                matched_row_idx <- which(filtered_star_alleles == matched_allele)
                
                if (length(matched_row_idx) > 0) {
                  if (matched_allele == "*1" | matched_allele == "*1A") {
                    # For *1 allele, score based on the number of SNP columns that differ from first row (i.e 0 difference for row 1)
                    score[i] <- 0 
                  } else {
                    # For other alleles, calculate the score based on SNP differences
                    score[i] <- sum(data_filled[matched_row_idx, ] != data_filled[1, ])
                  }
                }
                #save score if not 0
                if(score[i] != 0){
                  score[i] <- score[i]
                }
              }
            }
            # Calculate total score for both alleles
            total_score <- sum(score)
          }
          else{ #store total score 
            total_score <- total_score
          }
         
          #Dynamically create the diplotype column name based on the gene provided by the user
          diplotype_column_name <- paste0(gene, "_diplotype")
          df[[diplotype_column_name]] <- diplotypes
          } 
        }
      else{#check if data contains information for the input gene
        warning(paste("No matching SNPs found for gene:", gene))
        }
      }
  return(list(df, missing_geno, total_score)) 
  }

#' @title Star Allele to Activity Score Conversion
#'
#' @description Convert star allele to activity score 
#' @param star_alleles The star alleles you want to convert, input as a dataframe with one column containing the star alleles. For CYP2D6, a column with CNV_x9 should be included in input dataframe. If there isn't CNV provided, it will default to 2.
#' @param gene The specific CYP gene, as a string
#' @return The associated activity score, appended to the original dataframe input 
#' @examples 
#' score1 <- star_to_activity(df, "CYP2C9")
#' @export
star_to_activity <- function(df, gene) { 
  gene <- toupper(gene)
  data <- switch(gene,
                 "CYP2C9" = CYP2C9_data,
                 "CYP2C19" = CYP2C19_data,
                 "CYP2D6" = CYP2D6_data,
                 "CYP2B6" = CYP2B6_data,
                 "CYP3A4" = NULL,  # Allow CYP3A4 to proceed with NA
                 "CYP1A2" = NULL, # Allow CYP1A2 to proceed with NA
                 "CYP3A5" = NULL, # Allow CYP3A5 to proceed with NA
                 stop("Unsupported gene provided")  # Default case for unsupported genes
  )
  
  # Initialize empty lists to store activity scores
  activity_score1_list <- list()
  activity_score2_list <- list()
  duplication_list <- list() 
  comment_list <- list()
  
  # Check for genes that don't have activity scores per PharmGKB and assign NA to activity scores
  if (gene %in% c("CYP1A2", "CYP3A4", "CYP2C19", "CYP2B6", "CYP3A5")) {
    activity_score1_list <- rep(NA, nrow(df))
    comment_list <- rep("No activity scores available", nrow(df))
    
  } else if (gene == "CYP2D6") { #special use for cyp2d6
    for (i in 1:nrow(df)) {
      diplotype_col <- paste0(gene, "_diplotype")
      diplotype <- df[i, diplotype_col]
      
      if (diplotype == " " || is.na(diplotype) || is.null(diplotype)) {
        activity_score1_list[[i]] <- "Diplotype input missing"
        activity_score2_list[[i]] <- "Diplotype input missing"
        duplication_list[[i]] <- "No"
        comment_list[[i]] <- NA
      } else {
        cnv_x9_col <- grep("CNVx9", colnames(df), value = TRUE)
        if (length(cnv_x9_col) == 0) {
          cnv_x9_col <- grep("CNV", colnames(df), value = TRUE)
        }
        if (length(cnv_x9_col) == 0) {
          cnv_x9 <- 2
        } else {
          cnv_x9_value <- df[i, cnv_x9_col]
          cnv_x9 <- ifelse(is.na(cnv_x9_value), 2, as.numeric(strsplit(cnv_x9_value, " ")[[1]][1]))
        }
        
        duplication_status <- ifelse(cnv_x9 > 2, "Yes", "No")
        duplication_list[[i]] <- duplication_status
        
        alleles <- strsplit(diplotype, "/")[[1]]
        
        # Check if one of the alleles is in the form *X+*Y (e.g., *36+*10)
        plus_alleles <- grepl("\\+", alleles)
        
        if (any(plus_alleles)) {
          # Expand the *36+*10 into separate alleles
          expanded_alleles <- unlist(strsplit(alleles[plus_alleles], "\\+"))
          # Combine expanded alleles with the other non-plus allele
          all_alleles <- c(alleles[!plus_alleles], expanded_alleles)
          
          # Get activity scores for all alleles
          all_scores <- sapply(all_alleles, function(allele) {
            row_index <- which(data[, 1] == allele)
            if (length(row_index) == 0) {
              return(NA)
            }
            as.numeric(data[row_index, 2])
          })
          
          if (any(is.na(all_scores))) {
            activity_score1_list[[i]] <- "Allele not found"
            activity_score2_list[[i]] <- "Allele not found"
          } else {
            # Sum all scores regardless of CNV (because we already expanded the duplication)
            total_score <- sum(all_scores)
            activity_score1_list[[i]] <- total_score
            activity_score2_list[[i]] <- NA  # Usually, 3-allele diplotypes don't need a second AS
          }
        } else {
          # Handle standard diplotype (2 alleles)
          allele1 <- alleles[1]
          allele2 <- alleles[2]
  
          row_index1 <- which(data[, 1] == allele1)
          row_index2 <- which(data[, 1] == allele2)
          
          if (length(row_index1) == 0 || length(row_index2) == 0) {
            activity_score1_list[[i]] <- "Allele not found"
            activity_score2_list[[i]] <- "Allele not found"
          } else {
            score1 <- as.numeric(data[row_index1, 2])
            score2 <- as.numeric(data[row_index2, 2])
    
            if (!is.na(cnv_x9) && cnv_x9 > 2) { 
              activity_score1 <- score1 * (cnv_x9 - 1) + score2
              activity_score2 <- score2 * (cnv_x9 - 1) + score1
            } else {
              activity_score1 <- score1 + score2
              activity_score2 <- NA
            }
            
            activity_score1_list[[i]] <- activity_score1
            activity_score2_list[[i]] <- ifelse(duplication_status == "Yes", activity_score2, NA)
          }
        }
        comment_list[[i]] <- NA
      }
    }
    df$CYP2D6_Duplication <- unlist(duplication_list)
  }
else {
    # For genes other than CYP2D6, handle activity score 1 only
    for (i in 1:nrow(df)) {
      diplotype_col <- paste0(gene, "_diplotype")
      diplotype <- df[i, diplotype_col]
     #check if diplotype input is present
      if (diplotype == " " || is.na(diplotype) || is.null(diplotype)) {  
        activity_score1_list[[i]] <- "Diplotype input missing"
        comment_list[[i]] <- ""
      } else {
        row_index <- which(data[, 1] == diplotype)
        #check reverse diplotype if not found the first time
        if (length(row_index) == 0) {
          reversed_diplotype <- paste(rev(unlist(strsplit(diplotype, "/"))), collapse = "/")
          row_index <- which(data[, 1] == reversed_diplotype)
        }
        if (length(row_index) == 0) {  
          activity_score1_list[[i]] <- "Diplotype not found"
        } else {  
          score <- as.numeric(data[row_index, 2]) 
          activity_score1_list[[i]] <- score
        }
        comment_list[[i]] <- ""
      }
    }
  }
  
  # Append the activity scores and comments to the original dataframe
  df[[paste0(gene, "_AS_1")]] <- unlist(activity_score1_list)
  
  if (gene == "CYP2D6") {
    df[[paste0(gene, "_AS_2")]] <- unlist(activity_score2_list)
  }
  
  df$Comment <- unlist(comment_list)
  return(df)
}


#' @title Star Allele to Phenotype Conversion
#'
#' @description Convert star allele to metabolizer phenotype based on activity scores and external csv files from pharmgkb
#' @param star_alleles The star alleles you want to convert, input as a dataframe with one column containing star alleles. For CYP2D6, activity score columns need to be provided (i.e output from star_to_activity()).
#' @param gene The specific CYP gene, as a string
#' @return The associated metabolizer phenotype, appended to the original dataframe input
#' @examples 
#' pheno1 <- star_to_pheno(df, "CYP2C9")
#' @export
star_to_pheno <- function(df, gene) {
  gene <- toupper(gene)
  # Check if the dataframe is empty
  if (nrow(df) == 0) {
    warning("The dataframe is empty. No processing will be done.")
    return(df)
  }
  
  # Check if gene has phenotype data
  if (gene %in% c("CYP1A2", "CYP3A4")) {
    # Handle unsupported genes by adding NA for phenotype and a specific comment
    df[[paste0(gene, "_Phenotype")]] <- rep(NA, nrow(df))
    if (!"Comment" %in% colnames(df)) {
      df$Comment <- rep("No phenotype available", nrow(df))
    } else {
      df$Comment <- paste(df$Comment, "No phenotype available", sep = "; ")
    }
    return(df)
  }
  
  # Define phenotype data based on gene
  data <- switch(gene,
                 "CYP2C9" = CYP2C9_data,
                 "CYP2C19" = CYP2C19_data,
                 "CYP2B6" = CYP2B6_data,
                 "CYP2D6" = CYP2D6_data,
                 "CYP3A5" = CYP3A5_data,
                 NULL)  # Default case for unsupported genes
  
  # Initialize lists to store phenotypes and comments
  pheno_list_1 <- list()
  if ("Comment" %in% colnames(df)) {
    comment_list <- as.list(df$Comment)
  } else {
    comment_list <- rep("", nrow(df))
  }
  
  # Create the column name for the gene_diplotype dynamically
  diplotype_col <- paste0(gene, "_diplotype")
  star_alleles <- df[, diplotype_col, drop = FALSE]
  
  if (gene == "CYP2D6") {
    pheno_list_2 <- list()
    activity_score1 <- df[,"CYP2D6_AS_1", drop = FALSE]
    activity_score2 <- df[,"CYP2D6_AS_2", drop = FALSE]
    
    # Loop through activity scores to determine phenotype
    for (i in 1:nrow(activity_score1)) {
      score1 <- activity_score1[i, "CYP2D6_AS_1"]
      score2 <- activity_score2[i, "CYP2D6_AS_2"]
      diplotype_val <- df[i, diplotype_col]

      # Check for invalid or unmatched diplotype
      if (diplotype_val %in% c("no matching star alleles found", NA, "", " ")) {
        pheno_list_1[[i]] <- "Diplotype not found"
        pheno_list_2[[i]] <- NA
        comment_list[[i]] <- paste(comment_list[[i]], "No metabolizer status assigned due to unmatched diplotype", sep = "; ")
        next
      }
      
      # Check if scores are numeric
      if (!suppressWarnings(!is.na(as.numeric(score1)))) {
        pheno_list_1[[i]] <- "Invalid activity score"
        pheno_list_2[[i]] <- NA
        comment_list[[i]] <- paste(comment_list[[i]], "Non-numeric activity score; no metabolizer status assigned", sep = "; ")
        next
      }
      
      #assign phenotype from activity score 1
      score1 <- as.numeric(score1)
      pheno1 <- if (score1 >= 2.5) "Ultra Rapid" else if (score1 >= 1.25 && score1 <= 2.25) "Normal" else if (score1 >= 0.25 && score1 <= 1) "Intermediate" else if (score1 == 0) "Poor" else "Unknown"
      pheno_list_1[[i]] <- pheno1
      
      # Check if Activity_Score_2 exists and is not NA
      if (!is.na(score2) && score2 != "") {
        score2 <- as.numeric(score2)
        pheno2 <- if (score2 >= 2.5) "Ultra Rapid" else if (score2 >= 1.25 && score2 <= 2.25) "Normal" else if (score2 >= 0.25 && score2 <= 1) "Intermediate" else if (score2 == 0) "Poor" else "Unknown"
        pheno_list_2[[i]] <- pheno2
        
        # Handle ambiguous CNV cases
        comment_list[[i]] <- ifelse(pheno1 == pheno2, paste(comment_list[[i]], "Ambiguous CNV resulting in same phenotype", sep = "; "), paste(comment_list[[i]], "Ambiguous CNV, two potential phenotypes", sep = "; "))
      } else {
        pheno_list_2[[i]] <- NA
      }
    }
    
  } else {
    # For non-CYP2D6 genes, use external CSV to determine phenotype
    for (i in 1:nrow(star_alleles)) {
      diplotype <- star_alleles[i, diplotype_col]
      
      if (diplotype == " " || is.na(diplotype) || is.null(diplotype)) {  # Check if any input is missing
        pheno_list_1[[i]] <- "Diplotype input missing"
      } else {
        # Check in the external CSV file for the diplotype
        row_index <- which(data[,1] == diplotype)
        if (length(row_index) == 0) {  # If not found, check reverse diplotype
          reversed_diplotype <- paste(rev(unlist(strsplit(diplotype, "/"))), collapse = "/")
          row_index <- which(data[,1] == reversed_diplotype)
        }
        if (length(row_index) == 0) {  # Handle diplotype not found
          pheno_list_1[[i]] <- "Diplotype not found"
        } else {  # Diplotype found, assign phenotype
          pheno_list_1[[i]] <- data[row_index, 3]
        }
      }
    }
  }
  
  # Append the phenotype and comment columns to the original dataframe
  df[[paste0(gene, "_Phenotype")]] <- unlist(pheno_list_1)
  df$Comment <- unlist(comment_list)
  
  # For CYP2D6, also add the second phenotype column
  if (gene == "CYP2D6") {
    df[[paste0(gene, "_Phenotype_2")]] <- unlist(pheno_list_2)
  }
  
  return(df)
}


#' @title Filter For No Matching Star Alleles
#'
#' @description This function will filter the dataframe to only include rows where no star allele matches were found. It's used within the assign_diplotype function.
#' @param df The dataframe that contains the genotypes and the diplotype column
#' @param genes The specific CYP genes, as a list of strings
#' @return A dataframe only containing rows that have "no matching star alleles found" as values in the diplotype column. 
#' @examples 
#' filtered_df <- filter_no_match(df, c("CYP2C9", "CYP2B6"))
#' filtered_df2 <- filter_no_match(df, "CYP2C19")
#' @export
filter_no_match <- function(df, gene) {
  gene <- toupper(gene)
  #Create the column name dynamically based on the gene name
  diplotype_col <- paste0(gene, "_diplotype")
  
  #Filter the dataframe where the diplotype column has no star allele match
  filtered_df <- df[df[[diplotype_col]] == "No matching star alleles found", ]
  
  return(filtered_df)
}

#' @title Generate phased haplotype combinations
#' 
#' @description This function will produce all possible combinations with the genotypes provided.
#' @param snp_df The dataframe that contains the genotypes for the interested snps. As a haplotype, ex. "A" at rs1065852
#' @param snp_cols The names of the snps as the columns for the dataframe
#' @return A matrix array containing all possible haplotype combinations. 
#' @export
generate_phased_combinations <- function(snp_df, snp_cols) {
  haplotype_combinations <- list()
  
  # Generate phased haplotypes for each SNP column
  for (snp in snp_cols) {
    # Ensure the SNP values are characters
    snp_values <- as.character(snp_df[[snp]])
    # Split each SNP value into two alleles
    alleles <- strsplit(snp_values, " ")
    
    # For each SNP, create the two possible phased combinations in "A A" format
    haplotypes <- sapply(alleles, function(x) c(paste(x[1], x[2]), paste(x[2], x[1])))
    
    haplotype_combinations[[snp]] <- unique(haplotypes)  # Ensure unique combinations
  }
  
  # Generate all possible phased haplotype combinations
  all_combinations_df <- expand.grid(haplotype_combinations)
  
  # Convert the data.frame to a matrix of character vectors
  all_combinations_matrix <- as.matrix(all_combinations_df)
  
  return(all_combinations_matrix)
}

#' @title Function to format genotypes
#' @description Checks if genotypes contain a space between the alleles, if not, adds one (i.e AA to A A)
#' @param genotype the genotype from the current column and row
#' @return the newly formatted genotype 
#' @export
format_genotype <- function(genotype) {
  if (length(genotype) > 1) {
    stop("Unexpected vector input in format_genotype: ", paste(genotype, collapse = ", "))
  }
  
  # Handle NA or empty values
  if (is.na(genotype) || genotype == "") {
    return(genotype)
  }
  
  # Special case for genotypes with "-" or similar non-standard patterns
  if (grepl("-", genotype)) {
    if (nchar(genotype) > 1 && grepl("^-", genotype)) {
      # Case: "-TCT"
      return(paste0("- ", substr(genotype, 2, nchar(genotype))))
    } else if (nchar(genotype) > 1 && grepl("-$", genotype)) {
      # Case: "TCT-"
      return(paste0(substr(genotype, 1, nchar(genotype) - 1), " -"))
    }
  }
  
  # If no space and genotype length is even, split into two equal parts
  if (!grepl(" ", genotype) && nchar(genotype) %% 2 == 0) {
    half_length <- nchar(genotype) / 2
    return(paste0(substr(genotype, 1, half_length), " ", substr(genotype, half_length + 1, nchar(genotype))))
  }
  
  # If already formatted or cannot be split evenly, return as-is
  return(genotype)
}

#' @title Assigns more common cyp2d6 diplotype if tie in scoring
#' @description The function will check the diplotype and alternate diplotype call and ascertain the more likely diplotype.
#' @param df The dataframe that contains the patientID, the genotypes, and the diplotypes from assign_diplotype()
#' @return A dataframe with updated diplotype calls
#' @examples 
#' result <- adjust_diplotype(df)
#' 
#' @export
adjust_diplotype <- function(final_results) {
  if("CYP2D6_alternate_diplotype" %in% colnames(final_results)){
    # Define the pairs of diplotype replacements
    replace_pairs <- list(
      "*10/*4" = c("*1/*4", "*2/*4"),
      "*1/*3"  = c("*2/*3", "*10/*3"),
      "*1/*6"  = c("*2/*6", "*10/*6"),
      "*1/*9"  = c("*2/*9", "*10/*9"),
      "*1/*7"  = "*2/*7",
      "*10/*68+*4" = "*2/*68+*4"
    )
    
    # Loop through the pairs and apply replacements
    for (d in names(replace_pairs)) {
      for (a in replace_pairs[[d]]) {
        # Identify rows where CYP2D6_diplotype and CYP2D6_alternate_diplotype match
        idx <- with(final_results, CYP2D6_diplotype == d & CYP2D6_alternate_diplotype == a)
        # Swap the diplotype and set alternate_diplotype to NA
        final_results$CYP2D6_diplotype[idx] <- a
        final_results$CYP2D6_alternate_diplotype[idx] <- NA
      }
    }
    
    # Additional condition: If diplotype is *1/*4 and alternate is *10/*4, set alternate to NA
    idx_alt <- with(final_results, (CYP2D6_diplotype == "*1/*4" | CYP2D6_diplotype == "*2/*4") & CYP2D6_alternate_diplotype == "*10/*4")
    final_results$CYP2D6_alternate_diplotype[idx_alt] <- NA
    
    idx_alt <- with(final_results, (CYP2D6_diplotype == "*10/*14") & CYP2D6_alternate_diplotype == "*14/*2")
    final_results$CYP2D6_alternate_diplotype[idx_alt] <- NA
  }

  # Explicit checks for *4 and *10 assignments only if the diplotype is *10/*4 or *4/*10
  if("CYP2D6_rs1065852" %in% colnames(final_results) & "CYP2D6_rs3892097" %in% colnames(final_results)) {
    for (i in seq_len(nrow(final_results))) {
      if (final_results$CYP2D6_diplotype[i] %in% c("*10/*4", "*4/*10")) {
        rs1065852 <- final_results$CYP2D6_rs1065852[i]
        rs3892097 <- final_results$CYP2D6_rs3892097[i]
        
        if ((rs1065852 %in% c("AG", "GA")) & (rs3892097 %in% c("CT", "TC"))) {
          final_results$CYP2D6_diplotype[i] <- "*1/*4"
        } else if (rs1065852 == "AA" & (rs3892097 %in% c("CT", "TC"))) {
          final_results$CYP2D6_diplotype[i] <- "*4/*10"
        } else if (rs1065852 == "AA" & rs3892097 == "TT") {
          final_results$CYP2D6_diplotype[i] <- "*4/*4"
        }
      }
    }
  }

  return(final_results)
}

#' @title Reorder diplotypes
#' @description The function will reorder the diplotypes to make sure the format is small number/big number
#' @param diplotype The current diplotype to check 
#' @return A reformatted diplotype 
#' @export
reorder_diplotype <- function(diplotype) {
  # Function to handle alleles with '+'
  handle_plus_allele <- function(allele) {
    if (grepl("\\+", allele)) {
      parts <- strsplit(allele, "\\+")[[1]]
      # Sort parts based on numeric value, extracting the numeric part
      sorted_parts <- sort(parts, method = "radix", decreasing = TRUE)
      return(paste(sorted_parts, collapse = "+"))
    }
    return(allele)
  }
  
  # Split the diplotype on '/'
  alleles <- strsplit(diplotype, "/")[[1]]
  
  # Handle '+' and sort alleles numerically based on the numeric part
  sorted_alleles <- lapply(alleles, function(x) {
    # Extract numeric part and the allele (e.g., *4, *1)
    numeric_part <- as.numeric(gsub("\\D", "", x))  # Extract numeric part
    return(list(allele = x, numeric_part = numeric_part))  # Return both allele and numeric part
  })
  
  # Sort the alleles based on the numeric part
  sorted_alleles <- sorted_alleles[order(sapply(sorted_alleles, function(x) x$numeric_part))]
  
  # Reassemble the sorted alleles
  sorted_alleles <- sapply(sorted_alleles, function(x) x$allele)
  
  # Return the reordered diplotype
  paste(sorted_alleles, collapse = "/")
}

#' @title Assigns Star Allele Values, phased and unphased approach 
#' @description The function will first try to assign star alleles using a phased approach, any genotypes left without a matching diplotype will then undergo an unphased approach for assigning diplotypes.
#' @param df The dataframe that contains the patientID and the genotypes to convert
#' @param genes The specific CYP genes, as a list of strings
#' @param phased True/False to indicate if input data is phased or unphased 
#' @return A dataframe with added columns containing the corresponding diplotypes for the genotypes 
#' @examples 
#' result <- assign_diplotype(df, c("CYP2C9", "CYP2B6"), phased = FALSE)
#' result2 <- assign_diplotype(df, "CYP2C19", phased = TRUE)
#' @export
assign_diplotype <- function(df, genes, phased = FALSE, CYP1A2_name = "new") {
  # Convert genes to uppercase to ensure consistency
  genes <- toupper(genes)

  # Identify SNP columns (containing "rs" or "CNV" in their names) to format genotypes, X is included for the special case of *18 in CYP1A2 which does not have rsID
  snp_columns <- grep("rs|CNV|X", names(df), value = TRUE)
  df[snp_columns] <- apply(df[snp_columns], 2, function(col) {
    sapply(col, format_genotype)
  })

  # Initialize a list to store results for each gene
  all_gene_results <- vector("list", nrow(df))

  # Loop through all input genes
  for (gene in genes) {
    # Prepare columns for SNPs, get all relevant snp columns for the specified gene
    snp_cols <- grep(paste0(gene, "_(rs|CNV|X)"), names(df), value = TRUE)
    diplotype_col <- paste0(gene, "_diplotype")
    alternate_col <- paste0(gene, "_alternate_diplotype")

    # Initialize diplotype columns if not already in df
    if (!diplotype_col %in% names(df)) df[[diplotype_col]] <- NA
    if (!alternate_col %in% names(df)) df[[alternate_col]] <- vector("list", nrow(df))

    # Loop through each row in the dataframe
    for (row_idx in 1:nrow(df)) {
      # Check if all genotypes for this gene are NA, if yes assign NA and skip to next row
      if (all(is.na(df[row_idx, snp_cols, drop = FALSE]))) {
        df[row_idx, diplotype_col] <- NA
        df[[alternate_col]][[row_idx]] <- NA
        all_gene_results[[row_idx]] <- df[row_idx, ]
        next
      }
      # Identify the patient ID column dynamically
      patient_id_column <- names(df)[grep("ID", names(df), ignore.case = TRUE)][1]

      # Check if the patient ID column was found
      if (is.na(patient_id_column)) {
        stop("No patient ID column found. Ensure the dataframe contains a column with 'ID' in its name.")
      }

      # Include the patient ID column in the list of columns to select
      columns_to_select <- unique(c(patient_id_column, snp_cols))

      # Subset the dataframe
      snp_df <- df[row_idx, columns_to_select, drop = FALSE]
      patient_id <- snp_df[[patient_id_column]]
      snp_df <- snp_df[, !grepl("\\.1$", names(snp_df))]

      # Generate haplotype combinations or use snp_df directly
      if (phased) {
        haplotype_combinations <- snp_df
      } else {
        haplotype_combinations <- generate_phased_combinations(snp_df, snp_cols)
      }

      # Process haplotypes
      match_results <- data.frame(Diplotype = character(), Score = numeric())
      for (haplotype_idx in 1:nrow(haplotype_combinations)) {
        temp_df <- snp_df
        if (!phased) {
          for (snp_idx in 1:length(snp_cols)) {
            haplotype <- haplotype_combinations[haplotype_idx, snp_idx]
            temp_df[[snp_cols[snp_idx]]] <- haplotype
          }
        }

        if (gene == "CYP2D6") {
          # Define the original expected CNV columns
          cnv_cols <- c("CYP2D6_CNV5UTR", "CYP2D6_CNVInt6", "CYP2D6_CNVx9")

          # Check for original CNV columns first
          available_cnv_cols <- cnv_cols[cnv_cols %in% colnames(df)]

          # If no original CNV columns are found, check for CYP2D6_CNV
          if (length(available_cnv_cols) == 0 && "CYP2D6_CNV" %in% colnames(df)) {
            available_cnv_cols <- "CYP2D6_CNV"
          }

          if (length(available_cnv_cols) == 0) {
            # No CNV columns available
            warning("No CNV columns found in the dataframe. Proceeding without CNVs.")
          } else {

            # Extract CNV data for the matching patient
            cnv_data <- df[df[[patient_id_column]] == patient_id, available_cnv_cols, drop = FALSE]

            #format cnv data to match expected input in code
            if (nrow(cnv_data) > 0) {
              # Correct single-copy values (e.g., "2" becomes "2 2")
              cnv_data[] <- lapply(cnv_data, function(col) {
                ifelse(nchar(col) == 1, paste(col, col), col)
              })
             
              # Rename `CYP2D6_CNV` to `CYP2D6_CNVx9` if it's the only available column
              if (length(available_cnv_cols) == 1 && available_cnv_cols == "CYP2D6_CNV" && 
                  !"CYP2D6_CNVx9" %in% colnames(cnv_data)) {
                colnames(cnv_data) <- "CYP2D6_CNVx9"
              }
              
              # Ensure the updated formatted cnv cols in `cnv_data` replaces the ones in `temp_df`
              if ("CYP2D6_CNVx9" %in% colnames(temp_df) && "CYP2D6_CNVx9" %in% colnames(cnv_data)) {
                temp_df$CYP2D6_CNVx9 <- NULL  # Remove the existing column
              }
              if ("CYP2D6_CNVInt6" %in% colnames(temp_df) && "CYP2D6_CNVInt6" %in% colnames(cnv_data)) {
                temp_df$CYP2D6_CNVInt6 <- NULL  # Remove the existing column
              }
              if ("CYP2D6_CNV5UTR" %in% colnames(temp_df) && "CYP2D6_CNV5UTR" %in% colnames(cnv_data)) {
                temp_df$CYP2D6_CNV5UTR <- NULL  # Remove the existing column
              }
              # Add to temp_df
              temp_df <- cbind(temp_df, cnv_data)
            } else {
              warning("No matching rows for the patient in the CNV data.")
            }
          }
        }

        #call Find diplotype function
        find_diplotype_result <- find_diplotype(temp_df, gene, CYP1A2_name)
        
        #store results
        result_df <- find_diplotype_result[[1]]
        missing_geno <- find_diplotype_result[[2]]
        total_score <- find_diplotype_result[[3]]
        result_status <- result_df[[diplotype_col]]
        if (any(result_status != "No matching star alleles found")) {
          match_results <- rbind(match_results, data.frame(Diplotype = result_status, Score = total_score))
        }
      }

      # Select best match
      if (nrow(match_results) > 0) {
        max_score <- max(match_results$Score)
        best_matches <- match_results[match_results$Score == max_score, "Diplotype"]
        unique_matches <- unique(sapply(best_matches, function(dip) {
          components <- strsplit(dip, "/")[[1]]
          paste(sort(components), collapse = "/")
        }))
        df[row_idx, diplotype_col] <- unique_matches[1]
        if (length(unique_matches) > 1) {
          df[[alternate_col]][[row_idx]] <- unique_matches[-1]
        }
      } else {
        df[row_idx, diplotype_col] <- NA
        df[[alternate_col]][[row_idx]] <- NA
      }

      # Store the updated row's result in the list
      all_gene_results[[row_idx]] <- df[row_idx, ]
    }
  }

  # Combine the list into a single dataframe
  final_results <- do.call(rbind, all_gene_results)

  # Remove alternate diplotype columns where all values are NA or NULL
  for (gene in genes) {
    alternate_col <- paste0(gene, "_alternate_diplotype")
    if (all(sapply(final_results[[alternate_col]], function(x) is.null(x) || is.na(x)))) {
      final_results[[alternate_col]] <- NULL
    }
  }

  # Convert list columns to character strings before returning
  final_results <- convert_list_to_string(final_results)
  for(gene in genes){
    if(gene == "CYP2D6"){
      #adjust the potential alternate diplotypes 
      if("CYP2D6_alternate_diplotype" %in% colnames(final_results) || 
         ("CYP2D6_diplotype" %in% colnames(final_results) && 
          any(final_results$CYP2D6_diplotype %in% c("*10/*4", "*4/*10")))){
        final_results <- adjust_diplotype(final_results)
   
        #check again if alternate diplotype column is all NA, if yes then remove it
        if (all(sapply(final_results[["CYP2D6_alternate_diplotype"]], function(x) is.null(x) || is.na(x)))) {
          final_results[["CYP2D6_alternate_diplotype"]] <- NULL
        }
      }
    }
  }

  #reorder diplotypes so smallest number is first
  # Identify diplotype columns (columns containing "diplotype" in the name)
  diplotype_cols <- grep("diplotype", colnames(final_results), value = TRUE)
 
  # Apply the reordering function to all diplotype columns
  #final_results[diplotype_cols] <- lapply(final_results[diplotype_cols], reorder_diplotype)
  final_results[diplotype_cols] <- lapply(final_results[diplotype_cols], function(x) {
    # Convert to character if the column is factor
    x <- as.character(x)
    # Apply the reorder_diplotype function to each element
    sapply(x, reorder_diplotype)
  })

  return(list(final_results, missing_geno))
}


#' @title Helper function for assign_diplotype()
#' @description Convert any list columns to characters in order to write results to csv
#' @param df The dataframe to check
#' @return df The dataframe with no list columns
#' @export
# Function to convert list columns to character
convert_list_to_string <- function(df) {
  for (col in names(df)) {
    if (is.list(df[[col]])) {
      # Convert list to a character string, separating elements by a delimiter (e.g., ";")
      df[[col]] <- sapply(df[[col]], function(x) {
        if (is.null(x)) return(NA)  # Handle NULL values in the list
        paste(x, collapse = ";")
      })
    }
  }
  return(df)
}
        
#' @title Complete Summary Table Function 
#' @description Combine the dataframes from assign_diplotype, star_to_activity, and star_to_pheno into one large dataframe containing all the information in one place 
#' @param df The dataframe that contains the genotypes you want to convert
#' @param genes The specific CYP genes, as a list of strings
#' @return 1. The associated star alleles, metabolizer phenotypes, and activity scores for each genotype, formatted as a dataframe. 2. A message detailing missing genotypes. 
#' @examples 
#' summary1 <- all_geno_pheno(df, c("CYP2C9", "CYP2B6"))
#' summary2 <- all_geno_pheno(df2, "CYP2C19")
#' @export
all_geno_pheno <- function(df, genes, phased = FALSE, CYP1A2_name = "new") {
  genes <- toupper(genes)
  # Add star alleles to dataframe
  temp <- assign_diplotype(df, genes, phased, CYP1A2_name)
  results <- temp[[1]]
  missing_geno <- temp[[2]]

  # Identify the patient ID column dynamically
  patient_id_column <- names(df)[grep("ID", names(df), ignore.case = TRUE)][1]
  
  # Check if the patient ID column was found
  if (is.na(patient_id_column)) {
    stop("No patient ID column found. Ensure the dataframe contains a column with 'ID' in its name.")
  }
  
  for (gene in genes) {
    # Get the activity scores for the star alleles
    activity_df <- star_to_activity(results, gene)
    
    # Get the metabolizer phenotypes for the star alleles
    pheno_df <- star_to_pheno(activity_df, gene)
  
    # Append the newly calculated activity score and phenotype columns
    results[[paste0(gene, "_AS_1")]] <- pheno_df[[paste0(gene, "_AS_1")]]
    results[[paste0(gene, "_AS_2")]] <- pheno_df[[paste0(gene, "_AS_2")]]
    results[[paste0(gene, "_Phenotype")]] <- pheno_df[[paste0(gene, "_Phenotype")]]
    results[[paste0(gene, "_Phenotype_2")]] <- pheno_df[[paste0(gene, "_Phenotype_2")]]
    results[[paste0(gene, "_Comment")]] <- pheno_df$Comment
  }

  #Reorder dataframe columns so gene diplotypes, AS, phenotypes and comments are grouped together
  # Define the key categories we want to reorder
  categories <- c("diplotype", "alternate_diplotype", "AS_1", "AS_2", "Phenotype", "Phenotype_2", "Comment")
  
  # Identify the gene-specific columns based on pattern matching for each category
  gene_specific_columns <- grep(paste(categories, collapse = "|"), names(results), value = TRUE, ignore.case = TRUE)
  
  # Filter out columns that match the pattern but also contain a gene prefix, ensuring we only get relevant columns
  selected_gene_columns <- gene_specific_columns[grep("_(diplotype|alternate_diplotype|AS_1|AS_2|Phenotype|Phenotype_2|Comment)$", gene_specific_columns)]
  
  # Identify other columns to retain in their original order
  other_columns <- setdiff(names(results), selected_gene_columns)
  
  # Extract unique genes from selected columns by splitting on the underscore
  genes <- unique(sub("_(diplotype|alternate_diplotype|AS_1|AS_2|Phenotype|Phenotype_2|Comment)$", "", selected_gene_columns))
  
  # Initialize an empty vector to store the ordered gene-specific columns
  ordered_gene_columns <- c()
  
  # For each gene, find the relevant columns in the order of categories and add to the ordered list
  for (gene in genes) {
    for (category in categories) {
      col_name <- paste(gene, category, sep = "_")
      if (col_name %in% selected_gene_columns) {
        ordered_gene_columns <- c(ordered_gene_columns, col_name)
      }
    }
  }
  
  # Reorder the dataframe by combining other columns with the ordered gene-specific columns
  results <- results[, c(other_columns, ordered_gene_columns)]
  
  return(list(results, missing_geno))
}


#' @title Pie Chart Phenotype Summary Function
#' @description Creates a pie chart to summarize the metabolizer phenotypes for each gene
#' @param df The dataframe that contains the metabolizer phenotypes
#' @param genes The specific CYP genes, as a list of strings
#' @return A summary pie chart for each input gene
#' @examples 
#' CYP2D6_pie_chart <- pie_chart(df, "CYP2D6")
#' @export
pie_chart <- function(df, genes) {
  genes <- toupper(genes)
  # Initialize empty list to store charts
  charts <- list()
  
  for (gene in genes) {
    #get phenotype column
    phenotype_col <- paste0(gene, "_Phenotype")
    # Check if the phenotype column exists
    if (!phenotype_col %in% names(df)) {
      warning(paste("Phenotype column for", gene, "not found. Skipping."))
      next
    }
      
    # Count occurrences of each metabolizer phenotype
    slices <- table(df[[phenotype_col]])
    
    # Labels for the phenotypes
    labels <- names(slices)
    
    # Calculate percentages
    percentages <- round(100 * slices / sum(slices), 1)
    
    # Create labels with percentages
    labels_with_percentages <- paste(labels, percentages, "%", sep = " ")
    
    # Create pie chart
    summary_chart <- pie(slices, labels = labels_with_percentages, 
                         main = paste("Summary of", gene, "Metabolizer Phenotypes"), 
                         col = cm.colors(length(slices)))
    
    # Append the chart to the list
    charts[[gene]] <- summary_chart
  }
  
  return(charts)
}
