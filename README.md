# GenoStaR
An R package to convert genotypes to diplotypes and provide activity scores and metabolizer phenotypes.

# Package Overview
This package provides an efficient process to call diplotypes from genotypes based on allele definition tables provided by PharmGKB. You provide the genotypes, from .csv format and read in as a dataframe, and genes of interest and GenoStaR will call diplotypes, activity scores and metabolizer phenotypes. In the case where an alternate diplotype is returned, activity scores and phenotypes are based off the main diplotype call, the information is not applicable to the alternate diplotype. 

# Installation Instructions
#Ensure that devtools is installed:

install.packages("devtools")

#Install GenoStaR from GitHub:

devtools:: install_github(“GenoStaR-Genomics-Tools /GenoStaR”)

# Usage 
The package is easy to use. Most functions accept genotypes as a dataframe, with SNPs as column names and a list of genes for diplotype calling. The main function, all_geno_pheno(), generates diplotypes, metabolizer phenotypes, and activity scores. You can use functions individually, such as assign_diplotypes() for diplotypes only or star_to_pheno() for metabolizer phenotypes if diplotypes are already available. For CYP2D6, if copy number variant information is provided in the data, the package expects the column names to contain CNVx9 and CNVInt6. The phased parameter is used to mark if the input genotypes have been phased already (set to TRUE) or if they are in an unphased format (set to FALSE). Further function specific documentation is provided by running ?function_name(). 

It’s important to note that for CYP1A2, for any SNP that does not have an rsID, the SNP name will be referred to by its position, however the package expects special characters such as “-” and “>” to be removed. For example, the SNP -3954T>G should have the form X3954TG in the input data. Additionally, for CYP1A2, users have the option of specifying if they would like diplotypes to be returned following the old (pre-December 16th, 2024) or new (post Decemeber 16th, 2024) nomenclature as given by PharmVar. This can be done using the parameter ‘CYP1A2_name’ = “old” or “new”. The default is set to use the new nomenclature. 

#load the package

library(GenoStaR)

#sample call 

all_geno_pheno(genotypes, c(“CYP2D6”, “CYP2C9”), phased = FALSE, CYP1A2_name = “new”)

# Contributing
If any issues are found, please send an email to: alexcoulterr@gmail.com

# Authors
Alex Coulter, Arun Tiwari, Clement Zai

# References
M. Whirl-Carrillo1, R. Huddart1, L. Gong, K. Sangkuhl, C.F. Thorn, R. Whaley and T.E. Klein. "An evidence-based framework for evaluating pharmacogenomics knowledge for personalized medicine" Clinical Pharmacology & Therapeutics (2021) Sep;110(3):563-572. doi: 10.1002/cpt.2350. Epub 2021 Jul 22.


