#This script was used to sum the number of reads mapping to introns within a given region of interest. The region of interest are
#often regions that contain genes. 


library(stringr) 

#Load in the previous dataframe, saved from the intronsearch.R 
load("~/BRCA1Introns.Rda")

#change the names of the columns
names(BRCA1Introns)[1] <- "chr" 
names(BRCA1Introns)[2] <- "start" 
names(BRCA1Introns)[3] <- "end" 
names(BRCA1Introns)[4] <- "strand" 
names(BRCA1Introns)[5] <- "left_motif" 
names(BRCA1Introns)[6] <- "right_motif" 
names(BRCA1Introns)[7] <- "index" 
names(BRCA1Introns)[8] <- "number_of_reads_mapping" 

#Load in the index file that matches the numbers within the seventh column of the introns dataframe, with their 
#run accessions. 
index <- read.table(file="/dcl01/leek/data/sraintrons/index_to_SRA_accession.tsv",sep="\t",header=FALSE) 

#The first column corresponds to the number within the introns dataframe that represents the run_accession 
names(index)[1] <- "Index" 
names(index)[2] <- "study_accession"
names(index)[3] <- "sample_accession"
names(index)[4] <- "experiment_accession" 
names(index)[5] <- "run_accession" 

#The seventh and eighth column are comma separated. We have to separte the numbers and create new vectors containing the values. 
numbers <- as.vector(BRCA1Introns[,7])
indexnumbers <- vector() 
for(i in 1:length(numbers)){
  everynumber <- str_split(numbers[i],",") 
  indexnumbers <- append(indexnumbers,everynumber) %>% unlist() %>% as.numeric() 
}

reads <- as.vector(BRCA1Introns[,8]) 
readnumbers <- vector()
for(i in 1:length(reads)){
  everyreadcount <- str_split(reads[i],",") 
  readnumbers <- append(readnumbers,everyreadcount) %>% unlist() %>% as.numeric() 
} 

# Creating a dataframe with the indexnumbers and readnumbes vectors that were pulled out and unlisted from the original dataframe
indexandreads  <- data.frame(indexnumbers,readnumbers) 


#indexandreads contains the index number and the number of reads 
#This for loop gets the run_accessions for every index. 
runs <- vector() 
for(i in 1:nrow(indexandreads)){
  row <- grep(indexandreads[i,1],index[,1])[1] 
  run_accessions <- as.character(index[row,"run_accession"])
  runs <- append(runs,run_accessions) 
}

#Ensuring the loop worked
runs[1:10] 

#Adding the run_accessions to the indexandreads dataframe. 
indexandreads <- cbind(indexandreads,runs) 

#Here we split the data by the run_accessions (runs). 
#Splitting the data by run_accessions was necessary because indexandreads contained multiple rows that have the same run_accessopm 
#By splitting the data we get a list, where each component of the list is a unique run_accession, and their corresponding rows 
#consitsting of the number of reads mapping to each intron. 
splitdata <- split(indexandreads,indexandreads$runs) 
#Getting the run_accession names in the order they appear in on the list
splitdatanames <- names(splitdata) 

#Summing all of the inrons for each run_accession  
splitdatasums <- vector() 
for(i in 1:length(splitdata)){
 sums <- sum(splitdata[[i]]["readnumbers"]) 
 splitdatasums <- append(splitdatasums,sums) 
}

#Making the final dataframe that conists of the run_accessions and the summed number. 
BRCA1intronsums <- data.frame(splitdatanames,splitdatasums) 

#saving the dataframe. 
save(BRCA1intronsums,file="BRCA1intronsums.Rda") 

