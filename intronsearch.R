#This script was used to extract information regarding to the number of introns contained within a certain gene for all of the samples in SRA. This specific script was written for BRCA1.
#However, this method was used for every gene searched. 


#to begin load in the dataframe containing all of the introns for the SRA samples
all_SRA_introns <- read.table(file="/dcl01/leek/data/sraintrons/all_SRA_introns.tsv.gz",sep="\t", header = FALSE) 

#change the names of the columns, so that they are easier to search with
names(all_SRA_introns)[1] <- "chr"
names(all_SRA_introns)[2] <- "start" 
names(all_SRA_introns)[3] <- "end"
names(all_SRA_introns)[4] <- "strand" 
names(all_SRA_introns)[5] <- "left_motif" 
names(all_SRA_introns)[6] <- "right_motif" 
names(all_SRA_introns)[7] <- "index" 
names(all_SRA_introns)[8] <- "reads" 

#To extract the info on introns within a given gene region, we obtained chr,start,and end information from biogps.org. 
#Using subset() we were able to find all rows within all_SRA_introns that contained introns within our given region. 
BRCA1Introns <- subset(all_SRA_introns,subset=(chr=="chr17" & start>=43044295 & end<=43125483)) 

#The dataframe was then saved for use in another script. 
save(BRCA1Introns,file="BRCA1Introns.Rda") 


