#Specifying path used in the .sh file. This is the path to the 500 coverage bigwigs
#The path used in the .sh file = /dcl01/leek/data/sraonrail/sra_batch_0_sample_size_500_align/coverage_bigwigs
args <- commandArgs(trailingOnly = TRUE) 
path <- args[1] 

#loading in the readcount table. The file was copied from /dcl01/leek/data/sraonrail/sra_batch_0_sample_size_500_align/cross_sample_results/counts.tsv.gz
#This table contained information on the number of mapped reads for each chromsome, as well as the number of total mapped reads for each sample. 
counts <- read.table(file = '/home/other/rphillip/counts.tsv.gz',header = TRUE, sep = '\t', stringsAsFactors = FALSE)

#Creating a vector of the 500 bigwig file names. This will be inputted as the file parameter.  
f <-  scan("sra_samples.txt", what="", sep="\n")

#loading the libraries 
library(rtracklayer)
library(GenomicRanges)

normalized <- function(file,path,chr,start,end) {
    #function from Jean-Phillipe to extract the region of interest from the bigwig file. 
    extract.block <- function(files, chr, start, end, verbose = TRUE){
      rl <- IRanges::RangesList(IRanges::IRanges(start=start, end=end))
      names(rl) <- chr
      rles <- lapply(files, function(xx) {
        import(xx, as = "Rle", format = "bw", selection = BigWigSelection(rl))
      })
      megaMatrix <- do.call(cbind, lapply(rles, function(xx)      as.numeric(xx[[chr]][start:end])))
      megaMatrix
    }
  #creating an empty vector to input the normalized values into. length(file) is used here to ensure that the vector is the same length as the number of files inputted. 
  normalized_values <- rep(NA,length(file))
  if(length(file)>=1)
  {
    for(i in 1:length(file))
    {
      #the first column in counts.tsv.gz contains the name of the file. However, the bigwig files pulled from 
      #/dcl01/leek/data/sraonrail/sra_batch_0_sample_size_500_align/coverage_bigwigs have a .bw ending on them. 
      # To get the total number of mapped reads for each sample we have to use strsplit to get rid of the .bw and find the rile name
      #in the first column of counts. 
      c <- grep((strsplit(file[i],'.bw')[[1]][1]),counts$X)
      #The total.mapped.reads column in counts.tsv.gz is comma separated. The first value corresponds to the total mapped reads in the .bw files. The second value
      #is the total mapped reads for the unique.bw file. Therefore, I used strsplit and grapped the first value. 
      total <- as.numeric(strsplit(as.character(counts[c,'total.mapped.reads']),',')[[1]][1])
      print(paste0(path,file[i]))# Ensures that the for loop is running and that it is looking in the right place for the files. 
      #Extracting the block of interest using one of the files, and the chr,start, and end specified in the function.
      block <- extract.block(files = paste0(path, file[i]),
                             chr = chr,
                             start = start,
                             end = end,
                             verbose = TRUE)
      #Implementing the actual normalization factor. This divides each coverage value by the total number of mapped reads 
      #and mulitplies it by 40000000, and taking the mean
      x <- mean(block/(total))* 40000000
      #Plugging in the normalized values into the empty normalized_value vector created earlier. 
      normalized_values[i] <- x
    }
   # This script was used for running normalization locally. To change the name of the columns, and the actual saved dataframe, we just changed the name of the dataframe 
   # each time to the desired title.
   # Creating a data frame with the names of the files, and their corresponding normalized value
   USP9Y_SRA  <- data.frame(f,normalized_values)
   # Changing the name of the columns. 
   names(GDF3_SRA)[1] <- "file"
   names(GDF3_SRA)[2] <- "USP9Y"
   #Saving the dataframe. The dataframes were then imported into Rstudio using WinSCP amd the load() command. 
   save(USP9Y_SRA,file = "USP9Y_SRA.Rda")
  }
}


#implementing the normalized function. 
normalized(file = f, path = path, chr = "chrY", start = 14813160, end = 14972768)


