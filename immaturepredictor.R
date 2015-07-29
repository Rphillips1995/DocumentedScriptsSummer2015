library(stringr)


#These dataframe contain normalized gene expression counts for different genes 
#These genes were chosen because they are master body plan genes, or are overexpressed in immature cells.  
#Load in the HOXA1 expression data 
load("~/STUFF/HOXA1_SRA.Rda")
#Load in the ZEB2 expression data 
load("~/STUFF/ZEB2_SRA.Rda")
#Load in the SOX2 expression data 
load("~/STUFF/SOX2_SRA.Rda")
#Load in the OCT4 expression data 
load("~/STUFF/OCT4_SRA.Rda")
#Load in the PARK7 expression data 
load("~/STUFF/PARK7_SRA.Rda")
#Load in the RAN expression data 
load("~/STUFF/RAN_SRA.Rda")
#Load in the NANOG expression data 
load("~/STUFF/NANOG_SRA.Rda")
#Load in the FOXD3 expression data 
load("~/STUFF/FOXD3_SRA.Rda")
#Load in the BRAF expression data 
load("~/STUFF/BRAF_SRA.Rda")


#Build the immaturegenes data frame 
immature <- merge(HOXA1_SRA,ZEB2_SRA, by = "file")
immature <- merge(immature,SOX2_SRA,by = "file")
immature <- merge(immature, OCT4_SRA, by = "file")
immature <- merge(immature, PARK7_SRA, by  = "file")
immature <- merge(immature,RAN_SRA, by = "file")
immature <- merge(immature,NANOG_SRA, by = "file")
immature <- merge(immature,BRAF_SRA, by = "file")
immature <- merge(immature,FOXD3_SRA, by = "file")

#The names of the files are long strings of the accession numbers 
# We want to pull out the run accession,so that we can merge the metadata and the expresssion data by
#run_accession 
names <- rep(NA,500)
files <- as.vector(immature[,1])
for(i in 1:length(files)){
  x <- str_split(files[i],"_|-")[[1]][4]
  names[i] <- x
}
immature[,1]<- names
names(immature)[1] <- "run_accession"

#merge updatemetada and immature. 
#We only have gene expression values for 500 samples because we had 500 bigwig files for which we could pull gene expression values from. 
update_metadata_SRA500 <- merge(update_metadata,immature,all=FALSE)
#give update metadata an ID 
update_metadata_SRA500$ID <- 1:nrow(update_metadata_SRA500)

#Here we grep the update_metadata_SRA500 for rows that contain words that contain some sort of immature keyword. 
immaturesource <- grep("fibroblast|myoblast|neonatal|embryonic|HEK293|embryo|gestational|erythroblast|H1-HESC|TRIMESTER", update_metadata_SRA500$source_name, ignore.case = TRUE)
immaturecelltype <- grep("fibroblast|myoblast|neonatal|embryonic|HEK293|embryo|gestational|erythroblast|H1-HESC|TRIMESTER", update_metadata_SRA500$cell_type, ignore.case = TRUE)
immatureage <- grep("fibroblast|myoblast|neonatal|embryonic|HEK293|embryo|gestational|erythroblast|H1-HESC|TRIMESTER", update_metadata_SRA500$age, ignore.case = TRUE)
immatureorigin <- grep("fibroblast|myoblast|neonatal|embryonic|HEK293|embryo|gestational|erythroblast|H1-HESC|TRIMESTER", update_metadata_SRA500$sample_origin, ignore.case = TRUE)
immautreline <-grep("fibroblast|myoblast|neonatal|embryonic|HEK293|embryo|gestational|erythroblast|H1-HESC|TRIMESTER", update_metadata_SRA500$cell_line, ignore.case = TRUE)
immatureattribute <- grep("fibroblast|myoblast|neonatal|embryonic|HEK293|embryo|gestational|erythroblast|H1-HESC|TRIMESTER", update_metadata_SRA500$sample_attribute, ignore.case = TRUE)

#Putting all of the outputs from the grep searches above into one vecotr
immaturerows <- c(immaturesource,immaturecelltype,immatureage,immatureorigin,immatureline,immatureattribute)
#some rows are repeated so we only look for the unique row numbers
immaturerows <- unique(immaturerows)

#Creating the binary immature column within the update_metadata_SRA500 dataframe, so that we can apply the predictor. 
update_metadata_SRA500$Immature <- ifelse ((update_metadata_SRA500$ID %in% immaturerows)==TRUE, 
                                           'YES',
                                           "NO")


library(caret)
#creating the training  and testing datadata
inTrain <- createDataPartition(y=update_metadata_SRA500$Immature,p=.5,list=FALSE)
training <-update_metadata_SRA500[inTrain,]
testing <- update_metadata_SRA500[-inTrain,]
training$Immature<-as.factor(training$Immature)
testing$Immature <- as.factor(testing$Immature)

#Training the predictor
ctrl <- trainControl(method = "repeatedcv",
                     repeats = 3,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)
plsFit <- train(Immature ~ PARK7+FOXD3+HOXA1+ZEB2+NANOG+SOX2+GDF3,
                data = training,
                trControl = ctrl,
                method = "pls",
                metric = "ROC",
                preProc = c("center","scale"))

#Implementing the predictor on the testing data. 
plsProbs <- predict(plsFit, newdata = testing)
confusionMatrix(data=plsProbs, testing$Immature)


