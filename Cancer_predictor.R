# This script is extremely similar to the immature predictor script. 
# This script was used to try and determine whether or not a sample was cancerous based on gene expression values. 
# Reference immaturepredictor.R 

library(stringr)


#Normalized Gene Expression data for TP53 
load("~/STUFF/TP53_SRA.Rda")
#Normalized Gene Expression data for AR 
load("~/STUFF/AR_SRA.Rda")
#Normalized Gene Expression data for ACTB 
load("~/STUFF/ACTB_SRA.Rda")
#Normalized Gene Expression data for RASA1 
load("~/STUFF/RASA1_SRA.Rda")

#Build the dataframe with all of the normalized values 
cancergenes <- merge(AR_SRA,TP53_SRA,by = "file")
cancergenes <- merge(cancergenes,ACTB_SRA,by="file")
cancergenes <- merge(cancergenes,RASA1_SRA,by="file")

#The names of the files are long strings of the accession numbers 
# We want to pull out the run accession,so that we can merge the metadata and the expresssion data by
#run_accession 
names <- rep(NA,500)
files <- as.vector(cancergenes[,1])
for(i in 1:length(files)){
  x <- str_split(files[i],"_|-")[[1]][4]
  names[i] <- x
}

#change the first column to run_accession 
cancergenes[,1] <- names
names(cancergenes)[1] <- "run_accession"

#load in updated metadata
load("~/STUFF/update_metadata.Rda")


#Merge cancergenes and metadata
update_metadata_SRA500 <- merge(update_metadata,cancergenes,all=FALSE)
#ID 
update_metadata_SRA500$ID <- 1:nrow(update_metadata_SRA500)

#Get all of the rows that contain some kind of cancer
cancerdisease <- grep("CARCINOMA|MELANOMA|NEUROBLASTOMA|OSTEOBLASTOMA|OSTEOSARCOMA|LEUKEMIA|TUMOUR|BREAST CANCER|OVARIAN CANCER|ADENOCARCINOMA|ASTROCYTOMA|SARCOMA",update_metadata_SRA500$disease)
cancerabstract <- grep("CARCINOMA|MELANOMA|NEUROBLASTOMA|OSTEOBLASTOMA|OSTEOSARCOMA|LEUKEMIA|TUMOUR|BREAST CANCER|OVARIAN CANCER|ADENOCARCINOMA|ASTROCYTOMA|SARCOMA|CANCER|GLIOBLASTOMA",update_metadata_SRA500$study_abstract,ignore.case = TRUE)
canceratt <- grep("CARCINOMA|MELANOMA|NEUROBLASTOMA|OSTEOBLASTOMA|OSTEOSARCOMA|LEUKEMIA|TUMOUR|BREAST CANCER|OVARIAN CANCER|ADENOCARCINOMA|ASTROCYTOMA|SARCOMA|CANCER|GLIOBLASTOMA",update_metadata_SRA500$sample_attribute,ignore.case = TRUE)

#Get the rows with cancer
cancers <- c(cancerdisease,cancerabstract,canceratt)
cancer <- unique(cancers)

#Give update_metadata_SRA500 an id so we can match the cancers vector to the row numbers 
update_metadata_SRA500$Cancer <- ifelse ((update_metadata_SRA500$ID %in% cancer)==TRUE,
                                         'YES',
                                         "NO")


library(caret)

#Create the training and testing data
inTrain <- createDataPartition(y=update_metadata_SRA500$Cancer,p=.5,list=FALSE)
training <-update_metadata_SRA500[inTrain,]
testing <- update_metadata_SRA500[-inTrain,]
training$Cancer<-as.factor(training$Cancer)
testing$Cancer <- as.factor(testing$Cancer)

#Train the data
ctrl <- trainControl(method = "repeatedcv",
                     repeats = 3,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)
plsFit <- train(Cancer ~ AR+TP53+RASA1+ACTB,
                data = training,
                trControl = ctrl,
                method = "pls",
                metric = "ROC",
                preProc = c("center","scale"))

#Implement the model on the testing data. 
plsProbs <- predict(plsFit, newdata = testing)
confusionMatrix(data=plsProbs, testing$Cancer)
