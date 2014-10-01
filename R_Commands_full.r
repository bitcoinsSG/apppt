memory.limit(4095)
source("Load_Arsenal_APPPt_3_2_pept_mass.r")                                                        ##Load Fuctions
cols<-rep((brewer.pal(10,"Set3")),20)                                                               ##Cols for graphing

#####Load Training Data Set###########3

read.table("./Raw_Expression/Expression_Training",header=TRUE,row.names=1,na.string=0)->Ex.train                ##Load Expression Training
read.csv("./Raw_Expression/12672931360570Phenotype_Training.csv",header=TRUE)->Ph.train             ##Load Phenotype Training

names(Ex.train)<-gsub(pattern="[A-z]",replacement="",x=names(Ex.train),perl=TRUE)  ## Remove Xs from col.names in dataframe Ex.train

row.names(Ph.train)<-Ph.train[,2]                                                             ## rowname the phenotype dataframe with indentifier 
Ph.train<-Ph.train[order(as.numeric(row.names(Ph.train))),]                                              ## Sort the rows such that the identifier is in ascending order
Ex.train<-Ex.train[,order(as.numeric(names(Ex.train)))]                                       ## Sort the columns in the expression dataset to match order of phenotype dataframe identifier

Ex.train.na.gene<-as.data.frame(apply(apply(Ex.train,1,is.na),2,sum))              ## Calculate how many na's there are for each gene
colnames(Ex.train.na.gene)<-"NA_Count_Gene_Train"                                        ## Rename resultant data frame to reflect content
Ex.train.na.sample<-as.data.frame(t(as.data.frame(apply(apply(Ex.train,1,is.na),1,sum))))    ## Calculate how many na's there are for each sample
row.names(Ex.train.na.sample)<-"NA_Count_Sample_Train"                                             ## Rename Resultant data fram to reflect content

boxplot(log(Ex.train),las=2,col=cols)  #graph 

###########################################


####################Load Test######################


read.table("./Raw_Expression/Expression_Test_Corrected.txt",header=TRUE,row.names=1,na.string=0)->Ex.test                ##Load Expression Test
read.csv("./Raw_Expression/12672927740810Phenotype_Test.csv",header=TRUE)->Ph.test             ##Load Phenotype Test

names(Ex.test)<-gsub(pattern="[A-z]",replacement="",x=names(Ex.test),perl=TRUE)  ## Remove Xs from col.names in dataframe Ex.test

row.names(Ph.test)<-Ph.test[,2]                                                             ## rowname the phenotype dataframe with indentifier 
Ph.test<-Ph.test[order(as.numeric(row.names(Ph.test))),]                                              ## Sort the rows such that the identifier is in ascending order
Ex.test<-Ex.test[,order(as.numeric(names(Ex.test)))]                                       ## Sort the columns in the expression dataset to match order of phenotype dataframe identifier

Ex.test.na.gene<-as.data.frame(apply(apply(Ex.test,1,is.na),2,sum))              ## Calculate how many na's there are for each gene
colnames(Ex.test.na.gene)<-"NA_Count_Gene_Test"                                        ## Rename resultant data frame to reflect content
Ex.test.na.sample<-as.data.frame(t(as.data.frame(apply(apply(Ex.test,1,is.na),1,sum))))    ## Calculate how many na's there are for each sample
row.names(Ex.test.na.sample)<-"NA_Count_Sample_Test"                                             ## Rename Resultant data fram to reflect content

boxplot(log(Ex.test),las=2,col=cols)  #graph 

#######################################################

#####Example of sorting test dataset by col and then outputting in format ready for submission
toString(Ph.test[                   order(Ex.test.na.sample[1,])                               ,1])       ## order sample by how many missing values it has
Ex.train[(!is.na(match(row.names(Ex.train),c("Gene31193","Gene2536","Gene8854","Gene21492")))),]->pl ##Extracting only those rows that match the rownames in c(""""")