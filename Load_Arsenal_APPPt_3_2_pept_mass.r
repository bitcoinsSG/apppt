###TIP: Creating a dataframe with emperically defined row and col names look at the Normalization section

 # The command looks into the below mentioned directory and expects only three files namely...
#gct_file<-paste("./Parsed_Files_For_Bioconductor/",list.files("./Parsed_Files_For_Bioconductor/")[1],sep="")
#gmt_file<-paste("./Parsed_Files_For_Bioconductor/",list.files("./Parsed_Files_For_Bioconductor/")[2],sep="")
#redundant_proteins_file_name<-paste("./Parsed_Files_For_Bioconductor/",list.files("./Parsed_Files_For_Bioconductor/")[3],sep="")

                                #Library Loading Time    
##########################################################################################

#You may need to be root to excute access command, but these packages are installed so you should be fine.
#source("http://bioconductor.org/biocLite.R")
#biocLite("safe")
#biocLite("GSA")
#Command to install a package from command line
#R CMD INSTALL --configure-args='--enable-old-style-formdata' CGIwithR_0.72.tar.gz


## Must have package because of plots : http://bioconductor.org/packages/2.4/bioc/vignettes/HELP/inst/doc/HELP.pdf
## FDR for resampling :http://bioconductor.org/packages/2.4/bioc/html/fdrame.html
## Protein Annotation builder cool stuff :http://bioconductor.org/packages/2.4/bioc/vignettes/PAnnBuilder/inst/doc/PAnnBuilder.pdf
## pca MISSING VALUE IMPUTATION ETC http://bioconductor.org/packages/2.4/bioc/html/pcaMethods.html
### plgem Washburn http://bioconductor.org/packages/2.4/bioc/html/plgem.html
## John Storey Q Value http://bioconductor.org/packages/2.4/bioc/html/qvalue.html
## pre-processing functions http://bioconductor.org/packages/2.4/bioc/manuals/preprocessCore/man/preprocessCore.pdf
## Parametric GSEA :http://bioconductor.org/packages/2.4/bioc/html/PGSEA.html
## Non parametric  Gene set type : http://bioconductor.org/packages/2.4/bioc/vignettes/RankProd/inst/doc/RankProd.pdf
## Package ROC : http://bioconductor.org/packages/2.4/bioc/html/ROC.html
## Details abuout safe data types :http://bioconductor.org/packages/2.4/bioc/manuals/safe/man/safe.pdf
## Structured analysis of microarray, vrey extensive and useful : http://bioconductor.org/packages/2.4/bioc/manuals/stam/man/stam.pdf
## Visululization for gene set stuff : http://bioconductor.org/packages/2.4/bioc/vignettes/topGO/inst/doc/topGO.pdf
## Another gene set type graphing tool : http://bioconductor.org/packages/2.4/bioc/vignettes/SPIA/inst/doc/SPIA.pdf
## Cool heatmaps :http://bioconductor.org/packages/2.4/bioc/vignettes/Heatplus/inst/doc/Heatplus.pdf
##GSEA :http://bioconductor.org/packages/2.4/bioc/vignettes/GSEABase/inst/doc/GSEABase.pdf
## Linear models in GSEA http://bioconductor.org/packages/2.4/bioc/vignettes/GSEAlm/inst/doc/GSEAlm.pdf

## Really good graphing : http://bioconductor.org/packages/2.4/bioc/vignettes/affycomp/inst/doc/affycomp.pdf
### Tutorial on plots etc :http://lmdvr.r-forge.r-project.org/figures/figures.html
## Sophisticated Imputation : http://cran.r-project.org/web/packages/Amelia/vignettes/amelia.pdf
## Boot Strapping : http://cran.r-project.org/web/packages/boot/index.html
## Graphing for multivariate :http://cran.r-project.org/web/packages/denpro/index.html
## 3d plots http://cran.r-project.org/web/packages/vrmlgen/
## Opengl : http://cran.r-project.org/web/packages/rgl/rgl.pdf
##  R graphics main link : http://addictedtor.free.fr/graphiques/
library("safe")
library("MantelCorr")
library("minet")
library("infotheo")
library("GSA") #Tibshirani's Gene Set Enrichment
library("multtest")
library("HELP")
library("sigPathway")
library("lattice")
library("affy")
#library("Rcmdr")
library("limma")
library("marray")
library("RColorBrewer")
library("impute")
library("pcaMethods")
library("samr")
library("gplots")
library("vsn")
library("PAnnBuilder")
library("Heatplus")
library("hopach")
#library("pcot2")           # This packages masks pcaMethods PCA function
library("hu6800.db")
library("plgem")
##PAnnBuilder
if (.Platform$OS.type == "windows") {
library("org.Hs.ipi")
 } else {
 library("org.Hs.ipi.db")
 }




###Function to install all packages above if package is not already installed
install_all_packages_required<-function(){
source("http://www.bioconductor.org/biocLite.R")
biocLite("safe")
biocLite("GSA") #Tibshirani's Gene Set Enrichment
biocLite("multtest")
biocLite("HELP")
biocLite("sigPathway")
biocLite("lattice")
biocLite("affy")
biocLite("Rcmdr")
biocLite("limma")
biocLite("marray")
biocLite("RColorBrewer")
biocLite("impute")
biocLite("pcaMethods")
biocLite("samr")
biocLite("vsn")
biocLite("gplots")
biocLite("PAnnBuilder")
biocLite("Heatplus")
biocLite("hopach")
biocLite("pcot2")
biocLite("hu6800.db")
biocLite("plgem")
biocLite("GO.db")
biocLite("MantelCorr")
biocLite("minet")
biocLite("infotheo")



if (.Platform$OS.type == "windows") {
tmpFile <- paste("http://www.biosino.org/download/PAnnBuilder",
"org.Hs.ipi_1.0.0.zip", sep = "/")
 tmpPath <- file.path(tempdir(), "org.Hs.ipi_1.0.0.zip")
 } else {
 tmpFile <- paste("http://www.biosino.org/download/PAnnBuilder",
 "org.Hs.ipi_1.1.0.tar.gz", sep = "/")
 tmpPath <- file.path(tempdir(), "org.Hs.ipi_1.1.0.tar.gz")
 }
download.file(tmpFile, tmpPath)
install.packages(tmpPath, repos = NULL)
library(org.Hs.ipi)
}





##########################################################################################
                				#Functions for Pre-Processing 
              				  # Type 1 Filtering
              				  # Type 2 Transforming
              				  # Type 3 Normalization


# Type 1 Filtering
##########################################################################################				
##########################################################################################
##########################################################################################
# Type 2 Transforming
##########################################################################################	
Pre_Processing_S_2_All_Transformations.list<-function(Filtered_gct_file.list){
## Accepts Data frame of peptide level data with NAs 
## Transforms and outputs list non unlike Pre_Processeing_S_3_Normalization function
## Outputs List with datframes with various Transformations
Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter<-Filtered_gct_file.list[[1]]
## To Make sure the data is inputted in data.frame format
Log_2_Transformed<-log(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter,2)
Unlogged_raw<-Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter
All_Transformed_Forms_of_Input_Data_Frame_as_list<-list(Unlogged_raw,Log_2_Transformed)
names(All_Transformed_Forms_of_Input_Data_Frame_as_list)<-c("Raw","Log(2)")
return(All_Transformed_Forms_of_Input_Data_Frame_as_list)
}
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

##########################################################################################
# Type 2.5 Imputation
##########################################################################################	
Pre_Processing_S_2.5_All_Imputations.list<-function(Filtered_Tranformed_gct_file.list.with.dataframe,nPcs=5,colmax=0.5,row.max=0.5,k=10,rng.seed=345){
## Accepts list with only one data frame in it, i.e. output from Pre_Processing_S_2_All_Transformations[index], of peptide level data with NAs
## Imputes and outputs list non unlike Pre_Processeing_S_3_Normalization function
## Outputs List with datframes with various imputation 
Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter<-Filtered_Tranformed_gct_file.list.with.dataframe[[1]]

## To Make sure the data is inputted in data.frame format
None_Imputed<-Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter
None_Imputed[is.na(None_Imputed)]<-0    #### Replaces NAs with 0

Half_Min_Imputed<-Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter  
## Replacing all NA with dataset min/2
##Half_Min_Imputed[is.na(Half_Min_Imputed)]<- (min(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter,na.rm=TRUE)/2) # #### Replaces NAs with min/2

## Replacing NA sample wise with sample min/2
sample_wise_min_divided_by_2<-apply(Half_Min_Imputed,2,min,na.rm=TRUE)/2
for(i in 1:dim(Half_Min_Imputed)[2]){
Half_Min_Imputed[is.na(Half_Min_Imputed[,i]),i]<-as.double(sample_wise_min_divided_by_2[i])
}
#####
BPCA_Imputed<-as.data.frame(pca(as.matrix(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter), nPcs=nPcs, method="bpca")@completeObs)
KNN_Imputed<-as.data.frame(impute.knn(as.matrix(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter),k, row.max, colmax, maxp = 1500, rng.seed=rng.seed)$data)
All_Imputed_Forms_of_Input_Data_Frame_as_list<-list(None_Imputed,Half_Min_Imputed,BPCA_Imputed,KNN_Imputed)
names(All_Imputed_Forms_of_Input_Data_Frame_as_list)<-c(  
 (paste(names(Filtered_Tranformed_gct_file.list.with.dataframe)," + None",sep="")),
 (paste(names(Filtered_Tranformed_gct_file.list.with.dataframe)," + MIN/2",sep="")),
 (paste(names(Filtered_Tranformed_gct_file.list.with.dataframe)," + BPCA",sep="")),
 (paste(names(Filtered_Tranformed_gct_file.list.with.dataframe)," + KNN",sep=""))   
 
  )
return(All_Imputed_Forms_of_Input_Data_Frame_as_list)
}
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

##########################################################################################
##########################################################################################
# Type 3 Normalization
# This function is a one stop shop for all normalization of columns in a data frame
# !!!!!!!!!Has not been testing for NA or missing, so Imputation should be  handled before
# In this context it will probably be used primarily for peptide level information
# Accepts Type Data Frame, returns a List of Data Frames with all the Normalized values
# Uses Limma and affy
# By Default it uses 100 as number of iterations, this is bery arbitary just to make sure are using what I think is a large enough n u m b e r
##### Also Assumes the title of the inpute dataframe is "Log(2) + Imputed" if not specified
Pre_Processing_S_3_All_Normalizations.list<-function(Input_List_Peptide_Level_But_Does_Not_Really_Matter,Number_of_Iterations_For_Loess=100,Type_of_Input_Data_Frame=names(Input_List_Peptide_Level_But_Does_Not_Really_Matter)){
print(Type_of_Input_Data_Frame)
## From Package Limma
Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter<-Input_List_Peptide_Level_But_Does_Not_Really_Matter[[1]]
###Make sure the data is inputted in data.frame format
Quantile<-as.data.frame(normalizeBetweenArrays(as.matrix(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter),method="quantile"))
Global_Scaled<-as.data.frame(normalizeBetweenArrays(as.matrix(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter),method="scale"))
##Since Variance Stabalizing Normalization already logs the data, the data is unlogged(2) before normalization if indicated as logged(default), or else is left raw as in its logged form
# Length of output == 1 if true and 0 if False
if(  length(grep("Raw",Type_of_Input_Data_Frame))==1   ){
##Conduct VSN on raw unlogged vaules
Variance_Stabalizing<-as.data.frame(normalizeBetweenArrays(  (as.matrix(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter)),method="vsn"))
## Convert the VSN's logged output back to unlogged
Variance_Stabalizing<-as.data.frame(2^(Variance_Stabalizing))
}
if(  length(grep("Log",Type_of_Input_Data_Frame))==1   ){
##Conduct VSN on unlogged vaules by first unlogging the values if name has Log(2) in it i.e. if you write another tranformation, that should be accounted for in another if statement
Variance_Stabalizing<-as.data.frame(normalizeBetweenArrays(    2^(as.matrix(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter))   ,method="vsn"))
## Since VSN's output is already logged do nothing since values are back to log2 form as inputted
}


## From Package Affy
Cubic_Spline<-as.data.frame(matrix(          normalize.qspline(as.matrix(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter))               ,ncol=dim(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter)[2],dimnames=dimnames(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter)))
Cyclic_Loess_Smoothing<-as.data.frame(matrix(normalize.loess( log.it=FALSE,verbose = TRUE, maxit=Number_of_Iterations_For_Loess ,mat=as.matrix(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter))               ,ncol=dim(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter)[2],dimnames=dimnames(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter)))
####max iteration of cyclic loess is defined as 100
Homemades_R_Squared<-as.data.frame(t(apply((apply(as.matrix(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter),2,rank)),1,rank)))
Homemades_R_Cubed<-as.data.frame(apply(apply((apply(as.matrix(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter),2,rank)),1,rank),1,rank))

All_Normalized_Forms_of_Input_Data_Frame_as_list<-list(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter,Quantile,Global_Scaled,Variance_Stabalizing,Cubic_Spline,Cyclic_Loess_Smoothing,Homemades_R_Squared,Homemades_R_Cubed)
names(All_Normalized_Forms_of_Input_Data_Frame_as_list)<-c(
(paste(Type_of_Input_Data_Frame," + None",sep="")),
(paste(Type_of_Input_Data_Frame," + Quantile",sep="")),
(paste(Type_of_Input_Data_Frame," + Global Scaled",sep="")),
(paste(Type_of_Input_Data_Frame," + Variance Stab",sep="")),
(paste(Type_of_Input_Data_Frame," + Cubic Sm Splines",sep="")),
(paste(Type_of_Input_Data_Frame,        (paste(" + Cyclic Loess Itera=",as.character(Number_of_Iterations_For_Loess),sep=""))  ,sep="")),
(paste(Type_of_Input_Data_Frame," + Spearman Homemade R^2",sep="")),
(paste(Type_of_Input_Data_Frame," + Spearman Homemade R^3",sep=""))


#Type_of_Input_Data_Frame," + Quantile"," + Global Scaled"," + Variance Stab"," + Cubic Sm Splines",(paste(" + Cyclic Loess Itera=",as.character(Number_of_Iterations_For_Loess),sep=""))," + Spearman Homemade R^2"," + Spearman Homemade R^3"
)
return(All_Normalized_Forms_of_Input_Data_Frame_as_list)

}                       				
##########################################################################################				
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
































##########################################################################################
####################################################################################################################################################################################
##########################################################################################
			
##########################################################################################
##########################################################################################
# Type 3 Normalization
# This function is a one stop shop for all normalization of columns in a data frame
# !!!!!!!!!Has not been testing for NA or missing, so Imputation should be  handled before
# In this context it will probably be used primarily for peptide level information
# Accepts Type Data Frame, returns a List of Data Frames with all the Normalized values
# Uses Limma and affy
# By Default it uses 100 as number of iterations, this is bery arbitary just to make sure are using what I think is a large enough n u m b e r
##### Also Assumes the title of the inpute dataframe is "Log(2) + Imputed" if not specified
oldPre_Processing_S_3_All_Normalizations<-function(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter,Number_of_Iterations_For_Loess=100,Type_of_Input_Data_Frame="Log(2) + Imputed"){

## From Package Limma
Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter<-as.data.frame(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter)
###Make sure the data is inputted in data.frame format
Quantile<-as.data.frame(normalizeBetweenArrays(as.matrix(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter),method="quantile"))
Global_Scaled<-as.data.frame(normalizeBetweenArrays(as.matrix(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter),method="scale"))
##Since Variance Stabalizing Normalization already logs the data, the data is unlogged(2) before normalization if indicated as logged(default), or else is left raw as in its logged form
if(Type_of_Input_Data_Frame=="Log(2) + Imputed"){
Variance_Stabalizing<-as.data.frame(normalizeBetweenArrays(2^(as.matrix(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter)),method="vsn"))
}
if(!(Type_of_Input_Data_Frame=="Log(2) + Imputed")){
Variance_Stabalizing<-as.data.frame(2^(normalizeBetweenArrays(as.matrix(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter),method="vsn")))

}


## From Package Affy
Cubic_Spline<-as.data.frame(matrix(          normalize.qspline(as.matrix(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter))               ,ncol=dim(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter)[2],dimnames=dimnames(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter)))
Cyclic_Loess_Smoothing<-as.data.frame(matrix(normalize.loess( log.it=FALSE,verbose = TRUE, maxit=Number_of_Iterations_For_Loess ,mat=as.matrix(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter))               ,ncol=dim(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter)[2],dimnames=dimnames(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter)))
####max iteration of cyclic loess is defined as 100
Homemades_R_Squared<-as.data.frame(t(apply((apply(as.matrix(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter),2,rank)),1,rank)))
Homemades_R_Cubed<-as.data.frame(apply(apply((apply(as.matrix(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter),2,rank)),1,rank),1,rank))

All_Normalized_Forms_of_Input_Data_Frame_as_list<-list(Input_Data_Frame_Usually_Peptide_Level_But_Does_Not_Really_Matter,Quantile,Global_Scaled,Variance_Stabalizing,Cubic_Spline,Cyclic_Loess_Smoothing,Homemades_R_Squared,Homemades_R_Cubed)
names(All_Normalized_Forms_of_Input_Data_Frame_as_list)<-c(Type_of_Input_Data_Frame,"..+ Quantile","..+ Global Scaled","..+ Variance Stab","..+ Cubic Sm Splines",(paste("..+ Cyclic Loess Itera=",as.character(Number_of_Iterations_For_Loess),sep="")),"..+ Spearman Homemade R^2","..+ Spearman Homemade R^3"
)
return(All_Normalized_Forms_of_Input_Data_Frame_as_list)

}                       				
##########################################################################################				
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################



##########################################################################################
				#Functions to Accomodate external packages

### This function performs a generic two class unpaired sig Pathway on all proteins and returns sigpathway specific data type
##Must be peptide level data
## Information on output file:
## Outputfile[1] is gsList, indexing file for peptide-protein
## Outputfile[2]$list.NTk has lists: ngs, nsim, t.set, t.set.new,p.null, p.value, q.value on Protein Level
## Outputfile[3]$list.NEk has same as above with different permutation
## Outputfile[4]$df.pathways contains a dataframe with ndexG,Gene Set Category, Pathway Set,Size,PercentUp,NTk,Stat# This is an important dataframe
## Outputfile[5]$list.gPS is a list of protein length and stratified results
## Outputfile[6]$parameters is a list that contains  [1] "nprobes"     "nsamples"    "phenotype"   "ngroups"     "minNPS"
## [6] "maxNPS"      "ngs"         "nsim.NTk"    "nsim.NEk"    "weightType"
##[11] "annotpkg"    "npath"       "allpathways"
## This is one way of extracting stats from their default run sig Pathway
## cbind(result_si[2]$list.NTk$p.value,  result_si[2]$list.NTk$q.value,    result_si[3]$list.NEk$p.value, result_si[3]$list.NEk$p.value, (result_si[4]$df.pathways[order(result_si[4]$df.pathways[,1]),])   )
## This binds the output[2] which as NTK pvals, [3] with NEK pvals, and correctly ordered results on [4], you can check by match the qvalues in both [2] and [4]



##################################Running sigPathway Internal functions
#Running GSEA
#wwww<-calculate.GSEA(gct_file.dataframe,Simple_Binary_Class,asdf[1]$gsList,nsim=1000,verbose=FALSE,alwaysUseRandomPerm=FALSE)
#wwwwww<-rankPathways.NGSk(wwww,sig_pathway_gmt_load_file,asdf[1]$gsList,methodName="GSEA",npath=670)
#Output from GSEA parsed head(cbind(result_si$p.value,result_si$q.value,result_si$t.set,result_si$t.set.new))
#rankPathways.NGSk(wwww,sig_pathway_gmt_load_file,asdf[1]$gsList,methodName="GSEA",npath=670)
#(rankPathways.NGSk(result_si,sig_pathway_gmt_load_file,output_from_sigPathway[1]$gsList,methodName="GSEA",npath=result_si$ngs)

hypothesis.testing.sig.pathway.default<-function(input_data_frame,Simple_Binary_Class,Number_of_Permutations){
set.seed(1234)
Number_of_proteins_by_peptide_in_df<-length((unique(which.proteins.do.these.peptides.belong.to(row.names(input_data_frame)))))
results<-runSigPathway(sig_pathway_gmt_load_file, 1, (Max_Peptide_Number+1), input_data_frame, Simple_Binary_Class, nsim = Number_of_Permutations,weightType = "constant", ngroups = 2,npath=Number_of_proteins_by_peptide_in_df,verbose = FALSE,allpathways = TRUE, alwaysUseRandomPerm = FALSE)
return(results)
}



hypothesis.testing.sig.pathway.GSEA<-function(input_data_frame,Simple_Binary_Class,number_of_permutations){
set.seed(1234)
Number_of_proteins_by_peptide_in_df<-length((unique(which.proteins.do.these.peptides.belong.to(row.names(input_data_frame)))))
#We run the deafault sig path because we need the native annotation file that's included in the result to run sig path gsea
#Hence the low number of suggested permutations
results<-hypothesis.testing.sig.pathway.default(input_data_frame,Simple_Binary_Class,10)
GSEA_Results_from_Sig<-calculate.GSEA(input_data_frame,Simple_Binary_Class,results[1]$gsList,nsim=number_of_permutations,verbose=FALSE,alwaysUseRandomPerm=FALSE)
#wwwwww<-rankPathways.NGSk(wwww,sig_pathway_gmt_load_file,asdf[1]$gsList,methodName="GSEA",npath=670)
return(GSEA_Results_from_Sig)
}








################################Custom Modified Functions for Hypothesis Testing##########
####################Specifically for p-value reporting

#####This function uses sigPathway native functions to calculated custom p q values 
# Initially desgined for two groups two sided t
#custom.hypothesis.t.peptide.permutation<-function(input_data_frame,phenotype_information,number_of_permutation)
#{
###T is positive when mean of group 1 is greater
#t_value_peptide<-calcTStatFast(input_data_frame,phenotype_information,ngroups=2)$tstat

#}
#########################################################################################
#Modified External Hypothesis Testing for only P-Value extraction





#Wrapper for all P- Value Summarization per data frame
g.go.p.calculate.all<-function(input_data_frame,Simple_Binary_Class,Number_of_Permutations){

All_procedures<-hypothesis.p.protein.level.sig.pathway.GSEA(input_data_frame,Simple_Binary_Class,Number_of_Permutations)[1]
All_procedures[2]<-hypothesis.p.protein.level.sig.pathway.NEk(input_data_frame,Simple_Binary_Class,Number_of_Permutations)[1]
All_procedures[3]<-hypothesis.p.protein.level.sig.pathway.NTk(input_data_frame,Simple_Binary_Class,Number_of_Permutations)[1]
All_procedures[4]<-hypothesis.p.peptide.level.min.pvalue.1(input_data_frame,Simple_Binary_Class)[1]
All_procedures[5]<-hypothesis.p.protein.simple.sum.1(input_data_frame,Simple_Binary_Class)[1]
All_procedures[6]<-hypothesis.p.protein.level.Tibshirani.GSA(input_data_frame,Simple_Binary_Class,Number_of_Permutations)[1]
All_procedures[7]<-hypothesis.p.peptide.level.samr.min.pvalue(input_data_frame,Simple_Binary_Class,Number_of_Permutations)[1]
All_procedures[8]<-hypothesis.p.protein.level.APPPt.Default.sample.perm(input_data_frame,Simple_Binary_Class,Number_of_Permutations,break_ties="average")[1]
return(All_procedures)
}




#################################################################################################################3
hypothesis.p.protein.level.APPPt.Default.sample.perm<-function(input_data_frame,Simple_Binary_Clas,Number_of_Permutations=1000,Type_of_Comaparison="Two class unpaired",Non_para="n",Summarize_of_peptide_t="sum",break_ties="first"){
### break_ties=first,average
Protein_level_ts_from_perms<-APT.Calculate.Master(input_data_frame,Simple_Binary_Clas,Type_of_Comaparison,Number_of_Permutations,Non_para,Summarize_of_peptide_t)
APPPT_PVALUE<-(        ((   apply(abs(Protein_level_ts_from_perms),1,rank,ties.method=break_ties))[(dim(Protein_level_ts_from_perms)[2]),]  )   /  (dim(Protein_level_ts_from_perms)[2])               )
APPPT_PVALUE<-(1-APPPT_PVALUE)
Protein_level_ps<-as.data.frame(APPPT_PVALUE)
names(Protein_level_ps)<-paste("APPPt_SmP",dim(Protein_level_ts_from_perms)[2]-1,Non_para,Summarize_of_peptide_t,break_ties,sep=":")
row.names(Protein_level_ps)<-row.names(Protein_level_ts_from_perms)
return(Protein_level_ps)
}

#################################################################################################################3
hypothesis.p.protein.level.APPPt.Default.scramble.perm<-function(input_data_frame,Simple_Binary_Clas,Number_of_Permutations=1000,Type_of_Comaparison="Two class unpaired",Non_para="n",Summarize_of_peptide_t="sum",break_ties="first"){
### break_ties=first,average
Protein_level_ts_from_perms<-APT.Calculate.Master.with.Simultaneous.Peptide.Sample.Perms(input_data_frame,Simple_Binary_Clas,Type_of_Comaparison,Number_of_Permutations,Non_para,Summarize_of_peptide_t)
APPPT_PVALUE<-(        ((   apply(abs(Protein_level_ts_from_perms),1,rank,ties.method=break_ties))[(dim(Protein_level_ts_from_perms)[2]),]  )   /  (dim(Protein_level_ts_from_perms)[2])               )
APPPT_PVALUE<-(1-APPPT_PVALUE)
Protein_level_ps<-as.data.frame(APPPT_PVALUE)
names(Protein_level_ps)<-paste("APPPt_Src",dim(Protein_level_ts_from_perms)[2]-1,Non_para,Summarize_of_peptide_t,break_ties,sep=":")
row.names(Protein_level_ps)<-row.names(Protein_level_ts_from_perms)
return(Protein_level_ps)
}
####################################################################################################################
####################################################################################################################
###################Scramble perm with mean instead of sum for t value summarization from peptide to protein#########
###################Average for breaking ties#######################################################################
############################and mean renormalization after scramble#####################################################
hypothesis.p.protein.level.APPPt.2.0.scramble.perm<-function(input_data_frame,Simple_Binary_Clas,Number_of_Permutations=1000,Type_of_Comaparison="Two class unpaired",Non_para="n",Summarize_of_peptide_t="mean",break_ties="average",normalize_after_each_scramble=FALSE){
### break_ties=first,average
Protein_level_ts_from_perms<-APT.Calculate.Master.with.Simultaneous.Peptide.Sample.Perms(input_data_frame,Simple_Binary_Clas,Type_of_Comaparison,Number_of_Permutations,Non_para,Summarize_of_peptide_t,normalize_after_each_scramble)
APPPT_PVALUE<-(        ((   apply(abs(Protein_level_ts_from_perms),1,rank,ties.method=break_ties))[(dim(Protein_level_ts_from_perms)[2]),]  )   /  (dim(Protein_level_ts_from_perms)[2])               )
APPPT_PVALUE<-(1-APPPT_PVALUE)
Protein_level_ps<-as.data.frame(APPPT_PVALUE)
names(Protein_level_ps)<-paste("APPPt_Src",dim(Protein_level_ts_from_perms)[2]-1,Non_para,Summarize_of_peptide_t,break_ties,sep=":")
row.names(Protein_level_ps)<-row.names(Protein_level_ts_from_perms)
return(Protein_level_ps)
}















#####################BETA VERIOSN OF GSA ###################
hypothesis.p.protein.level.Tibshirani.GSA<-function(input_data_frame,Simple_Binary_Class,Number_of_Permutations){


subset.GSA.peptide.set.list.for.subseted.analysis<-function(subsetted_data_frame)
{
Protein_names_relavant_to_subsetted_dataframe<-unique(which.proteins.do.these.peptides.belong.to(row.names(subsetted_data_frame)))
Peptide_set_list_relavant_for_this_data_frame<-Proteins$genesets[!is.na(match(Proteins$geneset.names,Protein_names_relavant_to_subsetted_dataframe))]

names(Peptide_set_list_relavant_for_this_data_frame)<-c(New_ordered_Peptide_list_Names<-Proteins$geneset.names[!is.na(match(Proteins$geneset.names,Protein_names_relavant_to_subsetted_dataframe))])
return(Peptide_set_list_relavant_for_this_data_frame)
}

#input_data_frame<-Mock_gct_file.dataframe
#Number_of_Permutations<-1000


GSA.obj<-GSA(x=input_data_frame,y=Simple_Binary_Class,method="maxmean", genenames=row.names(input_data_frame), genesets=subset.GSA.peptide.set.list.for.subseted.analysis(input_data_frame), resp.type="Two class unpaired",knn.neighbors=10,s0=NULL, s0.perc=NULL,minsize=1,maxsize=Max_Peptide_Number,restand=TRUE,nperms=Number_of_Permutations)


yellow<-GSA.listsets(GSA.obj,geneset.names = names(subset.GSA.peptide.set.list.for.subseted.analysis(input_data_frame)), maxchar=300,FDRcut=10000)

yellow<-data.frame(Tibshirani_GSA=(rbind(yellow$negative,yellow$positive)[,4]),   row.names=rbind(yellow$negative,yellow$positive)[,2]   )	




row_names_for_ordering_gsa_output<-hypothesis.p.peptide.level.min.pvalue.1(input_data_frame,Simple_Binary_Class)

answer<-data.frame(
Tibshirani_GSA=(merge(yellow,row_names_for_ordering_gsa_output,by=0)[,2]),
row.names=(merge(yellow,row_names_for_ordering_gsa_output,by=0)[,1]))

#p-values are multiplied by two because of the one sided results
return((confirm.data.for.numerical(answer))*2)


}

###########################################################





hypothesis.p.protein.level.sig.pathway.NTk<-function(input_data_frame,Simple_Binary_Class,Number_of_Permutations){
Sig_Pathway_Default_Results<-hypothesis.testing.sig.pathway.default(input_data_frame,Simple_Binary_Class,Number_of_Permutations)
Extracted_P_Values<-cbind(Sig_Pathway_Default_Results[2]$list.NTk$p.value,  (Sig_Pathway_Default_Results[4]$df.pathways[order(Sig_Pathway_Default_Results[4]$df.pathways[,1]),2])  )

#return(data.frame(Prot_LV_Perm_NKt_Tian_et_al=as.numeric(Extracted_P_Values[,1]),row.names=Extracted_P_Values[,2]))
Results.NTK<-data.frame(Tian_NKt=as.numeric(Extracted_P_Values[,1]),row.names=Extracted_P_Values[,2])
names(Results.NTK)<-paste("Tian_NTk",Sig_Pathway_Default_Results[6]$parameters$nsim.NTk,sep=":")
return(Results.NTK)
}



#Same as above, just extracting different columns for NEk rather than NKT
hypothesis.p.protein.level.sig.pathway.NEk<-function(input_data_frame,Simple_Binary_Class,Number_of_Permutations){
Sig_Pathway_Default_Results<-hypothesis.testing.sig.pathway.default(input_data_frame,Simple_Binary_Class,Number_of_Permutations)
Extracted_P_Values<-cbind(Sig_Pathway_Default_Results[3]$list.NEk$p.value,  (Sig_Pathway_Default_Results[4]$df.pathways[order(Sig_Pathway_Default_Results[4]$df.pathways[,1]),2])  )

Results.NEk<-data.frame(Tian_NEk=as.numeric(Extracted_P_Values[,1]),row.names=Extracted_P_Values[,2])
names(Results.NEk)<-paste("Tian_NEk",Sig_Pathway_Default_Results[6]$parameters$nsim.NEk,sep=":")
return(Results.NEk)
}


hypothesis.p.protein.level.sig.pathway.GSEA<-function(input_data_frame,Simple_Binary_Class,number_of_permutations){
set.seed(1234)
results_for_gene_lists<-hypothesis.testing.sig.pathway.default(input_data_frame,Simple_Binary_Class,10)
GSEA_Results_from_Sig<-calculate.GSEA(input_data_frame,Simple_Binary_Class,results_for_gene_lists[1]$gsList,nsim=number_of_permutations,verbose=FALSE,alwaysUseRandomPerm=FALSE)
GSEA_Results_from_Sig_Part_2<-rankPathways.NGSk(GSEA_Results_from_Sig,sig_pathway_gmt_load_file,results_for_gene_lists[1]$gsList,methodName="GSEA",results_for_gene_lists[6]$parameters[12]$npath)
##GSEA_Results_from_Sig_Part_2[(order(GSEA_Results_from_Sig_Part_2[,1])),2 ] <- this is actually the protein names
return(data.frame(Subramanium_GSEA=GSEA_Results_from_Sig$p.value,row.names=(GSEA_Results_from_Sig_Part_2[(order(GSEA_Results_from_Sig_Part_2[,1])),2 ])))

}































##########################################################################################
#This function computes p-value for two group two sided t-test un equal variances
#Accpets a data_frame and returns data.frame with 1 column for p-values
hypothesis.p.row.level.simple.1<-function(input_data_frame,phenotype_information)
{
return(data.frame(Raw_P_tTest_Unequal_Var_2_Tailed=(calcTStatFast(input_data_frame,phenotype_information,ngroups=2)$pval),row.names=rownames(input_data_frame)))
}
#This function computes p-value for two unpaired group two sided t-test on the protein level
#Accepts a data_frame with peptide row names and summarizes on protein level
#Averaging the peptide intensities on the protein level and then computes a t-test
hypothesis.p.protein.simple.average.1<-function(input_data_frame,phenotype_information){
data_set_storing_peptide_intensities_summarized<-summarize.peptide.data.frame.on.protein.by.averaging(input_data_frame)


changed_column_name="Prot_LV_Pept_Intsty_Mean"
temporary_data_frame<-data.frame(hypothesis.p.row.level.simple.1(data_set_storing_peptide_intensities_summarized,phenotype_information))
names(temporary_data_frame)<-changed_column_name
return(temporary_data_frame)
}
#This function computes p-value for two unpaired group two sided t-test on the protein level
#Accepts a data_frame with peptide row names and summarizes on protein level
#Summing the peptide intensities on the protein level and then computes a t-test
hypothesis.p.protein.simple.sum.1<-function(input_data_frame,phenotype_information){
data_set_storing_peptide_intensities_summarized<-summarize.peptide.data.frame.on.protein.by.summing(input_data_frame)
changed_column_name="RSmith_PNNL"
temporary_data_frame<-data.frame(hypothesis.p.row.level.simple.1(data_set_storing_peptide_intensities_summarized,phenotype_information))
names(temporary_data_frame)<-changed_column_name
return(temporary_data_frame)
}


#This function computes p-value for two unpaired group two sided t-test on the peptide level than summarizes
# on the protein level by calculating the min of p-values
#Accepts a data_frame with peptide row names and summarizes on protein level
#To keep things consistent and keep stick with data frames, we have to trick the min function to think our p-value
# dataframe has two columns, otherwise the the rownames keep getting lost as R converts it to numeric, since it has 1 column
# This scenario is different from the above two because the peptide to summarization is done before the p-value on the
#above functions where the dataset still has more than one column
hypothesis.p.peptide.level.min.pvalue.1<-function(input_data_frame,phenotype_information){
temporary_data_frame<-data.frame(hypothesis.p.row.level.simple.1(input_data_frame,phenotype_information))

temporary_data_frame<-summarize.peptide.data.frame.on.protein.by.min(cbind(temporary_data_frame,temporary_data_frame))
##Here we get rid of the extra column and declare as a data.frame once again because of reason above
temporary_data_frame<-data.frame(temporary_data_frame[,1],row.names=row.names(temporary_data_frame))
changed_column_name="Common_Pep"
names(temporary_data_frame)<-changed_column_name
return(temporary_data_frame)
}

#This function computes p-value for two unpaired group two sided t-test on the peptide level 
# using sam from samr package then summarizes
# on the protein level by calculating the min of p-values
#Accepts a data_frame with peptide row names and summarizes on protein level
#To keep things consistent and keep stick with data frames, we have to trick the min function to think our p-value
# dataframe has two columns, otherwise the the rownames keep getting lost as R converts it to numeric, since it has 1 column
# This scenario is different from the above two because the peptide to summarization is done before the p-value on the
#above functions where the dataset still has more than one column
hypothesis.p.peptide.level.samr.min.pvalue<-function(input_data_frame,Simple_Binary_Class,number_of_permutations){
Sam_Obj_from_samr<-samr(list(x=as.matrix(input_data_frame),y=Simple_Binary_Class,logged2=TRUE),resp.type="Two class unpaired", nperms=number_of_permutations)
P_Values_From_Sam<-samr.pvalues.from.perms(Sam_Obj_from_samr$tt, Sam_Obj_from_samr$ttstar)
temporary_data_frame<-as.data.frame(as.numeric(P_Values_From_Sam))
row.names(temporary_data_frame)<-row.names(input_data_frame)

temporary_data_frame<-summarize.peptide.data.frame.on.protein.by.min(cbind(temporary_data_frame,temporary_data_frame))
##Here we get rid of the extra column and declare as a data.frame once again because of reason above
temporary_data_frame<-data.frame(temporary_data_frame[,1],row.names=row.names(temporary_data_frame))
changed_column_name="Sam_Pep_Min"
names(temporary_data_frame)<-changed_column_name
return(temporary_data_frame)
}










###################################################################################################
#This function summarizes peptide level info based on protein information, excepts any #of cols 
# But only a data.frame with peptides as row.names Returns data.frame on protein level with proteins as row.names
summarize.peptide.data.frame.on.protein.by.summing<-function(input_data_frame){
#Extract the redundant list of protein of the peptides using pre-defined function
redundant_protein_char<-which.proteins.do.these.peptides.belong.to(as.character(rownames(input_data_frame)))
#created new dataframe with redundant proteins info and numeric assigment to each unique protein
temp_data_frame_with_protein<-cbind((as.numeric(as.factor(redundant_protein_char))),redundant_protein_char,input_data_frame)
#The number of unique proteins are calculated by the calculating the max of the number assigned earlier now concatinated to a new data.frame
number_of_unique_proteins<-(max(as.numeric(temp_data_frame_with_protein[,1])))# This column stores the number assigned
#Now that we have the number of unique proteins we make a dataframe to store summarized info and reduce the number of row accordingly
Data_frame_to_store_protein_level_info<-input_data_frame[(1:number_of_unique_proteins),] # Create new dataframe and then reduce
for(i in 1:number_of_unique_proteins)
{
	temp_rows_info_for_protein<-input_data_frame[(temp_data_frame_with_protein[,1]==i),]
##########################Summarization by summing "colSum" method of rows
#############################
###TO SUM COLUMNS:-
Data_frame_to_store_protein_level_info[i,]<-(colSums(temp_rows_info_for_protein))
###TO MEAN COLUMS:-
#Data_frame_to_store_protein_level_info[i,]<-(colSums(temp_rows_info_for_protein))
###
##############################
#Assigning the New Dataframe the new values and Protein name for row name
#Get the protein infor from the 2 column since we concatinated it with redundant protein ids ealier
#row.names(Data_frame_to_store_protein_level_info[1,])<-as.character(unique(temp_data_frame_with_protein[temp_data_frame_with_protein[,1]==i,2]))
row.names(Data_frame_to_store_protein_level_info)[i]<-(as.character(unique(temp_data_frame_with_protein[(temp_data_frame_with_protein[,1]==i),2])))
}
return(Data_frame_to_store_protein_level_info)
}
###################################################################################################
#This function summarizes(mean) peptide level info based on protein information, excepts  >2of cols
#But only a data.frame with peptides as row.names:Identical to one previous but less information
summarize.peptide.data.frame.on.protein.by.averaging<-function(input_data_frame){
redundant_protein_char<-which.proteins.do.these.peptides.belong.to(as.character(rownames(input_data_frame)))
temp_data_frame_with_protein<-cbind((as.numeric(as.factor(redundant_protein_char))),redundant_protein_char,input_data_frame)
number_of_unique_proteins<-(max(as.numeric(temp_data_frame_with_protein[,1])))# This column stores the number assigned
Data_frame_to_store_protein_level_info<-input_data_frame[(1:number_of_unique_proteins),] # Create new dataframe and then reduce
for(i in 1:number_of_unique_proteins)
{
        temp_rows_info_for_protein<-input_data_frame[(temp_data_frame_with_protein[,1]==i),]
Data_frame_to_store_protein_level_info[i,]<-(mean(temp_rows_info_for_protein))
row.names(Data_frame_to_store_protein_level_info)[i]<-(as.character(unique(temp_data_frame_with_protein[(temp_data_frame_with_protein[,1]==i),2])))
}
return(Data_frame_to_store_protein_level_info)
}
###################################################################################################
###################################################################################################

###################################################################################################
###################################################################################################
summarize.peptide.data.frame.on.protein.by.peptide.cor<-function(input_data_frame){
####method computes pearson and spearmans correlation on peptides within a protein
#### then calculates the average of the values in the cor matrix after removing the diagonal 1s(actually subtracting from total) 
####
redundant_protein_char<-which.proteins.do.these.peptides.belong.to(as.character(rownames(input_data_frame)))
temp_data_frame_with_protein<-cbind((as.numeric(as.factor(redundant_protein_char))),redundant_protein_char,input_data_frame)
number_of_unique_proteins<-(max(as.numeric(temp_data_frame_with_protein[,1])))# This column stores the number assigned
Data_frame_to_store_protein_level_info<-input_data_frame[(1:number_of_unique_proteins),1:2] # Create new dataframe and then reduce
for(i in 1:number_of_unique_proteins)
{
        temp_rows_info_for_protein<-input_data_frame[(temp_data_frame_with_protein[,1]==i),]
#Data_frame_to_store_protein_level_info[i,]<-as.data.frame(t(scores(pca(t(temp_rows_info_for_protein),nPcs=1,method="svd"))))
Data_frame_to_store_protein_level_info[i,]<-

c(

 ((sum(cor(t(temp_rows_info_for_protein), method="spearman"))- dim(temp_rows_info_for_protein)[1]) / (  ((dim(temp_rows_info_for_protein)[1])^2) -   dim(temp_rows_info_for_protein)[1]))             ,
 ((sum(cor(t(temp_rows_info_for_protein), method="pearson"))- dim(temp_rows_info_for_protein)[1]) / (  ((dim(temp_rows_info_for_protein)[1])^2) -   dim(temp_rows_info_for_protein)[1]))
 
 

)
# c(
# sum(cor(t(temp_rows_info_for_protein),method="spearman")),    sum(cor(t(temp_rows_info_for_protein),method="pearson"))
# )

row.names(Data_frame_to_store_protein_level_info)[i]<-(as.character(unique(temp_data_frame_with_protein[(temp_data_frame_with_protein[,1]==i),2])))
names(Data_frame_to_store_protein_level_info)<-c("Spearmans_RHO_AVG","Pearsons_RHO_AVG")
}
return(Data_frame_to_store_protein_level_info)
}


###################################################################################################
#This function summarizes peptide level info based on protein information, excepts  >2of cols
#But only a data.frame with peptides as row.names:Identical to one previous but less information
summarize.peptide.data.frame.on.protein.by.1.pca<-function(input_data_frame,method="svd"){
####method can equal one of the following, by default is "svd"
#### pca, does restandardization before computation.
#### method ="nipals" good for missing values, "svd","bpca"
redundant_protein_char<-which.proteins.do.these.peptides.belong.to(as.character(rownames(input_data_frame)))
temp_data_frame_with_protein<-cbind((as.numeric(as.factor(redundant_protein_char))),redundant_protein_char,input_data_frame)
number_of_unique_proteins<-(max(as.numeric(temp_data_frame_with_protein[,1])))# This column stores the number assigned
Data_frame_to_store_protein_level_info<-input_data_frame[(1:number_of_unique_proteins),] # Create new dataframe and then reduce
for(i in 1:number_of_unique_proteins)
{
        temp_rows_info_for_protein<-input_data_frame[(temp_data_frame_with_protein[,1]==i),]
Data_frame_to_store_protein_level_info[i,]<-as.data.frame(t(scores(pca(t(temp_rows_info_for_protein),nPcs=1,method="svd"))))
row.names(Data_frame_to_store_protein_level_info)[i]<-(as.character(unique(temp_data_frame_with_protein[(temp_data_frame_with_protein[,1]==i),2])))
}
return(Data_frame_to_store_protein_level_info)
}
###################################################################################################
###################################################################################################


#This function summarizes peptide level p-values info based on protein information, excepts >2 cols
#But only a data.frame with peptides as row.names:Identical to one previous but this one is primary
# for summarizing by calculating the min more suitable for p-values and nagative t scores rather than intensities
summarize.peptide.data.frame.on.protein.by.min<-function(input_data_frame){
redundant_protein_char<-which.proteins.do.these.peptides.belong.to(as.character(rownames(input_data_frame)))
temp_data_frame_with_protein<-cbind((as.numeric(as.factor(redundant_protein_char))),redundant_protein_char,input_data_frame)
number_of_unique_proteins<-(max(as.numeric(temp_data_frame_with_protein[,1])))# This column stores the number assigned
Data_frame_to_store_protein_level_info<-input_data_frame[(1:number_of_unique_proteins),] # Create new dataframe and then reduce
for(i in 1:number_of_unique_proteins)
{
        temp_rows_info_for_protein<-input_data_frame[(temp_data_frame_with_protein[,1]==i),]
Data_frame_to_store_protein_level_info[i,]<-(apply(temp_rows_info_for_protein,2,min))
row.names(Data_frame_to_store_protein_level_info)[i]<-(as.character(unique(temp_data_frame_with_protein[(temp_data_frame_with_protein[,1]==i),2])))
}
return(Data_frame_to_store_protein_level_info)
}
###################################################################################################









###################################################################################################

				#Functions to Alter Data (Permute etc)
##########################################################################################
### This Function excepts a data.frame and randomizes the peptides, returns data.frame 

row.randomizer.data.frame<-function(original.data_frame){
number_of_peptides<-dim(original.data_frame)[1]
random_list<-as.data.frame(c(rnorm(number_of_peptides,mean=20000,sd=300)))
permuted.data_frame<-original.data_frame[order(random_list[,1]),]
row.names(permuted.data_frame)<-row.names(original.data_frame)
return(permuted.data_frame)
}


###  This Function excepts a data.frame and randomizes the classes, returns data.frame 
class.randomizer.data.frame<-function(original.data_frame){
number_of_classes<-dim(original.data_frame)[2]
random_list<-as.data.frame(c(rnorm(number_of_classes,mean=20000,sd=300)))
permuted.data_frame<-original.data_frame[,order(random_list[,1])]
names(permuted.data_frame)<-names(original.data_frame)
return(permuted.data_frame)
}


###  This Function excepts a data.frame and scrambles all data, returns data.frame 
scramble.data.frame<-function(original.data_frame){
number_of_peptides<-dim(original.data_frame)[1]
number_of_classes<-dim(original.data_frame)[2]
one_dimensional_data_storage<-as.vector(as.matrix(original.data_frame))
random_order_data_storage<-as.vector(rnorm((number_of_peptides*number_of_classes),24344,sd=45.6))
two_dimensional_data_frame<-as.data.frame(cbind(random_order_data_storage,one_dimensional_data_storage))
scrambled_data_frame<-as.data.frame(matrix((two_dimensional_data_frame[order(two_dimensional_data_frame[,1]),2]),number_of_peptides,ncol=number_of_classes))
names(scrambled_data_frame)<-names(original.data_frame)
row.names(scrambled_data_frame)<-row.names(original.data_frame)
return(scrambled_data_frame)
}


###  This Function excepts a data.frame and creates a Mock Data set with their respective options, returns data.frame

## Creates a Data.frame with the same overall mean and sd normaly distributed
create.simple.normal.mock.data.frame<-function(original.data_frame){
number_of_peptides<-dim(original.data_frame)[1]
number_of_classes<-dim(original.data_frame)[2]
Input_Dataset_Mean<-mean(mean(original.data_frame,na.rm=TRUE))
Input_Dataset_Sd<-mean(sd(original.data_frame,na.rm=TRUE))
Mock_Dataset<-as.data.frame(matrix(c((rnorm((number_of_peptides*number_of_classes),mean=Input_Dataset_Mean, sd=Input_Dataset_Sd))),number_of_peptides,ncol=number_of_classes))
names(Mock_Dataset)<-names(original.data_frame)
row.names(Mock_Dataset)<-row.names(original.data_frame)
return(Mock_Dataset)
}


create.simple.gamma.mock.data.frame<-function(original.data_frame){
number_of_peptides<-dim(original.data_frame)[1]
number_of_classes<-dim(original.data_frame)[2]
Input_Dataset_Mean<-mean(mean(original.data_frame,na.rm=TRUE))
Input_Dataset_Sd<-mean(sd(original.data_frame,na.rm=TRUE))
Mock_Dataset<-as.data.frame(matrix(c((rgamma((number_of_peptides*number_of_classes),shape=100,scale=Input_Dataset_Mean))),number_of_peptides,ncol=number_of_classes))
names(Mock_Dataset)<-names(original.data_frame)
row.names(Mock_Dataset)<-row.names(original.data_frame)
return(Mock_Dataset)
}



create.variable.gamma.mock.data.frame<-function(original.data_frame,shape){
number_of_peptides<-dim(original.data_frame)[1]
number_of_classes<-dim(original.data_frame)[2]
Input_Dataset_Mean<-mean(mean(original.data_frame,na.rm=TRUE))
Input_Dataset_Sd<-mean(sd(original.data_frame,na.rm=TRUE))
Mock_Dataset<-as.data.frame(matrix(c((rgamma((number_of_peptides*number_of_classes),shape=shape,scale=Input_Dataset_Mean))),number_of_peptides,ncol=number_of_classes))
names(Mock_Dataset)<-names(original.data_frame)
row.names(Mock_Dataset)<-row.names(original.data_frame)
return(Mock_Dataset)
}


########################################Graphing Functions################################
########################################Graphing Functions################################
########################################Graphing Functions################################
########################################Graphing Functions################################
########################################Graphing Functions################################
########################################Graphing Functions################################
#This function takes in a dataframe, a character of group names == to column number, and plots ma plots for columns with a group to a pdf
graph.ma.plots.for.all.comparison.within.group<-function(Input_Data_Frame_,Character_Class_Designation_Groups,Add_Character_To_Title_of_Graph){
Number_of_Levels_<-levels(as.factor(Character_Class_Designation))
for(i in 1:length(Number_of_Levels_)){
mva.pairs(Input_Data_Frame_[,(Character_Class_Designation_Groups==Number_of_Levels_[i])],log.it=FALSE,cex=1.5,main=paste("MA Plots of ",Number_of_Levels_[i],Add_Character_To_Title_of_Graph))
}
}


#####This function graphs the output from Pre.Processing_S_3_All_Normalizations
##### Assumes at least 6 different dataframes in list , but can be modified below

graph.output.from.Pre.Processing_S_3_All_Normalizations<-function(List_of_Data_Frames){
cols<-brewer.pal((dim(List_of_Data_Frames[[1]])[2]),"Set3")
Types_of_Pre_Processing<-names(List_of_Data_Frames)
par(mfrow=c(2,4))
boxplot(List_of_Data_Frames[[1]],col=cols,main=Types_of_Pre_Processing[1], las=2,cex=0.8)
boxplot(List_of_Data_Frames[[7]],col=cols,main=Types_of_Pre_Processing[7], las=2,cex=0.8)
boxplot(List_of_Data_Frames[[8]],col=cols,main=Types_of_Pre_Processing[8], las=2,cex=0.8)
boxplot(List_of_Data_Frames[[3]],col=cols,main=Types_of_Pre_Processing[3], las=2,cex=0.8)
meanSdPlot(as.matrix(List_of_Data_Frames[[1]]),col=cols,las=2,cex=0.8)
meanSdPlot(as.matrix(List_of_Data_Frames[[7]]),col=cols,las=2,cex=0.8)
meanSdPlot(as.matrix(List_of_Data_Frames[[8]]),col=cols,las=2,cex=0.8)
meanSdPlot(as.matrix(List_of_Data_Frames[[3]]),col=cols,las=2,cex=0.8)


boxplot(List_of_Data_Frames[[4]],col=cols,main=Types_of_Pre_Processing[4], las=2,cex=0.8)
boxplot(List_of_Data_Frames[[6]],col=cols,main=Types_of_Pre_Processing[6], las=2,cex=0.8)
boxplot(List_of_Data_Frames[[5]],col=cols,main=Types_of_Pre_Processing[5], las=2,cex=0.8)
boxplot(List_of_Data_Frames[[2]],col=cols,main=Types_of_Pre_Processing[2], las=2,cex=0.8)
meanSdPlot(as.matrix(List_of_Data_Frames[[4]]),col=cols,las=2,cex=0.8)
meanSdPlot(as.matrix(List_of_Data_Frames[[6]]),col=cols,las=2,cex=0.8)
meanSdPlot(as.matrix(List_of_Data_Frames[[5]]),col=cols,las=2,cex=0.8)
meanSdPlot(as.matrix(List_of_Data_Frames[[2]]),col=cols,las=2,cex=0.8)
par(mfrow=c(1,1))
}
#########################################################################################
graph.heatmap.peptide.data.frame<-function(input_peptide_data_frame,main_title=""){
heatmap.2((as.matrix(input_peptide_data_frame)),las=1,dendrogram="col",main=main_title,cexRow=0.5,cexCol=1,tracecol="black",labRow=NA, density.info="density",col=bluered,breaks=160,
ColSideColors=as.character(  (rainbow(60)[c(6,16,29,11,3,36,46)])[Simple_Binary_Class]  ),
colsep=c(1:length(Simple_Binary_Class)),sepcolor="white",sepwidth=c(0.05,0),xlab="Sample",ylab="Peptide")
#colsep=c(1:length(Simple_Binary_Class)),sepcolor="white",sepwidth=c(0.05,0),xlab="Sample",ylab="Peptide",scale="column")
}      

#############################################################################################

graph.missing.data.basic.statistics.for.data.frame.within.group<-function(INPUT_DATA_FRAME_ , simple_binary_class=Simple_Binary_Class,character_class_designation=Character_Class_Designation,na_string="NA")
{
cols<-rainbow(dim(INPUT_DATA_FRAME_)[2])
INPUT_DATA_FRAME_[INPUT_DATA_FRAME_==na_string]<-"NA"
confirm.data.for.numerical(INPUT_DATA_FRAME_)->INPUT_DATA_FRAME_
par(mfrow=c(1,3))
barplot(calculate.missing.value.percentage.for.input.data.frame(INPUT_DATA_FRAME_),col=cols,ylab="% of Missing Values",las=2,axes=TRUE,cex.names=0.8)
# plot any correlation
Associate_Values_Per_Group_For_Missing_Values<-calculate.missing.value.associated.abundace.for.input.data.frame(INPUT_DATA_FRAME_)
Group_Number_for_plot=1
hist(Associate_Values_Per_Group_For_Missing_Values[is.na(apply(Associate_Values_Per_Group_For_Missing_Values[,(simple_binary_class==Group_Number_for_plot)], 1, mean)), (simple_binary_class==Group_Number_for_plot)],col=cols[Group_Number_for_plot],xlab="Percentile of Non-Missing Adjacent Data",cex=0.4,,xlim=c(0,100),main=(unique(character_class_designation)[Group_Number_for_plot]))
Group_Number_for_plot=2
hist(Associate_Values_Per_Group_For_Missing_Values[is.na(apply(Associate_Values_Per_Group_For_Missing_Values[,(simple_binary_class==Group_Number_for_plot)], 1, mean)), (simple_binary_class==Group_Number_for_plot)],col=cols[Group_Number_for_plot],xlab="Percentile of Non-Missing Adjacent Data",cex=0.4,xlim=c(0,100),main=(unique(character_class_designation)[Group_Number_for_plot]))

par(mfrow=c(1,1))
}

#########################################################################3

#########################################################################################
graph.whole.data.frame.simple<-function(input.data_frame){
x11()
par(mfrow = c(1, 3))
hist(as.matrix(input.data_frame),main=paste("Histogram",""),xlab="Abundance")
stripchart(as.matrix(input.data_frame),add=TRUE,at=15.5)
qqnorm(as.matrix(input.data_frame))
qqline(input.data_frame)

boxplot(input.data_frame,las=2)

}











quickly.plot.p.values.from.all.methods<-function(input_data_frame,Simple_Binary_Class,Number_of_Permutations){
answer<-g.go.p.calculate.all(input_data_frame,Simple_Binary_Class,Number_of_Permutations)
Number_of_Peptide_for_Protein<-how.many.peptides.belong.to.each.of.these.proteins(row.names(answer))
answer<-cbind(answer,Number_of_Peptide_for_Protein)
x11()
stripplot( Prot_LV_Perm_GSEA_Subramanium_et_al ~ factor(Number_of_Peptide_for_Protein), answer, 
           jitter.data = TRUE, alpha = 0.6,
           xlab = "Number of Peptides for Protein", ylab = "P-Values")
x11()
stripplot( P_Pep_LV_Min_P_Protein_Classical ~ factor(Number_of_Peptide_for_Protein), answer, 
           jitter.data = TRUE, alpha = 0.6,
           xlab = "Number of Peptides for Protein", ylab = "P-Values")

}

########################################################################################




###################################Functions For Loading Data etc####################


### This function allows you to visually confirm correct comparison variable designation for tests and graphs
report.global.variables.details<-function()
{

print(paste("Binary_Simple_Class_Binary",toString(Binary_Simple_Class_Binary),sep="         ...         "))
print(paste("Binary_Simple_Class_Binary",toString(Binary_Simple_Class_Binary),sep="         ...         "))
print(paste("Character_Type_of_Class_Comparison",toString(Character_Type_of_Class_Comparison),sep="         ...         "))
print(paste("Character_Class_Designation",toString(Character_Class_Designation),sep="         ...         "))
print(paste("gct_file",toString(gct_file),sep="         ...         "))
print(paste("gmt_file",toString(gmt_file),sep="         ...         "))
print(paste("redundant_proteins_file_name",toString(redundant_proteins_file_name),sep="         ...         "))
print(paste("Number_of_Peptides",toString(Number_of_Peptides),sep="         ...         "))
print(paste("Number_of_Proteins",toString(Number_of_Proteins),sep="         ...         "))
print(paste("Max_Peptide_Number",toString(Max_Peptide_Number),sep="         ...         "))
#print(paste("",toString(),sep="         ...         "))
}


#Splits the Column Header for gct file by character "_" and attempts to convert to 2s and 1s if binary class designation
#Designed to print desgination for confirmation, accepts data.frame, returns 1 dimensional numeric
extract.class.designatio.of.binary.designation<-function(original.data_frame){
extracted_list_for_designation<-strsplit(names(original.data_frame),split="_")
single_binary_factors<-factor(as.character(as.matrix((as.data.frame(extracted_list_for_designation)[1,]))))
print(cbind(names(original.data_frame), as.numeric(single_binary_factors)))
return(as.numeric(single_binary_factors))
}


#Reads in the GSA gmt native file, and determines the protein and corresponding number of peptides
# NOTE: right now there is always one more peptide than should be because of the unix script that 
# adds one extra tab at the end of every gmt file
extract.protein.name.and.number.of.peptides<-function(Proteins_the_GSA_native_data_format){
Number_of_Proteins<-length(Proteins_the_GSA_native_data_format$genesets)
Proteins_and_Peptide_Numbers<-data.frame(Number_of_Peptides=c(rep(0,Number_of_Proteins)))
row.names(Proteins_and_Peptide_Numbers)<-Proteins_the_GSA_native_data_format$geneset.names
for(i in 1:Number_of_Proteins){
Proteins_and_Peptide_Numbers[i,1]<-length(Proteins_the_GSA_native_data_format$genesets[[i]])
}
return(Proteins_and_Peptide_Numbers)
}

## This function return character list of Proteins with a given threshold of peptides from original list
## Accepts a number and outputs a character type
which.proteins.have.atleast.x.peptides<-function(min_number_of_peptides){
return(row.names(Proteins_and_Peptide_Numbers)[Proteins_and_Peptide_Numbers>(min_number_of_peptides-1)])
}

## This function returns a character list of peptides that belong to proteins from original list
## Accepts character list i.e. outputted from which.proteins.have.at...
which.peptides.belong.to.these.proteins<-function(list_of_proteins_of_character_type){
return(Peptides[!is.na(match(as.character(Redundant_Proteins[,1]),list_of_proteins_of_character_type))])
}


## This Functions returns a as.character type of redundant proteins matching the peptides
#used as the input
which.proteins.do.these.peptides.belong.to<-function(list_of_peptides_of_character_type){
return(as.character(Redundant_Proteins[!is.na(match(Peptides,list_of_peptides_of_character_type)),]))
}

## This functions returns a numeric vector of equal length of the inputed character list of proteins
## it returns the number of peptides for each of the proteins in the list
how.many.peptides.belong.to.each.of.these.proteins<-function(list_of_proteins_of_character_type){
Temp_Number_Recorder<-0
for(i in 1:length(list_of_proteins_of_character_type)){
Temp_Number_Recorder[i]<-length(which.peptides.belong.to.these.proteins(list_of_proteins_of_character_type[i]))
}
return(as.numeric(Temp_Number_Recorder))
}


#Utility Function to make sure the data frame has numerics instead of factors, happens when the dataframe is on
#column and the column has identical values it's considered to bea  factor with levels < number of rows
# Saw problems plottinf the conglomerate dataframe of pvalues from different t-tests (perm)
confirm.data.for.numerical<-function(input_data_frame){
temp.var<-data.frame(apply(apply(input_data_frame,2,as.character),2,as.numeric))
row.names(temp.var)<-row.names(input_data_frame)
colnames(temp.var)<-names(input_data_frame)
return(temp.var)
}








## This function accepts a dataframe with peptide as rownames, and creates a subset of it by filtering for
## only those peptides that meet the criteria of belong to a protein that have atleast the # peptides
## Accepts a number, and data.frame and outputs data.frame
subset.peptides.row.name.data.frame.that.belong.to.a.protein.with.atleast.x.peptides<-function(min_number_of_peptides,input_data_frame){
all_peptides_with_criteria_from_original_data_set<-which.peptides.belong.to.these.proteins(which.proteins.have.atleast.x.peptides(min_number_of_peptides))
return(input_data_frame[!is.na(match(row.names(input_data_frame),all_peptides_with_criteria_from_original_data_set)),])
}

## This function accepts a dataframe with protein as rownames, and creates a subset of it by filtering for
## only those proteins that meet the criteria of having atleast the # peptides inputted
## Accepts a number, and data.frame and outputs data.frame
subset.protein.row.name.data.frame.with.proteins.with.least.x.peptides<-function(min_number_of_peptides,input_data_frame){
all_proteins_with_criteria_from_original_data_set<-which.proteins.have.atleast.x.peptides(min_number_of_peptides)
##Below the match only matches those rows with the proteins that match the proteins that have met the criteria
return(input_data_frame[!is.na(match(row.names(input_data_frame),all_proteins_with_criteria_from_original_data_set)),])
}



## This function accepts a dataframe as input with NAs (usually NAs are identified by a 0), but this is done before the dataframe comes into
# this function. The function then computes the percentage of missing values per column for the inputted data set.
# returns a dataframe.
calculate.missing.value.percentage.for.input.data.frame<-function(input_data_frame){
missing_values_percentage_per_sample<-(apply((apply(input_data_frame,2,is.na)),2,sum)*100)/dim(input_data_frame)[1]
return(missing_values_percentage_per_sample)
}


## This function accepts a dataframe as input with NAs (usually NAs are identified by a 0), but this is done before the dataframe comes into
# this function. The function then returns the percentile of the non-missing data adjacent to the missing data in the same group
#Returns a data.frame
calculate.missing.value.associated.abundace.for.input.data.frame<-function(Input_Data_Frame){
# Order the abundaces by rank, ignore missing values
Input_data_set_missing_associated_ranks<-apply(Input_Data_Frame,2,rank,na.last="keep")
Maxes_for_each_column<-apply(Input_data_set_missing_associated_ranks,2,max,na.rm=TRUE)
for(i in 1:(length(Maxes_for_each_column)))
{
	Input_data_set_missing_associated_ranks[,i]<-Input_data_set_missing_associated_ranks[,i]/Maxes_for_each_column[i]
}
Input_data_set_missing_associated_ranks<-Input_data_set_missing_associated_ranks*100
Extract_all_rows_with_missing_values<-Input_data_set_missing_associated_ranks[(is.na(apply(Input_data_set_missing_associated_ranks,1,mean))),]
#Get the unique class designations for calulating the mean percentile of abundance for missing value rows 
groups_designations_unique<-unique(Simple_Binary_Class)
Extract_all_rows_with_missing_values
#data_for_graphs<-Extract_per_group[(is.na(apply(Extract_per_group,1,mean))),]
return(Extract_all_rows_with_missing_values)
}



##This function returns the imputed data frame using bpca
## Several Alternatives are also listed 
##Accepts Dataframe returns data frame
impute.data.frame.using.pca<-function(input_data_frame,method,colmax=0.5,row.max=0.5,k=10,rng.seed=345){
## Methods: “svd”, “nipals”, “bpca”, “ppca”, “svdImpute”, “nlpca”, “robustPca”
##Reference for Imputations Methods:
#http://www.bioconductor.org/packages/2.2/bioc/vignettes/pcaMethods/inst/doc/pcaMethods.pdf
##Reference of Bpca outperforming knn:
#Troyanskaya O. and Cantor M. and Sherlock G. and Brown P. and Hastie T. and Tibshirani R. and Botstein D. and Altman RB. Missing
#value estimation methods for DNA microarrays. Bioinformatics. 2001
imputed_data_frame <- as.data.frame(pca(as.matrix(input_data_frame), nPcs = 5, method = method)@completeObs)
return(imputed_data_frame)
}

Impute.data.frame.using.x.method<-function(input_data_frame,method,nPcs=5,colmax=0.5,row.max=0.5,k=10,rng.seed=345){
## Methods: “svd”, “nipals”, “bpca”, “ppca”, “svdImpute”, “nlpca”, “robustPca”,"knn"
if(method=="knn"){
return(as.data.frame(impute.knn(as.matrix(input_data_frame),k, row.max, colmax, maxp = 1500, rng.seed=rng.seed)$data))
}
imputed_data_frame <- as.data.frame(pca(as.matrix(input_data_frame), nPcs = 5, method = method)@completeObs)
return(imputed_data_frame)
}




calculate_t_scores_for_peptides_and_report_on_protein_level<-function(input_data_frame,Simple_Binary_Class){
T_values_for_peptides<-as.data.frame(mt.teststat(input_data_frame,(Simple_Binary_Class)-1,test="t",nonpara="n"))
row.names(T_values_for_peptides)<-row.names(input_data_frame)
Number_of_Proteins_Represented<-length(unique(which.proteins.do.these.peptides.belong.to(row.names(input_data_frame))))
Data_Frame_for_T_Values_For_Proteins<-as.data.frame((matrix(data=(rep("NA",Number_of_Proteins_Represented*Max_Peptide_Number)),nrow=Number_of_Proteins_Represented,ncol=Max_Peptide_Number,byrow=FALSE)))
Data_Frame_for_T_Values_For_Proteins<-as.data.frame(apply(Data_Frame_for_T_Values_For_Proteins,2,as.numeric))
row.names(Data_Frame_for_T_Values_For_Proteins)<-unique(which.proteins.do.these.peptides.belong.to(row.names(input_data_frame)))
for(i in 1:(dim(Data_Frame_for_T_Values_For_Proteins)[1])){
Protein=row.names(Data_Frame_for_T_Values_For_Proteins[i,])
T_Values<-extract.peptide.level.singular.value.for.a.given.protein(Protein,T_values_for_peptides)
Data_Frame_for_T_Values_For_Proteins[i,1:length(T_Values)]<-T_Values
}

###Order the dataframe such that it is ordered by the row name to maintain consistency with other functions
Data_Frame_for_T_Values_For_Proteins<-Data_Frame_for_T_Values_For_Proteins[order(row.names(Data_Frame_for_T_Values_For_Proteins)),]
return(Data_Frame_for_T_Values_For_Proteins)
}
extract.peptide.level.singular.value.for.a.given.protein<-function(One_Protein_Name,Data_Frame_With_Singular_Value_Per_Peptide)
{
peptides_for_protein<-which.peptides.belong.to.these.proteins(One_Protein_Name)
Values_For_Peptides_Belonging_to_Protein<-as.numeric(Data_Frame_With_Singular_Value_Per_Peptide[match(peptides_for_protein,row.names(Data_Frame_With_Singular_Value_Per_Peptide)),])
return(Values_For_Peptides_Belonging_to_Protein)
}



#########################This function is identical to the one above with the addition of more flexibility of type of t
calculate_t_scores_for_peptides_and_report_on_protein_level_with_more_options<-function(input_data_frame,Binary_Class,Type_of_T="t",Non_para="n"){
T_values_for_peptides<-as.data.frame(mt.teststat(input_data_frame,Binary_Class,test=Type_of_T,nonpara=Non_para))
row.names(T_values_for_peptides)<-row.names(input_data_frame)
Number_of_Proteins_Represented<-length(unique(which.proteins.do.these.peptides.belong.to(row.names(input_data_frame))))
Data_Frame_for_T_Values_For_Proteins<-as.data.frame((matrix(data=(rep("NA",Number_of_Proteins_Represented*Max_Peptide_Number)),nrow=Number_of_Proteins_Represented,ncol=Max_Peptide_Number,byrow=FALSE)))
Data_Frame_for_T_Values_For_Proteins<-as.data.frame(apply(Data_Frame_for_T_Values_For_Proteins,2,as.numeric))
row.names(Data_Frame_for_T_Values_For_Proteins)<-unique(which.proteins.do.these.peptides.belong.to(row.names(input_data_frame)))
for(i in 1:(dim(Data_Frame_for_T_Values_For_Proteins)[1])){
Protein=row.names(Data_Frame_for_T_Values_For_Proteins[i,])
T_Values<-extract.peptide.level.singular.value.for.a.given.protein(Protein,T_values_for_peptides)
Data_Frame_for_T_Values_For_Proteins[i,1:length(T_Values)]<-T_Values
}


###Order the dataframe such that it is ordered by the row name to maintain consistency with other functions
Data_Frame_for_T_Values_For_Proteins<-Data_Frame_for_T_Values_For_Proteins[order(row.names(Data_Frame_for_T_Values_For_Proteins)),]
return(Data_Frame_for_T_Values_For_Proteins)
}
extract.peptide.level.singular.value.for.a.given.protein<-function(One_Protein_Name,Data_Frame_With_Singular_Value_Per_Peptide)
{
peptides_for_protein<-which.peptides.belong.to.these.proteins(One_Protein_Name)
Values_For_Peptides_Belonging_to_Protein<-as.numeric(Data_Frame_With_Singular_Value_Per_Peptide[match(peptides_for_protein,row.names(Data_Frame_With_Singular_Value_Per_Peptide)),])
return(Values_For_Peptides_Belonging_to_Protein)
}


 #########################This function is identical to the one above with the addition of more flexibility of type of t
 #####Value for trim =0 is mean 1 is median, and 0.5 is 50 percent etc
Calculate.Fold.Change.1.minus.0.for.t.or.pairt.with.trim.and.report.on.protein<-function(input_data_frame,Binary_Class,Type_of_T="t",trimm=0){

if(Type_of_T=="pairt")
{
T_values_for_peptides<-(input_data_frame[,Binary_Class>0.5]-input_data_frame[,Binary_Class<0.05])
T_values_for_peptides<-as.data.frame(apply(T_values_for_peptides,1,mean,trim=trimm))
}

if(Type_of_T=="t")
{
T_values_for_peptides<-as.data.frame(apply(input_data_frame[,Binary_Class>0.5],1,mean,trim=trimm))-(apply(input_data_frame[,Binary_Class<0.5],1,mean,trim=trimm))
}


#T_values_for_peptides<-as.data.frame(mt.teststat(input_data_frame,Binary_Class,test=Type_of_T,nonpara=Non_para))
row.names(T_values_for_peptides)<-row.names(input_data_frame)
Number_of_Proteins_Represented<-length(unique(which.proteins.do.these.peptides.belong.to(row.names(input_data_frame))))
Data_Frame_for_T_Values_For_Proteins<-as.data.frame((matrix(data=(rep("NA",Number_of_Proteins_Represented*Max_Peptide_Number)),nrow=Number_of_Proteins_Represented,ncol=Max_Peptide_Number,byrow=FALSE)))
Data_Frame_for_T_Values_For_Proteins<-as.data.frame(apply(Data_Frame_for_T_Values_For_Proteins,2,as.numeric))
row.names(Data_Frame_for_T_Values_For_Proteins)<-unique(which.proteins.do.these.peptides.belong.to(row.names(input_data_frame)))
for(i in 1:(dim(Data_Frame_for_T_Values_For_Proteins)[1])){
Protein=row.names(Data_Frame_for_T_Values_For_Proteins[i,])
T_Values<-extract.peptide.level.singular.value.for.a.given.protein(Protein,T_values_for_peptides)
Data_Frame_for_T_Values_For_Proteins[i,1:length(T_Values)]<-T_Values
}


###Order the dataframe such that it is ordered by the row name to maintain consistency with other functions
Data_Frame_for_T_Values_For_Proteins<-Data_Frame_for_T_Values_For_Proteins[order(row.names(Data_Frame_for_T_Values_For_Proteins)),]
return(Data_Frame_for_T_Values_For_Proteins)
}














help.me.plot<-function(start_n,start_plus_n,Decreasing_peptide_number_Results.first,Decreasing_peptide_number_T_Scores){

num_methods<-dim(Decreasing_peptide_number_Results.first)[2]-1



stat.graph1<-dotplot(as.matrix(Decreasing_peptide_number_T_Scores[start_n:(start_n+start_plus_n),]),xlim=c(-7,7),col="blue",xlab="T scores for  peptides",auto.key=FALSE,


panel = function(...) {
           panel.grid(h = 0, v = -50)
	   panel.refline(v=c(-2,2),col="yellow")
           panel.refline(v=c(0),col="blue")
           panel.refline(v=c(-3,3),col="dark orange")
           panel.refline(v=c(-4,4),col="red")
           panel.refline(v=c(-5,5),col="dark red")
           panel.dotplot(...)
       }
)
#######################################################################################################################
stat.graph2<-dotplot(as.matrix(Decreasing_peptide_number_Results.first[start_n:(start_n+start_plus_n),1:num_methods],rownames.force=FALSE),xlim=c(0,0.06),xlab="Uncorrected P- Values",
key=simpleKey(pch=1:10,names(Decreasing_peptide_number_Results.first[,1:num_methods]),space="right",cex=0.6)
,pch=1:10,

panel = function(...) {
           panel.grid(h = 0, v = -50)
	   panel.refline(v=0.05,col="black")
          # panel.barchart(...)
           panel.dotplot(...,cex=1.2,xlabels=FALSE,auto.key=FALSE,axis=FALSE,scales=list(draw=FALSE,name=FALSE))
       }

)




#stat.graph2.1<-dotplot(as.matrix(Decreasing_peptide_number_Results.first[start_n:(start_n+start_plus_n),1:num_methods]),xlim=c(0,0.06),xlab="Uncorrected P- Values",key=simpleKey(pch=1:10,names(Decreasing_peptide_number_Results.first[,1:6]),space="right",cex=0.8),pch=1:10,

stat.graph2.1<-dotplot(as.matrix(Decreasing_peptide_number_Results.first[start_n:(start_n+start_plus_n),1:num_methods]),xlim=c(0,0.06),xlab="Uncorrected P- Values",key=simpleKey(pch=1:10,names(Decreasing_peptide_number_Results.first[,1:(dim(Decreasing_peptide_number_Results.first)[2])]),space="right",cex=0.8),pch=1:10,
panel = function(...) {
           panel.grid(h = 0, v = -50)
	   panel.refline(v=0.05,col="black")
          panel.barchart(...)
           #panel.stripplot(...)
           panel.dotplot(...)
       }

)
#######################################################################################################################

plot(stat.graph1,split=c(1,1,2,1))

plot(stat.graph2,split=c(2,1,2,1),newpage=FALSE)
#plot(stat.graph2.1,split=c(2,1,2,1),newpage=FALSE)



}


###Function to convert Simple_Binary_Class to 1 and 0
APT.Simple.Class.Binary.to.Binary.Simple.Class.Binary<-function(Simple,Character_Type_of_Class_Comparison="Two class unpaired")
{
if(Character_Type_of_Class_Comparison=="Two class unpaired")
##This conversion does handle class comparisons in this form (0,0,0,0,1,1,1,1) along with (1,1,1,1,2,2,2,2), infact any two numbers
## see the min function below
{
Binary_Simple_Class_Binary<-as.numeric(Simple>min(Simple))

}
if(Character_Type_of_Class_Comparison=="Two class paired")
{
#This converts the negative numbers to 0,non-negative to 1
Binary_Simple_Class_Binary<-as.numeric(Simple>0)
}
return(Binary_Simple_Class_Binary)
}



#Function1:
# Important ! Paired Analysis should always have the pairs after each other as in: (-1,1,-2,2,-3,3,-4,4) than (-1,-2,-3,-4,1,2,3,4)
#Function that returns all Class Label Permutations given type of class labels as a matrix
# Number of Permutations
# This function handles and has been tested for paired and unpaired:
# Algorithm, 1) get matrix, 2) subset to two rows to increase speed for Samr 3) call samr to conduct hypothesis testing for 20 * number of perm
# 4) extract the [17] object from sam object which indicates the actual class permutations
# Why? Samr does redundant permutations for unparied for some reason, so we call the perm *20 and extract the number of perm unique permutations
# sam excepts class order in 1,1,1,1,2,2,2,2
APT.Return.Class.Permutations.From.Samr<-function(x,Character_Type_of_Class_Comparison,Simple_Binary_Class_Lab,Number_ofPermutations)
{
x<-scramble.data.frame(x)
x<-x[1:2,]
Sam_Obj_from_samr<-samr(list(x=as.matrix(x),y=Simple_Binary_Class_Lab,logged2=TRUE),resp.type=Character_Type_of_Class_Comparison, nperms=(20*Number_ofPermutations),random.seed=1234)
Class_Label_Perms<-unique(as.matrix(as.data.frame(Sam_Obj_from_samr[17])))
correct_order_query<-apply((apply(Class_Label_Perms,1,'!=',Simple_Binary_Class_Lab)),2,sum)
Class_Label_Perms<-Class_Label_Perms[correct_order_query!=0,]
if(Character_Type_of_Class_Comparison=="Two class unpaired")
##This conversion does handle class comparisons in this form (0,0,0,0,1,1,1,1) along with (1,1,1,1,2,2,2,2), infact any two numbers
## see the min function below
{
Class_Label_Perms<-apply(Class_Label_Perms>(min(Simple_Binary_Class_Lab)),2,as.numeric)
}
if(Character_Type_of_Class_Comparison=="Two class paired")
{
#This converts the negative numbers to 0,non-negative to 1
Class_Label_Perms<-apply(Class_Label_Perms>0,2,as.numeric)
}
return(head(Class_Label_Perms,Number_ofPermutations))
}


#Function1.NON UNIQUE FOR SMALL SAMPLE SIZES:
# Important ! Paired Analysis should always have the pairs after each other as in: (-1,1,-2,2,-3,3,-4,4) than (-1,-2,-3,-4,1,2,3,4)
#Function that returns all Class Label Permutations given type of class labels as a matrix
# Number of Permutations
# This function handles and has been tested for paired and unpaired:
# Algorithm, 1) get matrix, 2) subset to two rows to increase speed for Samr 3) call samr to conduct hypothesis testing for 20 * number of perm
# 4) extract the [17] object from sam object which indicates the actual class permutations
# Why? Samr does redundant permutations for unparied for some reason, so we call the perm *20 and extract the number of perm unique permutations
# sam excepts class order in 1,1,1,1,2,2,2,2
APT.Return.Class.Non.Unique.Permutations.From.Samr<-function(x,Character_Type_of_Class_Comparison,Simple_Binary_Class_Lab,Number_ofPermutations)
{
x<-scramble.data.frame(x)
x<-x[1:2,]
Sam_Obj_from_samr<-samr(list(x=as.matrix(x),y=Simple_Binary_Class_Lab,logged2=TRUE),resp.type=Character_Type_of_Class_Comparison, nperms=(20*Number_ofPermutations),random.seed=1234)
Class_Label_Perms<-as.matrix(as.data.frame(Sam_Obj_from_samr[17]))
#correct_order_query<-apply((apply(Class_Label_Perms,1,'!=',Simple_Binary_Class_Lab)),2,sum)
#Class_Label_Perms<-Class_Label_Perms[correct_order_query!=0,]
if(Character_Type_of_Class_Comparison=="Two class unpaired")
##This conversion does handle class comparisons in this form (0,0,0,0,1,1,1,1) along with (1,1,1,1,2,2,2,2), infact any two numbers
## see the min function below
{
Class_Label_Perms<-apply(Class_Label_Perms>(min(Simple_Binary_Class_Lab)),2,as.numeric)
}
if(Character_Type_of_Class_Comparison=="Two class paired")
{
#This converts the negative numbers to 0,non-negative to 1
Class_Label_Perms<-apply(Class_Label_Perms>0,2,as.numeric)
}
return(head(Class_Label_Perms,Number_ofPermutations))
}






##this function handles t statistic generation on peptide level and then report on protein by listing
# all the t scores on the rows for the protein, this way you can create a separate summarization function
##Function 2
APT.Return.Test.Statistic.Protein.Level.t.or.paired.t.parameteric.and.non.parametric<-function(input_data_frame,Character_Type_of_Class_Comparison,Simple_Binary_Class_,parametric_y_or_n){
if(Character_Type_of_Class_Comparison=="Two class unpaired")
{
Character_Type_of_Class_Comparison<-"t"
}
if(Character_Type_of_Class_Comparison=="Two class paired")
{
#This converts the negative numbers to 0,non-negative to 1
Character_Type_of_Class_Comparison<-"pairt"
}
T_values_for_peptides<-as.data.frame(mt.teststat(input_data_frame,Simple_Binary_Class_,test=Character_Type_of_Class_Comparison,nonpara=parametric_y_or_n))
row.names(T_values_for_peptides)<-row.names(input_data_frame)
Number_of_Proteins_Represented<-length(unique(which.proteins.do.these.peptides.belong.to(row.names(input_data_frame))))
Data_Frame_for_T_Values_For_Proteins<-as.data.frame((matrix(data=(rep("NA",Number_of_Proteins_Represented*Max_Peptide_Number)),nrow=Number_of_Proteins_Represented,ncol=Max_Peptide_Number,byrow=FALSE)))
Data_Frame_for_T_Values_For_Proteins<-as.data.frame(apply(Data_Frame_for_T_Values_For_Proteins,2,as.numeric))
row.names(Data_Frame_for_T_Values_For_Proteins)<-unique(which.proteins.do.these.peptides.belong.to(row.names(input_data_frame)))
for(i in 1:(dim(Data_Frame_for_T_Values_For_Proteins)[1])){
Protein=row.names(Data_Frame_for_T_Values_For_Proteins[i,])
T_Values<-extract.peptide.level.singular.value.for.a.given.protein(Protein,T_values_for_peptides)
Data_Frame_for_T_Values_For_Proteins[i,1:length(T_Values)]<-T_Values
}

###Order the dataframe such that it is ordered by the row name to maintain consistency with other functions
Data_Frame_for_T_Values_For_Proteins<-Data_Frame_for_T_Values_For_Proteins[order(row.names(Data_Frame_for_T_Values_For_Proteins)),]
return(Data_Frame_for_T_Values_For_Proteins)
}
extract.peptide.level.singular.value.for.a.given.protein<-function(One_Protein_Name,Data_Frame_With_Singular_Value_Per_Peptide)
{
peptides_for_protein<-which.peptides.belong.to.these.proteins(One_Protein_Name)
Values_For_Peptides_Belonging_to_Protein<-as.numeric(Data_Frame_With_Singular_Value_Per_Peptide[match(peptides_for_protein,row.names(Data_Frame_With_Singular_Value_Per_Peptide)),])
return(Values_For_Peptides_Belonging_to_Protein)
}


#######This function summarized Peptide level Test Statistic to the Protein level 
# Take the output of Data   APT.Return.Test.Statistic.Protein.Level.t.or.paired.t.parameteric.and.non.parametric
# which already includes protein level statistics but un summarized, the input data.frame kinda looks like this
#  Protein1 -1,-2.34,-2,NA,NA
# The dimensions of the data frame is the number of proteins and columns is number of max peptides for that protein data frame
# hence the NAs
# This function outputs a dataframe with the same number of rows as the input but only one column since it summarizes
# This function can summarize in many ways, the default it sum, this method is optimal for peptide to protein since
# you can to control for the random variability around 0 t score because of probability
# Max mean would be tibshirani GSA method, assuming the test scores are Z-Values, something this function does have to care about
# GSEA would not be be implementable if the R correlation values are first ranked and then a Smirnov Fologomogorov like test statistics is calculated
# Remember there is also a t score option available in GSEA windows edition, so try to compare the results from here and from the program
# There is a normalization step after the summarization is conducted.

APT.Summarize.Peptide.Level.Statistic.To.Protein<-function(OutPut_From_APT.Return.Test.Statistic_etc,Type="sum")
{
if(Type=="sum"){
OutPut_From_APT.Return.Test.Statistic_etc_new<-apply(OutPut_From_APT.Return.Test.Statistic_etc,1,sum,na.rm=TRUE)
}

if(Type=="mean"){
OutPut_From_APT.Return.Test.Statistic_etc_new<-apply(OutPut_From_APT.Return.Test.Statistic_etc,1,mean,na.rm=TRUE)
}

if(Type=="median"){
OutPut_From_APT.Return.Test.Statistic_etc_new<-apply(OutPut_From_APT.Return.Test.Statistic_etc,1,median,na.rm=TRUE)
}


OutPut_From_APT.Return.Test.Statistic_etc_new<-data.frame(Summarized_Test_Statistic=as.numeric(OutPut_From_APT.Return.Test.Statistic_etc_new))
row.names(OutPut_From_APT.Return.Test.Statistic_etc_new)<-row.names(OutPut_From_APT.Return.Test.Statistic_etc)
return(OutPut_From_APT.Return.Test.Statistic_etc_new)
}


### Normalization Function for Protein level test statistics, created to mimic GSEA, Sigpathway, GSA Algorithm
### As of now, this function does not normalize the data in any way, but is a good practise to include in the workflow
#for future implementation
APT.Normalize.Protein.Level.Summarized.Test.Statistic<-function(APT.Summarize.Peptide.Level.Statistic.To.Protein_)
{

return(APT.Summarize.Peptide.Level.Statistic.To.Protein_)
}


####################

APT.Calculate.Master<-function(Input_data_frame_peptide_level,Simple_Binary_Class_,Type_of_Comparison,Number_of_Permutations,parametric_or_not,Type_of_summarization)
{
#Binary_Converted_Dynamic<-Binary_Simple_Class_Binary<-APT.Simple.Class.Binary.to.Binary.Simple.Class.Binary(Simple_Binary_Class_)

Binary_Converted_Dynamic<-APT.Simple.Class.Binary.to.Binary.Simple.Class.Binary(Simple_Binary_Class_,Type_of_Comparison)

Original_T_Scores_Without_Perms<-(APT.Normalize.Protein.Level.Summarized.Test.Statistic(APT.Summarize.Peptide.Level.Statistic.To.Protein((APT.Return.Test.Statistic.Protein.Level.t.or.paired.t.parameteric.and.non.parametric(Input_data_frame_peptide_level,Type_of_Comparison,Binary_Converted_Dynamic,parametric_or_not)),Type=Type_of_summarization)))


Original_T_Scores_With_Perms<-Original_T_Scores_Without_Perms



Matrix_with_Permuted_Class_Labels<-APT.Return.Class.Permutations.From.Samr(Input_data_frame_peptide_level,Type_of_Comparison,Simple_Binary_Class_,Number_of_Permutations)


for(i in 1:(dim(Matrix_with_Permuted_Class_Labels)[1]) ){
Original_T_Scores_With_Perms[i]<-(APT.Normalize.Protein.Level.Summarized.Test.Statistic(APT.Summarize.Peptide.Level.Statistic.To.Protein((APT.Return.Test.Statistic.Protein.Level.t.or.paired.t.parameteric.and.non.parametric(Input_data_frame_peptide_level,Type_of_Comparison,Matrix_with_Permuted_Class_Labels[i,],parametric_or_not)),Type=Type_of_summarization)))
print(i)
}


return(cbind(Original_T_Scores_With_Perms,Original_T_Scores_Without_Perms))

for(i in 1:dim(er)[1]){
pvalue[i,]<-((1+dim(er)[2])-((order(cbind(abs(yes_orig[i,]),abs(yes[i,])))[1])))/(1+dim(yes)[2])
}


}


APT.Calculate.Master.with.Redundant.Sample.Perms<-function(Input_data_frame_peptide_level,Simple_Binary_Class_,Type_of_Comparison,Number_of_Permutations,parametric_or_not,Type_of_summarization)
{
#Binary_Converted_Dynamic<-Binary_Simple_Class_Binary<-APT.Simple.Class.Binary.to.Binary.Simple.Class.Binary(Simple_Binary_Class_)

Binary_Converted_Dynamic<-APT.Simple.Class.Binary.to.Binary.Simple.Class.Binary(Simple_Binary_Class_,Type_of_Comparison)

Original_T_Scores_Without_Perms<-(APT.Normalize.Protein.Level.Summarized.Test.Statistic(APT.Summarize.Peptide.Level.Statistic.To.Protein((APT.Return.Test.Statistic.Protein.Level.t.or.paired.t.parameteric.and.non.parametric(Input_data_frame_peptide_level,Type_of_Comparison,Binary_Converted_Dynamic,parametric_or_not)),Type=Type_of_summarization)))


Original_T_Scores_With_Perms<-Original_T_Scores_Without_Perms


##### THis is only change from the original
Matrix_with_Permuted_Class_Labels<-APT.Return.Class.Non.Unique.Permutations.From.Samr(Input_data_frame_peptide_level,Type_of_Comparison,Simple_Binary_Class_,Number_of_Permutations)


for(i in 1:(dim(Matrix_with_Permuted_Class_Labels)[1]) ){
Original_T_Scores_With_Perms[i]<-(APT.Normalize.Protein.Level.Summarized.Test.Statistic(APT.Summarize.Peptide.Level.Statistic.To.Protein((APT.Return.Test.Statistic.Protein.Level.t.or.paired.t.parameteric.and.non.parametric(Input_data_frame_peptide_level,Type_of_Comparison,Matrix_with_Permuted_Class_Labels[i,],parametric_or_not)),Type=Type_of_summarization)))
print(i)
}


return(cbind(Original_T_Scores_With_Perms,Original_T_Scores_Without_Perms))

for(i in 1:dim(er)[1]){
pvalue[i,]<-((1+dim(er)[2])-((order(cbind(abs(yes_orig[i,]),abs(yes[i,])))[1])))/(1+dim(yes)[2])
}


}


APT.Calculate.Master.with.Simultaneous.Peptide.Sample.Perms<-function(Input_data_frame_peptide_level,Simple_Binary_Class_,Type_of_Comparison,Number_of_Permutations,parametric_or_not,Type_of_summarization,normalize_after_each_scramble=FALSE)
{
#Binary_Converted_Dynamic<-Binary_Simple_Class_Binary<-APT.Simple.Class.Binary.to.Binary.Simple.Class.Binary(Simple_Binary_Class_)
Binary_Converted_Dynamic<-APT.Simple.Class.Binary.to.Binary.Simple.Class.Binary(Simple_Binary_Class_,Type_of_Comparison)
Original_T_Scores_Without_Perms<-(APT.Normalize.Protein.Level.Summarized.Test.Statistic(APT.Summarize.Peptide.Level.Statistic.To.Protein((APT.Return.Test.Statistic.Protein.Level.t.or.paired.t.parameteric.and.non.parametric(Input_data_frame_peptide_level,Type_of_Comparison,Binary_Converted_Dynamic,parametric_or_not)),Type=Type_of_summarization)))
Original_T_Scores_With_Perms<-Original_T_Scores_Without_Perms
##### THis is only change from the original
#Matrix_with_Permuted_Class_Labels<-APT.Return.Class.Non.Unique.Permutations.From.Samr(Input_data_frame_peptide_level,Type_of_Comparison,Simple_Binary_Class_,Number_of_Permutations)
data_frame_with_original_values<-Input_data_frame_peptide_level
for(i in 1:Number_of_Permutations ){
Input_data_frame_peptide_level<-scramble.data.frame(data_frame_with_original_values)
if(normalize_after_each_scramble==TRUE){
Input_data_frame_peptide_level<-as.data.frame(normalizeBetweenArrays(as.matrix(Input_data_frame_peptide_level),method="scale"))
}
Original_T_Scores_With_Perms[i]<-(APT.Normalize.Protein.Level.Summarized.Test.Statistic(APT.Summarize.Peptide.Level.Statistic.To.Protein((APT.Return.Test.Statistic.Protein.Level.t.or.paired.t.parameteric.and.non.parametric(Input_data_frame_peptide_level,Type_of_Comparison,Binary_Converted_Dynamic,parametric_or_not)),Type=Type_of_summarization)))
print(i)
}
return(cbind(Original_T_Scores_With_Perms,Original_T_Scores_Without_Perms))
}




APT.Calculate.Master.with.redundant.Peptide.Perms<-function(Input_data_frame_peptide_level,Simple_Binary_Class_,Type_of_Comparison,Number_of_Permutations,parametric_or_not,Type_of_summarization)
{
#Binary_Converted_Dynamic<-Binary_Simple_Class_Binary<-APT.Simple.Class.Binary.to.Binary.Simple.Class.Binary(Simple_Binary_Class_)

Binary_Converted_Dynamic<-APT.Simple.Class.Binary.to.Binary.Simple.Class.Binary(Simple_Binary_Class_,Type_of_Comparison)

Original_T_Scores_Without_Perms<-(APT.Normalize.Protein.Level.Summarized.Test.Statistic(APT.Summarize.Peptide.Level.Statistic.To.Protein((APT.Return.Test.Statistic.Protein.Level.t.or.paired.t.parameteric.and.non.parametric(Input_data_frame_peptide_level,Type_of_Comparison,Binary_Converted_Dynamic,parametric_or_not)),Type=Type_of_summarization)))


Original_T_Scores_With_Perms<-Original_T_Scores_Without_Perms


##### THis is only change from the original
#Matrix_with_Permuted_Class_Labels<-APT.Return.Class.Non.Unique.Permutations.From.Samr(Input_data_frame_peptide_level,Type_of_Comparison,Simple_Binary_Class_,Number_of_Permutations)


for(i in 1:Number_of_Permutations ){
Input_data_frame_peptide_level<-scramble.data.frame(Input_data_frame_peptide_level)
Original_T_Scores_With_Perms[i]<-(APT.Normalize.Protein.Level.Summarized.Test.Statistic(APT.Summarize.Peptide.Level.Statistic.To.Protein((APT.Return.Test.Statistic.Protein.Level.t.or.paired.t.parameteric.and.non.parametric(Input_data_frame_peptide_level,Type_of_Comparison,Binary_Converted_Dynamic,parametric_or_not)),Type=Type_of_summarization)))
print(i)
}


return(cbind(Original_T_Scores_With_Perms,Original_T_Scores_Without_Perms))

for(i in 1:dim(er)[1]){
pvalue[i,]<-((1+dim(er)[2])-((order(cbind(abs(yes_orig[i,]),abs(yes[i,])))[1])))/(1+dim(yes)[2])
}


}






##this function handles t statistic generation on peptide level and then report on protein by listing
# all the t scores on the rows for the protein, this way you can create a separate summarization function
##Function 2
APT.Return.Test.Statistic.Protein.Level.t.or.paired.t.parameteric.and.non.parametric<-function(input_data_frame,Character_Type_of_Class_Comparison,Simple_Binary_Class_,parametric_y_or_n){
if(Character_Type_of_Class_Comparison=="Two class unpaired")
{
Character_Type_of_Class_Comparison<-"t"
}
if(Character_Type_of_Class_Comparison=="Two class paired")
{
#This converts the negative numbers to 0,non-negative to 1
Character_Type_of_Class_Comparison<-"pairt"
}
T_values_for_peptides<-as.data.frame(mt.teststat(input_data_frame,Simple_Binary_Class_,test=Character_Type_of_Class_Comparison,nonpara=parametric_y_or_n))
row.names(T_values_for_peptides)<-row.names(input_data_frame)
Number_of_Proteins_Represented<-length(unique(which.proteins.do.these.peptides.belong.to(row.names(input_data_frame))))
Data_Frame_for_T_Values_For_Proteins<-as.data.frame((matrix(data=(rep("NA",Number_of_Proteins_Represented*Max_Peptide_Number)),nrow=Number_of_Proteins_Represented,ncol=Max_Peptide_Number,byrow=FALSE)))
Data_Frame_for_T_Values_For_Proteins<-as.data.frame(apply(Data_Frame_for_T_Values_For_Proteins,2,as.numeric))
row.names(Data_Frame_for_T_Values_For_Proteins)<-unique(which.proteins.do.these.peptides.belong.to(row.names(input_data_frame)))
for(i in 1:(dim(Data_Frame_for_T_Values_For_Proteins)[1])){
Protein=row.names(Data_Frame_for_T_Values_For_Proteins[i,])
T_Values<-extract.peptide.level.singular.value.for.a.given.protein(Protein,T_values_for_peptides)
Data_Frame_for_T_Values_For_Proteins[i,1:length(T_Values)]<-T_Values
}

###Order the dataframe such that it is ordered by the row name to maintain consistency with other functions
Data_Frame_for_T_Values_For_Proteins<-Data_Frame_for_T_Values_For_Proteins[order(row.names(Data_Frame_for_T_Values_For_Proteins)),]
return(Data_Frame_for_T_Values_For_Proteins)
}
extract.peptide.level.singular.value.for.a.given.protein<-function(One_Protein_Name,Data_Frame_With_Singular_Value_Per_Peptide)
{
peptides_for_protein<-which.peptides.belong.to.these.proteins(One_Protein_Name)
Values_For_Peptides_Belonging_to_Protein<-as.numeric(Data_Frame_With_Singular_Value_Per_Peptide[match(peptides_for_protein,row.names(Data_Frame_With_Singular_Value_Per_Peptide)),])
return(Values_For_Peptides_Belonging_to_Protein)
}



##########################################################################################				
##########################################################################################
#####################Annotation and Systems Biology#######################################				
##########################################################################################
##########################################################################################

Fetch.and.add.All.Basic.Protein.Annotation.for.these.IPI.proteins.row.named.data.frame.microsoft.windows<-function(Input_Data_Frame){
library("org.Hs.ipi")
APPPT_Protein_Annotation<-as.data.frame(Input_Data_Frame[,1])
row.names(APPPT_Protein_Annotation)<-row.names(Input_Data_Frame)
APPPT_Protein_Annotation[,1]<-gsub("_.*","",row.names(Input_Data_Frame))
names(APPPT_Protein_Annotation)[1]<-"IPI_ID"
APPPT_Protein_Annotation[,2]<-match(APPPT_Protein_Annotation[,1],  gsub("[.].*","",names(as.list(org.Hs.ipiGENEID)))  )
names(APPPT_Protein_Annotation)[2]<-"org.Hs.ipi.match.Index"
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(Gene_Alias=as.character(as.list(org.Hs.ipiSYMBOL))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(Description=as.character(as.list(org.Hs.ipiDE))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(Chromosome=as.character(as.list(org.Hs.ipiCHR))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(MW=as.character(as.list(org.Hs.ipiMW))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(Length=as.character(as.list(org.Hs.ipiLEN))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(Orientation=as.character(as.list(org.Hs.ipiORIENT))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(START=as.character(as.list(org.Hs.ipiSTART))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(END=as.character(as.list(org.Hs.ipiEND))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(Sequence=as.character(as.list(org.Hs.ipiSQ))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(UNIGENE=as.character(as.list(org.Hs.ipiUNIGENE))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(GI=as.character(as.list(org.Hs.ipiGI))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(RefSeq=as.character(as.list(org.Hs.ipiREFSEQ))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
return(APPPT_Protein_Annotation)
}

##############################################################################
###Function to Return same as above, but designed to take peptide level intensity data.frame with samples separated, i.e. Real_Processed_Data_Frame, type. It returns a data frame with Annotations, and also Protein Intensity Ranks, and # of Peptides per Protein############
###############################################
Fetch.peptide.stats.and.All.Basic.Protein.Annotation.for.these.IPI.peptides.row.named.data.frame.microsoft.windows<-function(peptide__intensity_Input_Data_Frame){
###This line below first summarizes peptide intensity to protein by averaging the intesities, then 
### with the protein level intensity table, it calculates the maximum intensity each protein among all the samples, then
### it converts the protein level intensities to it's ranks, converts it to a percentile, by dividing it by the max rank * 100
#Input_Data_Frame<-cbind(100*(rank(apply(summarize.peptide.data.frame.on.protein.by.averaging(peptide__intensity_Input_Data_Frame),1,max))/(max(rank(apply(summarize.peptide.data.frame.on.protein.by.averaging(peptide__intensity_Input_Data_Frame),1,max))))),how.many.peptides.belong.to.each.of.these.proteins(row.names(summarize.peptide.data.frame.on.protein.by.averaging(peptide__intensity_Input_Data_Frame))))
#################
####Instead of above complication, summarize peptide intensity by averaging to the protein level, and then get the max of the protein intensity acroass the samples
Input_Data_Frame<-cbind(apply(summarize.peptide.data.frame.on.protein.by.averaging(peptide__intensity_Input_Data_Frame),1,max),how.many.peptides.belong.to.each.of.these.proteins(row.names(summarize.peptide.data.frame.on.protein.by.averaging(peptide__intensity_Input_Data_Frame))))
Input_Data_Frame<-as.data.frame(Input_Data_Frame)
names(Input_Data_Frame)<-c("Prot_Ints","#Peptides")
library("org.Hs.ipi")
APPPT_Protein_Annotation<-as.data.frame(Input_Data_Frame[,1])
row.names(APPPT_Protein_Annotation)<-row.names(Input_Data_Frame)
APPPT_Protein_Annotation[,1]<-gsub("_.*","",row.names(Input_Data_Frame))
names(APPPT_Protein_Annotation)[1]<-"IPI_ID"
APPPT_Protein_Annotation[,2]<-match(APPPT_Protein_Annotation[,1],  gsub("[.].*","",names(as.list(org.Hs.ipiGENEID)))  )
names(APPPT_Protein_Annotation)[2]<-"org.Hs.ipi.match.Index"
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(Gene_Alias=as.character(as.list(org.Hs.ipiSYMBOL))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(Description=as.character(as.list(org.Hs.ipiDE))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(Chromosome=as.character(as.list(org.Hs.ipiCHR))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(MW=as.numeric(as.character(as.list(org.Hs.ipiMW))[c(as.numeric(APPPT_Protein_Annotation[,2]))])))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(Length=as.numeric(as.character(as.list(org.Hs.ipiLEN))[c(as.numeric(APPPT_Protein_Annotation[,2]))])))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(Orientation=as.character(as.list(org.Hs.ipiORIENT))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(START=as.numeric(as.character(as.list(org.Hs.ipiSTART))[c(as.numeric(APPPT_Protein_Annotation[,2]))])))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(END=as.numeric(as.character(as.list(org.Hs.ipiEND))[c(as.numeric(APPPT_Protein_Annotation[,2]))])))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(Sequence=as.character(as.list(org.Hs.ipiSQ))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(UNIGENE=as.character(as.list(org.Hs.ipiUNIGENE))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(GI=as.character(as.list(org.Hs.ipiGI))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
APPPT_Protein_Annotation<-cbind(APPPT_Protein_Annotation,data.frame(RefSeq=as.character(as.list(org.Hs.ipiREFSEQ))[c(as.numeric(APPPT_Protein_Annotation[,2]))]))
return(cbind(Input_Data_Frame,APPPT_Protein_Annotation))
}
#################################################################
################################################################



############################################3
#################
Convert.NCBI.Interaction.csv.file.based.on.column.Human_Gene_Alias<-function(Name.of.csv.file){
NCBI_Hiv<-read.csv(Name.of.csv.file,header=TRUE)
NCBI_Hiv<-NCBI_Hiv[!is.na(NCBI_Hiv$Human_Gene_Alias),]
NCBI_Hiv<-NCBI_Hiv[!(NCBI_Hiv$Human_Gene_Alias=="NA"),]
NCBI_HIV_INTERACTIONS<-as.data.frame(rep(0,length(unique(NCBI_Hiv$Human_Gene_Alias))),row.names=as.character(unique(NCBI_Hiv$Human_Gene_Alias)))

#row.names(NCBI_HIV_INTERACTIONS)<-as.character(NCBI_HIV_INTERACTIONS[,1])
for(i in 1:dim(NCBI_HIV_INTERACTIONS)[1]){
NCBI_HIV_INTERACTIONS[i,]<-sum(NCBI_Hiv$Human_Gene_Alias==row.names(NCBI_HIV_INTERACTIONS)[i])
}
NCBI_HIV_INTERACTIONS<-cbind(NCBI_HIV_INTERACTIONS,((dim(NCBI_HIV_INTERACTIONS)[1]+1)-rank(NCBI_HIV_INTERACTIONS[,1],ties.method="first")))
names(NCBI_HIV_INTERACTIONS)<-c("Cited HIV Interactions","Rank")
#Index<-match(NCBI_Hiv$Human_Gene_Alias,row.names(NCBI_HIV_INTERACTIONS))
#Genes<-row.names(NCBI_HIV_INTERACTIONS)
List_of_Genes_and_Ncbi_File<-list(row.names(NCBI_HIV_INTERACTIONS),NCBI_HIV_INTERACTIONS,NCBI_Hiv)
names(List_of_Genes_and_Ncbi_File)<-c("Genes","Basic","Details")
return(List_of_Genes_and_Ncbi_File)
}
#####################


######################################################Advanced Function##########################

#############################Clustering HOPACH####################


conduct.hopach.clustering.analysis.on.samples.for.peptide.or.protein.level.data<-function(input_data_frame_peptide_or_protein,distance_of_choice="cor",clusterdirection="sample",perms=1000){
# This function conducts a significance analysis on clusters created via bootstrapping, this function
# clusters samples by default because of the time required
# The currently available options are "cosangle" (cosine angle or uncentered correlation distance),
# "abscosangle" (absolute cosine angle or absolute uncentered correlation distance), "euclid" (Euclidean distance), 
#"abseuclid" (absolute Euclidean distance), "cor" (correlation distance), and "abscor" (absolute correlation distance)

if(clusterdirection=="sample"){
input_data_frame_peptide_or_protein=as.data.frame(t(input_data_frame_peptide_or_protein))
}

temp_matrix_to_store_pair_wise_dist<-distancematrix(as.matrix(input_data_frame_peptide_or_protein),distance_of_choice,na.rm=TRUE)
print("Performing Hopach...")
temp_hopach.output<-hopach(as.matrix(input_data_frame_peptide_or_protein),dmat=temp_matrix_to_store_pair_wise_dist)
#dist.hiv.hopach$clust$k
#table(dist.hiv.hopach$clust$labels)
#table(dist.hiv.hopach$clust$sizes)
dplot(temp_matrix_to_store_pair_wise_dist, temp_hopach.output, ord = "final", main = paste("Clustering from matrix M: #Clusters = ",temp_hopach.output$clust$k,sep=""))

temp_hopach_bobj <- boothopach(as.matrix(input_data_frame_peptide_or_protein),temp_hopach.output, B = perms)
print("Performing Bootstraps...")
bootplot(temp_hopach_bobj,temp_hopach.output,main=paste("Bootstrap Fuzzy Clustering. B = ",perms,sep=""))

makeoutput(input_data_frame_peptide_or_protein, temp_hopach.output, temp_hopach_bobj, file = paste(standardized_output_file_name_apppt,"_HOPACH_Output.txt",sep=""),gene.names = row.names(input_data_frame_peptide_or_protein))
print(paste("Writing results to file: ",standardized_output_file_name_apppt,"_HOPACH_Output.txt",sep=""))

#hopach2tree(input_data_frame_peptide_or_protein, file = "GolubTree", hopach.genes = temp_hopach.output,hopach.arrays = NULL, dist.genes = temp_matrix_to_store_pair_wise_dist, gene.names = row.names(input_data_frame_peptide_or_protein))
#Output for Michael Einsen's Treeview
}


convert.peptide.sequence.to.mass_avg<-function(peptide,calculate_begin_end=TRUE,calculate_modification_addition=TRUE)
{
# Calculator Details
#http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator-notes.asp
modification_addition<-0
if(calculate_modification_addition==TRUE){
modification_addition<-sum(as.numeric(strsplit(gsub("\\[","",gsub("[A-Z]","", peptide,perl=TRUE),perl=TRUE),split="\\]")[[1]]))
}

#gsub("\\[","",gsub("[A-Z]","", asdff,perl=TRUE),perl=TRUE)
peptide<-gsub("\\[.+\\]","", peptide,perl=TRUE)
##Remove the Brackets and numbers within for rest of program
temp_mass<-strsplit(peptide,split="")[[1]]
temp_mass<-gsub("A",71.0788,temp_mass)


temp_mass<-gsub("N",114.1039,temp_mass)

temp_mass<-gsub("C",103.1448,temp_mass)

temp_mass<-gsub("Q",128.1308,temp_mass)

temp_mass<-gsub("H",137.1412,temp_mass)

temp_mass<-gsub("L",113.1595,temp_mass)

temp_mass<-gsub("M",131.1986,temp_mass)

temp_mass<-gsub("P",97.1167,temp_mass)

temp_mass<-gsub("T",101.1051,temp_mass)

temp_mass<-gsub("Y",163.1760,temp_mass)

temp_mass<-gsub("R",156.1876,temp_mass)

temp_mass<-gsub("D",115.0886,temp_mass)

temp_mass<-gsub("E",129.1155,temp_mass)

temp_mass<-gsub("G",57.052,temp_mass)

temp_mass<-gsub("I",113.1595,temp_mass)

temp_mass<-gsub("K",128.1742,temp_mass)

temp_mass<-gsub("N",114.1039,temp_mass)

temp_mass<-gsub("S",87.0782,temp_mass)

temp_mass<-gsub("W",186.2133,temp_mass)

temp_mass<-gsub("V",99.1326,temp_mass)


temp_mass<-gsub("F",147.1766,temp_mass)

if(calculate_begin_end==FALSE){
total<-sum(as.numeric(temp_mass))
}
if(calculate_begin_end==TRUE){
total<-sum(as.numeric(temp_mass)) + 16.00026 + (1.0079*2) 

}
total<-total+modification_addition
return(total)
}
##############################################################################################
#############################################################################################3
##############################################################################################3
Master.Convert.Peptide.String.OR.Data.Frame.With.Row.Name<-function(data.frame.input.with.row.names.or.string.with.peptides,data_frame_toggle=TRUE,calculate_begin_end=TRUE,calculate_modification_addition=TRUE)
{
if(data_frame_toggle==TRUE)
{
Peptide_Masses<-data.frame(Peptide_MW=rep(0,(dim(data.frame.input.with.row.names.or.string.with.peptides)[1])))  
row.names(Peptide_Masses)<-row.names(data.frame.input.with.row.names.or.string.with.peptides)
  for(i in 1:dim(data.frame.input.with.row.names.or.string.with.peptides)[1])
  {
              Peptide_Masses[i,]<-convert.peptide.sequence.to.mass_avg( (row.names(data.frame.input.with.row.names.or.string.with.peptides)[i])    ,calculate_begin_end,calculate_modification_addition)
   }
}

if(data_frame_toggle==FALSE)
{
}

return(Peptide_Masses)
}


##############################################################################################
#############################################################################################3
##############################################################################################3
 
Master.Calculate.Peptide.Masses.Peptide.Mean.Protein<-function(Input_data_frame_with_peptide_row_names,calculate_begin_end=TRUE,calculate_modification_addition=TRUE)
{
Resulting_Data_Frame<-Master.Convert.Peptide.String.OR.Data.Frame.With.Row.Name(Input_data_frame_with_peptide_row_names,data_frame_toggle=TRUE,calculate_begin_end=TRUE,calculate_modification_addition=TRUE)
Resulting_Data_Frame<-cbind(Resulting_Data_Frame,Redundant_Proteins,as.data.frame(apply(gct_file.dataframe,1,mean,na.rm=TRUE)))
names(Resulting_Data_Frame)[2:3]<-c("Protein","Average_Peptide_Abundance")
return(Resulting_Data_Frame)

}
##############################################################################################
#############################################################################################3
##############################################################################################3

convert.data.frame.to.lattice<-function(input_data_frame_)
{
  Resulting_lattice_format<-data.frame(Orig_Col_Names=rep(0,(dim(input_data_frame_)[1]*dim(input_data_frame_)[2])), Data=as.numeric(as.matrix(input_data_frame_)), Orig_Row_Names=
  rep(row.names(input_data_frame_),(dim(input_data_frame_)[2]))

    )

for(i in 1:dim(input_data_frame_)[2]){
 #print(i)
 #print( 1+((dim(input_data_frame_)[1])*(i-1)))
 #print(  dim(input_data_frame_)[1] *i)
   Resulting_lattice_format[(1+((dim(input_data_frame_)[1])*(i-1))): ((dim(input_data_frame_)[1]) *i),1]<- rep(names(input_data_frame_)[i],(dim(input_data_frame_)[1]))
}

return(Resulting_lattice_format)
}
