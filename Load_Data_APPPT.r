# ste("./Parsed_Files_For_Bioconductor/",list.files("./Parsed_Files_For_Bioconductor/")[2],sep="")
#########################################################################################

                                        #Important Variables
##########################################################################################
#                                 # Assign Variables and Load Data
# The command looks into the below mentioned directory and expects only three files namely...
# Constructs a name for the pdf file by concatinating it with the date and word APPPT

#####################################################Naming the output files based on input files
standardized_output_file_name_apppt<-paste(
strsplit((list.files("./Parsed_Files_For_Bioconductor/")[1]),split=".gct.bioconductor")[[1]][1],
"APPPT",
strsplit(date(),split=" ")[[1]][2],
strsplit(date(),split=" ")[[1]][3],
strsplit(date(),split=" ")[[1]][5],
#strsplit(date(),split=" ")[[1]][4],    # Minutes and Seconds
sep="_")
####################################################

gct_file<-paste("./Parsed_Files_For_Bioconductor/",list.files("./Parsed_Files_For_Bioconductor/")[1],sep="")
gmt_file<-paste("./Parsed_Files_For_Bioconductor/",list.files("./Parsed_Files_For_Bioconductor/")[2],sep="")
redundant_proteins_file_name<-paste("./Parsed_Files_For_Bioconductor/",list.files("./Parsed_Files_For_Bioconductor/")[3],sep="")


gct_file<-paste("./Parsed_Files_For_Bioconductor/",(list.files("./Parsed_Files_For_Bioconductor/")[ (grepl(("[.]gct[.]bioconductor$"),  list.files("./Parsed_Files_For_Bioconductor/"))) ]), sep="")
gmt_file<-paste("./Parsed_Files_For_Bioconductor/",(list.files("./Parsed_Files_For_Bioconductor/")[ (grepl(("[.]gmt$"),  list.files("./Parsed_Files_For_Bioconductor/"))) ]), sep="")
redundant_proteins_file_name<-paste("./Parsed_Files_For_Bioconductor/",(list.files("./Parsed_Files_For_Bioconductor/")[grepl(("_redundant[.]Proteins$"),  list.files("./Parsed_Files_For_Bioconductor/"))]), sep="")

#list.files("./Parsed_Files_For_Bioconductor/")[  (grepl((".gct."),  list.files("./Parsed_Files_For_Bioconductor/")) )], sep="")

#########################################################################################
				  #Global Variables used by everyone
gct_file.dataframe<-read.table(gct_file,header=TRUE,na.string=0)   #Read the data.frame in
gct_file.matrix<-as.matrix(gct_file.dataframe)         #Convert to Matrix for GSA program
Peptides<-row.names(gct_file.dataframe)           #Extract Peptide.names from row
Number_of_Peptides<-length(Peptides)              #Calculate number of peptides
Proteins<-GSA.read.gmt(gmt_file)                  #Get the Protein Sets Using GSA function, this file is specific to the package GSA
Number_of_Proteins<-length(Proteins$genesets)     #Calculate the Number of Proteins
Redundant_Proteins<-read.table(redundant_proteins_file_name,header=FALSE) #To View Peptides from Same Protein
Simple_Binary_Class<-extract.class.designatio.of.binary.designation(gct_file.dataframe)
Binary_Simple_Class_Binary<-APT.Simple.Class.Binary.to.Binary.Simple.Class.Binary(Simple_Binary_Class)
Character_Class_Designation<-(as.character(as.matrix(as.data.frame(strsplit(names(gct_file.dataframe),split="_"))[1,])))
Character_Type_of_Class_Comparison="Two class unpaired"
Proteins_and_Peptide_Numbers<-extract.protein.name.and.number.of.peptides(Proteins)
Max_Peptide_Number<-max(Proteins_and_Peptide_Numbers)

############################################################################################
				#Variable Used by Specific Packages
sig_pathway_gmt_load_file<-(gmtToG(gmt_file, verbose = TRUE))