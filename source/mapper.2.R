#! /usr/bin/Rscript
##! /usr/lib64/R/bin/Rscript
##! /usr/bin/env Rscript 


#########    /usr/lib64/R/bin/Rscript #ATLAS
#########    /usr/bin/env Rscript  #local
#########    /usr/bin/Rscript #azure

#indeces that need to be counted given the run specifics
data.type.position <- 8# the data type (train, test or srd)
value.positions <- c(3:7,9:13)
REP_N <- as.numeric(Sys.getenv("REP_N"))


input <- file( "stdin" , "r" )


while( TRUE ){
	 currentLine <- readLines( input , n=1 ) 
	 if( 0 == length( currentLine ) ){
	 	break
	 }
	 currentFields <- unlist( strsplit( currentLine , "," ) ) 
	 ## SRD needs date in the Key!
	 if(currentFields[8]=="srd"){
		 key <- paste(currentFields[data.type.position],currentFields[5],currentFields[2],sep="-") # <type> - < date> - row.id
	 }else{
		 key <- paste(currentFields[data.type.position],currentFields[2],sep="-") # <type> - row.id
	 }

	 fold.id <- as.numeric(strsplit(currentFields[1],"-")[[1]][1])
	 sample.id <- ceiling(fold.id/REP_N)

	 value <- paste(currentFields[value.positions],collapse=",")
	 value <- paste(as.character(sample.id),value,sep=",")

	 result <- paste(key,value,sep="\t")
	 cat(result,"\n",sep="")
}
close( input )




