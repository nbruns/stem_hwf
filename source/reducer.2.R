#! /usr/bin/Rscript
#! /usr/lib64/R/bin/Rscript 
##! /usr/bin/env Rscript  

#########    /usr/lib64/R/bin/Rscript #ATLAS
#########    /usr/bin/env Rscript  #local
######### /usr/bin/Rscript


##reduce task functionality defined here
	# (ie, what is computed on each group of objects assigned the same key)
options(warn=-1) # no warnings will be printed! 
### these castings are tricky if there is NA mixed in! Because the prediction is now a factor... so, must do an involved casting...


##aha! Ok, it should be expected that everything coming out of reducer 1 should be of fixed length.


process_bucket <- function(point_name,cur_bucket){
	# 4-9,215,0.75,0.9,0.375,TRUE,srd,0.8732796
	# if(length(cur_bucket) == 0) next  ##DEBUG
	# print(point_name)
	# print(str(cur_bucket)) ##DEBUG
	# if(ncol(cur_bucket)==7) cur_bucket <- cbind(cur_bucket,NA,NA,NA,NA) ## this is a kludgey fix for now, can be removed later!
	# cur_bucket <- read.csv('ttt.bucket.csv',header=F) ## for local, my machine, exploring!!
	names(cur_bucket) <- c("sample.id","lon","lat","date","obs","bs.sample","pred.pi","pred.tau","pred.mu",
		"pred.truncated","pred.pi.mu")
	cur_bucket <- data.frame(cur_bucket)
	cur_bucket$date <- as.character(as.numeric(cur_bucket$date))
	cur_bucket$pred.pi <- as.numeric(as.character(cur_bucket$pred.pi)) #step is recquired is NA's are mixed into the preds, in which case the whole thing is coded as a factor instead of numeric. Be sure to caste to charcter beofre numeric, else you just get the factor numbers!
	cur_bucket$pred.tau <- as.numeric(as.character(cur_bucket$pred.tau)) 
	cur_bucket$pred.mu <- as.numeric(as.character(cur_bucket$pred.mu)) 
	cur_bucket$pred.truncated <- as.numeric(as.character(cur_bucket$pred.truncated)) 
	cur_bucket$pred.pi.mu <- as.numeric(as.character(cur_bucket$pred.pi.mu)) 

	sample_avs <- aggregate(cur_bucket,by=list(cur_bucket$sample.id),mean,na.rm=T)
	sample_avs$bs.sample <- as.logical(sample_avs$bs.sample)
	pred.mat <- sample_avs[,c("pred.pi","pred.tau","pred.mu","pred.truncated","pred.pi.mu")]

#mean,median,lower.10kupper.90
	# sample_avs <- na.exclude(sample_avs)
	# cat(sample_avs,"\n")
	name.tuple <- strsplit(point_name[1],"-")[[1]]
	point.location <- paste(cur_bucket[1,c("lon","lat","date","obs")],collapse=",")
	if(name.tuple[1]=="train"){	
		#do in and out of bag seperately
		#a lot happens in this one liner:
			# columns are compuonents of zi model, rows are different points
			# summarize each column, paying note to "sample membership"
		in.bag_summary.mat <- apply(pred.mat[sample_avs$bs.sample,],2,col_summary_helper) #5 * 4 matrix
		#now parse and report out!
		in.bag_summary.string <- parse_summary_mat(in.bag_summary.mat) #comma seperated string of those values
		#these are [mu's mean,tau's mean,  .... ] 
		cat(name.tuple[1],".in.bag,",name.tuple[2],",",point.location,",",
			in.bag_summary.string,"\n",sep="")	 

		#repeat on out of bag summaries
		out.bag_summary.mat <- apply(pred.mat[!sample_avs$bs.sample,],2,col_summary_helper)
		out.bag_summary.string <- parse_summary_mat(out.bag_summary.mat)
		cat(name.tuple[1],".out.bag,",name.tuple[2],",",point.location,",",
			out.bag_summary.string,"\n",sep="")	 	
	}else{
		#first, remember srd is 3 part name, vs test with has only 2
		if(name.tuple[1]=="srd"){
			row.id <- name.tuple[3]
		}else{
			row.id <- name.tuple[2]
		}
		summary.mat <- apply(pred.mat,2,col_summary_helper)
		## put option to skip NA's from output, here!
		summary.string <- parse_summary_mat(summary.mat) # flattens matrix insto standarized string!
		cat(name.tuple[1],",",row.id,",",point.location,",",summary.string,"\n",sep="")	 	 	 
	}
	## all the rest is string parsing and packing now pac	
}


#return summary string
## input, point.n * 5 matrix
###			we take column summaries! 4 column summaries
			#  with 5 columns
##output, 5 * 4 matrix
## 		(components of zi model) * measurments in sample
###   (pi,tau,c,trun,smooth) * (mean,median,q10,90)

## then last, not in this function but above you unfurl that into a one long row, by columns
		##yeilding the 20 element vector
		# mean(pi), mean(tau), ... q10(smooth.pred)



col_summary_helper <- function(cur_col){ 	
	if(all(is.na(cur_col))){
		return(rep(NA,4))
	}
	quants <- unlist(quantile(cur_col,probs=c(.1,.5,.9),na.rm=T)) 	#get median, 10th and 90th quatile
	return(c(mean(cur_col,na.rm=T),quants)) #prepend the mean, return the value
}

parse_summary_mat <- function(summary_mat){
	summary_vec <- as.vector(summary_mat) # column order
	paste(summary_vec,collapse=",")
}



#-------------------------------------------------------------------------------
# reducing routine
#-------------------------------------------------------------------------------
input <- file( "stdin" , "r" )
lastKey <- ""
tempFile <- tempfile( pattern="hadoop-mr-demo-" , fileext="csv" ) 
tempHandle <- file( tempFile , "w" )
while( TRUE ){
	 currentLine <- readLines( input , n=1 )
	 if( 0 == length( currentLine ) ){
	 	break
	 }
	 tuple <- unlist( strsplit( currentLine , "\t" ) ) 
	 currentKey <- tuple[1]
	 currentValue <- tuple[2]
	 if( ( currentKey != lastKey ) ){
		 if( lastKey != "" ){ 
			 close( tempHandle )
			#####all logic is here!! ie, critical section-- ie, the computation of note
			bucket <- read.csv( tempFile , header=FALSE ) 	
			process_bucket(lastKey,bucket)
			tempHandle <- file( tempFile , "w" ) 
		 }
		 lastKey <- currentKey 
	 }
	 cat( currentValue , "\n" , file=tempHandle ) #write to temp bin!
}
close( tempHandle ) 

#-------------------------------------------------------------------------------
# process last bucket
#-------------------------------------------------------------------------------
bucket <- read.csv( tempFile , header=FALSE )
process_bucket(currentKey,bucket)
unlink( tempFile )

close( input )
