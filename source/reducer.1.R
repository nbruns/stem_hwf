#! /usr/bin/Rscript
##! /usr/lib64/R/bin/Rscript 
##! /usr/bin/env Rscript  


#########    /usr/lib64/R/bin/Rscript #ATLAS
#########    /usr/bin/env Rscript  #local
#########  /usr/bin/Rscript  #azure


##define the functions you need.

# .libPaths("/home/fs01/neb76/R_libs.new/") ##not seeing these packages! test now.
source("reducer.1.modeling.library.abundance.R")

library(gbm)
options(warn=-1)

#-------------------------------------------------------------------------------
# declar fit critical logic
#-------------------------------------------------------------------------------
##these are the paramters specied in our base code
WEIGHT.TYPE_ZI <- "sqrt"
WEIGHT.TYPE_CALIB <- "equal"
WEIGHT.TYPE_COUNT <- "sqrt"

SHRINKAGE.PARAM <- .01
MIN.STIXEL.SAMPLE.N <- 40
MIN.STIXEL.POSITIVE.N <- 10
GBM_N_TRAIN.TREES <- 1000
GBM_N_PREDICT.TREES <- GBM_N_TRAIN.TREES



trim<-function(p_string){
        return(gsub("(^ +)|( +$)", "", p_string))
}


##touch here when you want to swap out a new basemodel
fit_pred_bucket <- function(train.data){
			bs.index <- as.logical(train.data$bs.sample)
			sampled_data <- train.data[bs.index,]
			cur_model <- abund_ZI_GBM_train(sampled_data,
				min.data.size = MIN.STIXEL.SAMPLE.N,
				min.positive.obs=MIN.STIXEL.POSITIVE.N,
				bag.fraction = .80,
				# gbm complexity
				shrinkage = SHRINKAGE.PARAM,
				n.trees = GBM_N_TRAIN.TREES,
				cv.folds = 0,
				weights.zi.type=WEIGHT.TYPE_ZI,
				weights.calib.type=WEIGHT.TYPE_CALIB,
				weights.count.type=WEIGHT.TYPE_COUNT,
				verbose = F,
				keep.data = T)
				# gbm Sample size: abs min = 25
			 return( cur_model) 
		}

## below appending steps are no longer necessary, because now sending along all NA predictions
bucket_prediction <- function(prediction.data,bucket.model){
	preds <- abund_ZI_GBM_predict(prediction.data,bucket.model,verbost=F,n.trees=GBM_N_PREDICT.TREES)
	## Right here there needs to be an expansion, enforece the 5 column reponse. Or, inside actually the library.
	##tentatively:
	# if(ncol(preds)==1) preds <- cbind(preds,NA,NA,NA,NA)
	return(preds)
}

process_bucket <- function(stixel_id,cur_bucket){
	# names(cur_bucket) <- c("row_id","x","y","t","r","z","e","obs","type","bs.sample")
	loc_info_names <- c("row_id",
		"lon",
		"lat",
		"date")
	covar_names <- c("I.STATIONARY",
		"YEAR",
		"DAY",
		"TIME",
		"EFFORT_HRS",
		"EFFORT_DISTANCE_KM",
		"NUMBER_OBSERVERS",
		"ASTER2011_DEM",
		"UMD2011_LANDCOVER",
		"Water", 
		"Evergreen_Needleleaf", 
		"Evergreen_Broadleaf",
		"Deciduous_Needleleaf", 
		"Deciduous_Broadleaf", 
		"Mixed_Forest", 
		"Woodland",
		"Wooden_Grassland",
		"Closed_Shrubland",
		"Open_Shrubland",
		"Grassland",
		"Cropland",
		"Urban_Built",
		"Bare")

	names(cur_bucket) <- c(loc_info_names,covar_names,"obs","data.type","bs.sample")



	cur_bucket <- data.frame(cur_bucket)
	cur_bucket$bs.sample <- as.logical(as.character(trim(cur_bucket$bs.sample))) #very involved casting of text to logical
	fixed_date <- cur_bucket$date
	fixed_date <- as.numeric(as.character(fixed_date))
	cur_bucket$date <- fixed_date - floor(fixed_date)
	#-------------------------------------------------------------------------------
	# data sorting
		## testing for data existence: 
	#-------------------------------------------------------------------------------
	bucket_by_type <- split(cur_bucket,cur_bucket$data.type)
	#-------------------------------------------------------------------------------
	# model fitting and prediction
	#    data quantity checks
			# training, uses a minimun hoped for amount, 
			# test/srd just need some data, else skip the reporting (otherwise gives strange NA reports)
	#-------------------------------------------------------------------------------
	## count data of each type
	train.n <-nrow(bucket_by_type$train) 
	if(is.null(train.n)) train.n <- 0
	test.n <- nrow(bucket_by_type$test)
	if(is.null(test.n)) test.n <- 0
	srd.n <- nrow(bucket_by_type$srd)
	if(is.null(srd.n)) srd.n <- 0

	if(sum(bucket_by_type$train$bs.sample,na.rm=T) < MIN.STIXEL.SAMPLE.N){
		train.preds <- matrix(NA,nrow=train.n,ncol=5)
		test.preds <- matrix(NA,nrow=test.n,ncol=5)
		srd.preds <- matrix(NA,nrow=srd.n,ncol=5)
	}else{	
		cur.model <- fit_pred_bucket(bucket_by_type$train)	# *** Model fitting is here!! ***
		train.preds <- bucket_prediction(bucket_by_type$train,cur.model)
		if(test.n > 0) test.preds <- bucket_prediction(bucket_by_type$test,cur.model)
		if(srd.n > 0)  srd.preds <- bucket_prediction(bucket_by_type$srd,cur.model)
	}
	#-------------------------------------------------------------------------------
	# report results
	## remember, with streaming, we're just printing out. 
		##so here, you need for format each row, point at a time, hence the seemingly tedeious writing out of each column within the printing of each row...
	#-------------------------------------------------------------------------------
	##choose reporting keys
	train.key.frame <- bucket_by_type$train[,c("row_id","lon","lat","date","obs")] #row_id, lat,lon, time, obs
	train.key.frame <- cbind(train.key.frame,bucket_by_type$train$bs.sample) #add bs.sample for reasembly
	#choose your reporting values along with above keys
	##train has a 5th column, wether in sample or not
		## test and srd have dummy, true values
	if(train.n > 0){
		for(row.iii in 1:nrow(train.key.frame)){
		 	cat(stixel_id,",", train.key.frame[row.iii,1],",",train.key.frame[row.iii,2],",",
		 	train.key.frame[row.iii,3],",",train.key.frame[row.iii,4],",",train.key.frame[row.iii,5],
		 	",",train.key.frame[row.iii,6],",train,",sep="")
		 	cat(train.preds[row.iii,] , sep=",")  
		 	cat("\n")
		}	
	}
	if(test.n > 0){
		test.key.frame <- bucket_by_type$test[,c("row_id","lon","lat","date","obs")]
		for(row.iii in 1:nrow(test.key.frame)){
			cat(stixel_id,",", test.key.frame[row.iii,1],",",test.key.frame[row.iii,2],",",
				test.key.frame[row.iii,3],",", test.key.frame[row.iii,4],",", test.key.frame[row.iii,5],
				",TRUE,test,",sep="")
			cat(test.preds[row.iii,] , sep=",")
		 	cat("\n")  
		}
	}	
	if(srd.n > 0){
		srd.key.frame <- bucket_by_type$srd[,c("row_id","lon","lat","date","obs")]
		for(row.iii in 1:nrow(srd.key.frame)){
			cat(stixel_id,",", srd.key.frame[row.iii,1],",",srd.key.frame[row.iii,2],",",
				srd.key.frame[row.iii,3],",", srd.key.frame[row.iii,4],",", srd.key.frame[row.iii,5],
				",TRUE,srd,",sep="")
			cat(srd.preds[row.iii,],sep=",")
			cat("\n")
		}
	}
}

#-------------------------------------------------------------------------------
# reducing routine
#-------------------------------------------------------------------------------
input <- file( "stdin" , "r" )
# input <- file("post.map.1.txt","r")

# library(debug)
# mtrace(process_bucket)

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
			 bucket <- try(read.csv( tempFile , header=FALSE ),silent=T)
			 if(class(bucket) != "try-error"){
			         process_bucket(lastKey,bucket)
			 }
			# bucket <- read.csv( tempFile , header=FALSE ) 	
			# process_bucket(lastKey,bucket)
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
bucket <- try(read.csv( tempFile , header=FALSE ),silent=T)
##wrap in try, mainly in case the tempFile is empty!  
if(class(bucket) != "try-error"){
        process_bucket(lastKey,bucket)
}
unlink( tempFile )
# bucket <- read.csv( tempFile , header=FALSE )
# process_bucket(currentKey,bucket)
# unlink( tempFile )

close( input )
