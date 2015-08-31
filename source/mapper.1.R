#! /usr/bin/Rscript
##! /usr/lib64/R/bin/Rscript
##! /usr/bin/env Rscript 


#########    /usr/lib64/R/bin/Rscript #ATLAS
#########    /usr/bin/env Rscript  #local
#########   #! /usr/bin/Rscript   #Azure
options(warn=-1)

input <- file( "stdin" , "r" )

#-------------------------------------------------------------------------------
# mapping params! (experimental design)
#-------------------------------------------------------------------------------
REP_N <- as.numeric(Sys.getenv("REP_N"))
SAMPLE_N <- as.numeric(Sys.getenv("SAMPLE_N"))
FOLD_N <- REP_N * SAMPLE_N
GRID_DIMS <- c(10,10,40/365)
# REP_N <- 2
# SAMPLE_N <- 2
SUB_SAMPLE_P <- .7
data.type_position <- 29
data.date_position <- 4

all_fold_offsets <- read.csv('fold_offsets.csv',header=F)
##PUT TEST FOR FOLD OFFSETS HERE! or just initialize to 200 always
cur_fold_offsets <- all_fold_offsets[1:FOLD_N,]

DEBUG <- F
if(DEBUG){
	SRD_DATE_VEC <- c(.28,0.30,0.32,0.34, 0.36, 0.38,0.40,0.42,.44)  ## profile case! These are late april early may
	# SRD_DATE_VEC <- c(.24,.26)  ##2 weeks for the debug case! exactly a week apart!	
}else{
	SRD_DATE_VEC <- seq(from = 0, to= 1, length= 52 +1)
	SRD_DATE_VEC <- (SRD_DATE_VEC[1:52] + SRD_DATE_VEC[2:(52+1)])/2	
	SRD_DATE_VEC <- round(SRD_DATE_VEC,digits=2)
}



srd_expand_mat <- cur_fold_offsets[rep(1:FOLD_N,length(SRD_DATE_VEC)),]
# srd_repped_vec <- rep(SRD_DATE_VEC,each=FOLD_N)
srd_expand_mat <- cbind(srd_expand_mat,rep(SRD_DATE_VEC,each=FOLD_N))
names(srd_expand_mat) <- c()

#    FRIDAY NOTE::: figure out this apply step! the expanded matrix looks good now!
raw_stixel_id_vec_hash <- function(loc_info,grid_dims,grid_shift){
	abs(floor((loc_info + grid_shift) / grid_dims))
}


##just put an NA in the loc vec! that way conditionals are at the top!
## output, for each point, matrix:
	# (dimensions + 1) x (fold number * srd_week number)
srd_offset_processor <- function(srd_offset_row,loc_vec){
	cur_point <- c(loc_vec[1],loc_vec[2],srd_offset_row[5])
	row_offset <- as.vector(srd_offset_row[2:4])
	return_vec <- raw_stixel_id_vec_hash(cur_point,GRID_DIMS,row_offset)
	return_vec <- c(srd_offset_row[1],return_vec)
	return(return_vec)
}

obs_offset_processor <- function(obs_offset_row,loc_vec){
	return_vec <- raw_stixel_id_vec_hash(loc_vec,GRID_DIMS,obs_offset_row[2:4]) 
	return(c(obs_offset_row[1],return_vec))
}
#-------------------------------------------------------------------------------
# mapping outine!
#-------------------------------------------------------------------------------
while( TRUE ){
	 currentLine <- readLines( input , n=1 ) 
	 if( 0 == length( currentLine ) ){
	 	break
	 }
	 ## past for debug purposes!
			 # currentLine<- "5,-121.95376750389,37.0115581240276,NA,0,2013,NA,7,1,1,1,41,8,0,16.6667,2.7778,0,0,0,0,0,75,0,0,2.7778,2.7778,0,NA,\"srd\""
			 # currentLine<- "9,-121.1816165,37.0860762,2013.26229508197,0,2013,96,15.17,0.333,0.322,1,234,10,5.2632,0,0,0,0,0,0,0,0,13.1579,73.6842,7.8947,0,0,0,\"train\""
	 currentFields <- unlist( strsplit( currentLine , "," ) ) 
	 #remember, these come through as characters, so caste into numeric!
	 ##currentFields
	 	# "row_id","x","y","t","r","z","e","obs","type"
	 # "1 row_id" , "2 x" , "3 y" , "4 t" , "5 r" , "6 z" , "7 e" , "8 obs" , "9 type"

	 if(currentFields[data.type_position]=="\"srd\""){
	 	hash_func <- srd_offset_processor
	 	hash_mat <- srd_expand_mat
	 	date_vec <- srd_expand_mat[,5]
	 	loc_vec <- c(lon=as.numeric(currentFields[2]),
	 		 	as.numeric(currentFields[3]),
	 		 	NA) 
	 }else{
	 	hash_func <- obs_offset_processor
	 	hash_mat <- cur_fold_offsets
	 	circle_date <- as.numeric(currentFields[4])
	 	circle_date <- circle_date - floor(circle_date) #here, enforcing year agreggation	
	 	date_vec <- rep(circle_date,FOLD_N)
	 	loc_vec <- c(lon=as.numeric(currentFields[2]),
	 		 	as.numeric(currentFields[3]),
	 		 	circle_date) 
	 }

	 ## do the subsampling if a training point! else use a dummer
	 if(currentFields[data.type_position]=="\"train\"") {
	 	##get specific sampling functions!
	 	raw_draws <- as.logical(rbinom(SAMPLE_N,1,SUB_SAMPLE_P))
	 	sample_vec <- rep(raw_draws,1,each=REP_N)
	 } else{
	 	sample_vec <- T
	 }

	stixel_ids <- t(apply(hash_mat,1,hash_func,loc_vec=loc_vec))

	#now format bringing the other values along for the fold/srd repeates
	#template
	template_row <- data.frame(t(currentFields),stringsAsFactors=F)
	value_frame <- template_row[rep(1,nrow(stixel_ids)),]
	#now do the two special replacements, special dates for SRD, special samples for training
	value_frame[,4] <- date_vec
	value_frame[,ncol(value_frame) +1 ] <- sample_vec
	# ok, all ingredients are together! now parse it into string, report out!	
	## first, convert keys from mat to "-"" seperated string vector
	stixel_id_strings <- apply(stixel_ids,1,function(x){paste(x,collapse="-")})
	## convert values from mat to ","" seperated string vector
	value_row_strings <- apply(value_frame,1,function(x){paste(x,collapse=",")})
	##combine the two string vectors with tabs
	final_string_form<- paste(stixel_id_strings,value_row_strings,sep="\t")
	cat(final_string_form,sep="\n")

} ## done processing all lines!
close( input )
