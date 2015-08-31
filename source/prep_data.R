#! /usr/lib64/R/bin/Rscript
##! /usr/bin/env Rscript 
work.station_clue <- Sys.getenv("HOME")
if(work.station_clue== "/Users/nicholasbruns"){
	work.station_root <- "/Users/nicholasbruns/repos/work_repos/"
}else{
	work.station_root <- "/home/fs01/neb76/"
	}

ebird_abund_dir <- paste(work.station_root,"hadoop_sandbox/ebird_abund/",sep="")
source(paste(ebird_abund_dir,"source/mapper.1.partitioning.library.R",sep=""))


# run.name <- "TRES_very.small.cali"
run.name <- "TRES_nth.amer"
# run.name <- as.numeric(Sys.getenv("RUN_NAME"))
run.dir <- paste(ebird_abund_dir,"runs/",run.name,"/",sep="")
dir.create(run.dir)
dir.create(paste(run.dir,"/data",sep=""))
dir.create(paste(run.dir,"/results",sep=""))
dir.create(paste(run.dir,"/results/hadoop_results",sep=""))




#-------------------------------------------------------------------------------
# mapping params!
#-------------------------------------------------------------------------------
REP_N <- as.numeric(Sys.getenv("REP_N"))
SAMPLE_N <- as.numeric(Sys.getenv("SAMPLE_N"))
# REP_N <- 50
# SAMPLE_N <- 4
SUB_SAMPLE_P <- .7

##need 3 data types
## spatial extent stated (grab from the parameter file)
###  (or just do mins and maxs here!) I like that, minimize the munging


setwd(work.station_root)
# stem_v3_munge_string <- "stem_v3/runs/TRES_very.small.cali/data/small.cali_Tree_Swallow"  #this is a kludgey way to grab the necessary information to quickly grab the datat I want
stem_v3_munge_string <- "stem_abund/runs/TRES_nth.amer_smooth.method.run.full/data/nth_amer_Tree_Swallow"  #this is a kludgey way to grab the necessary information to quickly grab the datat I want
# stem_v3_munge_string <- "stem_v3/runs/Canada_warbler_debug.erd.2013/data/small_center_US_Canada_Warbler"  #this is a kludgey way to grab the necessary information to quickly grab the datat I want

#-------------------------------------------------------------------------------
# munge in the train, test and srd data!
#row_name,loc-position,covariates,observation,data.type
#-------------------------------------------------------------------------------
load(paste(stem_v3_munge_string,'.erd.train.data.RData',sep=""))
options(stringsAsFactors=F)
train.frame <- data.frame(
	row_id=1:length(erd.train$dates),
	erd.train$locs,
	erd.train$dates,
	erd.train$X,
	erd.train$y,
	type.vev=as.character("train") 
	)
write.table(train.frame,file=paste(run.dir,"data/ebird.abund_",run.name,"_erd.train.data.csv",sep=""),
	row.names=F,
	col.names=F,
	sep=",")

write.table(train.frame[1:600,],file=paste(run.dir,"data/debug.size_ebird.abund_",run.name,"_erd.train.data.csv",sep=""),
	row.names=F,
	col.names=F,
	sep=",")

##  ##  ##
## test set!
##	##  ##
load(paste(stem_v3_munge_string,'.erd.test.data.RData',sep=""))
test.frame <- data.frame(
	row_id=1:length(erd.test$dates),
	erd.test$locs,
	erd.test$dates,
	erd.test$X,
	erd.test$y,
	type.vev=as.character("test") 
	)
write.table(test.frame,file=paste(run.dir,"data/ebird.abund_",run.name,"_erd.test.data.csv",sep=""),
	row.names=F,
	col.names=F,
	sep=",")
test.n <- dim(test.frame)[1]
write.table(test.frame[sample(1:test.n,60,replace=F),],file=paste(run.dir,
	"data/debug.size_ebird.abund_",run.name,"_erd.test.data.csv",sep=""),
	row.names=F,
	col.names=F,
	sep=",")


load(paste(stem_v3_munge_string,'.srd.3km.data.RData',sep=""))
##use the srd to get the extent
##this will be filled in to be same size as others
srd.frame <- data.frame(
	row_id=1:length(srd$locs$x),
	srd$locs,
	dates=NA,
	#fill in the observer information here on prep side-- could be done in, which we'll fill in on the 
	I.STATIONARY=0,
	YEAR=2013,
	DAY=NA,
	TIME=7.0,
	EFFORT_HRS=1.0,
	EFFORT_DISTANCE_KM=1.0,
	NUMBER_OBSERVERS=1.0,
	srd$X,
	obs=NA,
	type.vev=as.character("srd") 
	)
write.table(srd.frame,file=paste(run.dir,"data/ebird.abund_",run.name,"_srd.data.csv",sep=""),
	row.names=F,
	col.names=F,
	sep=",")

srd.n <- dim(srd.frame)[1]
write.table(srd.frame[sample(1:srd.n,1200,replace=F),],file=paste(run.dir,"data/debug.size_ebird.abund_",run.name,"_srd.data.csv",sep=""),
	row.names=F,
	col.names=F,
	sep=",")




##now, erd i suppose we'll have to do carefully.
	## i think, first get the train/test mapping out
	## and then? get the srd mapping out, where every week is given an id
		#I should think a little bit more about how to best organize the srd predictions...

grid_dims <- c(10,10,40/365)
REP_N <- 50
SAMPLE_N <- 4
# REP_N <- as.numeric(Sys.getenv("REP_N"))
# SAMPLE_N <- as.numeric(Sys.getenv("SAMPLE_N"))
FOLD_N <- REP_N * SAMPLE_N

## build up the grid shift
fold_offsets <- matrix(rep(grid_dims , FOLD_N),nrow=FOLD_N,byrow=T)  ## set up a blank matrix with the dimension sized
fold_offsets <- apply(fold_offsets,1,function(x) x*runif(n=3, min=0, max=1) )  ## multiply by random 0-1 piece to get a random offset
## output of apply gives rows=dimensions, columns=folds
fold_offsets <- t(fold_offsets) ## store the transpose, rows as folds, columns as dimensions
fold_offsets <-cbind(1:FOLD_N,fold_offsets)
write.table(fold_offsets,
	file=paste(run.dir,"data/fold_offsets.csv",sep=""),
	row.names=F,
	col.names=F,
	sep=","
	)
