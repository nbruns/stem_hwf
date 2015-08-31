#use library and not require, to force exit!


## prediction functions:
	# now give garuntee of same sized output


library(PresenceAbsence)
library(verification)
library(gbm)
library(scam)
library(splines)

abund_ZI_GBM_train <- function(	D,
					min.data.size = 40,
					model.formula = NULL,
					weights.zi.type = "equal", weights.calib.type = "equal", weights.count.type = "equal",
					count.dist = "poisson", int.depth_zi = 10, int.depth_count = 10,
					min.gbmnode_zi = 10, min.gbmnode_count=10,
					min.positive.obs=10, percent.cutoff=0.0175, ss.cutoff=3000,
					add_zero = FALSE,calib.mod.epsilon=.01,
					seperation.case_calibrated.threshold.default=0.5, 
					uncalibrated.preds.sub.epsilon.case_default.occurence.thresh_raw=1.1,
					plot.calibration=FALSE, plot.calibration.name=NULL,
					...){


	med.count <- d.gbm.zi <- d.gbm.count <- calib.mod <- occurence_threshold_calibrated <- list(NA)

	# Check all data requirements are sufficient:
	# ---------------------------------------
	suff.data <- TRUE
	if(nrow(D) < min.data.size) suff.data <- FALSE
	D$count <- D$obs
	D$zi <- as.numeric(D$obs > 0)

	if(!suff.data) print("Models not run: Insufficient data to run ZI model")
	if (suff.data){
		# ---------------------------------------
		# Dealing with very uneven datasets. Set all to 0 or 1 to prevent extreme skews:
		# ---------------------------------------
		uneven_data <- FALSE
		if(nrow(D)>=ss.cutoff & any(table(D$zi) < (percent.cutoff*ss.cutoff))) uneven_data <- TRUE
		if(nrow(D)<ss.cutoff & any(table(D$zi) < floor(percent.cutoff*nrow(D)))) uneven_data <- TRUE

		if ( uneven_data ) {
			# Check for cases that are considered to be an
			# Overwhelming evidence of zero in stixel
			# so we force zero.
			if(sum(D$zi==1)<sum(D$zi==0)){
				D$zi <- 0
			} else {
				D$zi <- 1
			}
		} # close if(uneven_data)

		if(is.null(model.formula)) {model.formula <- list(get_model_formula_abund(D, response_var="zi"), 0)}
		if(length(model.formula)!=2) { model.formula <- NULL }
		if(is.null(model.formula)) {model.formula <- list(get_model_formula_abund(D, response_var="zi"), 0)}

		# ---------------------------------------
		# NA Predictor Logic
		# ---------------------------------------
		D <- D[ , apply(is.na(D), 2, sum) != nrow(D)]
		# ---------------------------------------

		# Produce default model formulae using all variables, if model.formula not provided or not sufficient
		# To be 'sufficient', the model.formula object supplied, must be a list of length 2.
		# if no weights given produce equal weight vectors
		# add the weights as an additional column in the dataset to deal
		# with bug in gbm when calling from inside a function

		weights.zi <- create.weight.vector(D$count, weights.zi.type)
		weights.calib <- create.weight.vector(D$count, weights.calib.type)
		weights.count <- create.weight.vector(D$count, weights.count.type)

		D <- cbind(D, weights.zi, weights.calib, weights.count)

		# ---------------------------------------
		# Run binomial zero-inflated model
		# ---------------------------------------
		d.gbm.zi <- tryCatch(gbm(formula=as.formula(model.formula[[1]]),
						data=D,
						distribution="bernoulli",
						weights = weights.zi,
						interaction.depth = int.depth_zi,
						n.minobsinnode = min.gbmnode_zi,
						...),
						error=function(e) e
						)
		if(!inherits(d.gbm.zi, "error")){

			# ---------------------------------------
			# Do calibration. Uses gaussian GAM
			# ---------------------------------------
			if(d.gbm.zi$cv.folds>0) best.iter.zi <- gbm.perf(d.gbm.zi, method="cv", plot.it=FALSE)
			if(d.gbm.zi$cv.folds==0) best.iter.zi <- 1000
			pred.zi <- predict(d.gbm.zi, D, best.iter.zi, type="response")
			obs.zi <- D$zi
			calibration_run <- abund_ZI_calibration_kappa.max(pred.zi=pred.zi, 
				obs.zi=obs.zi,
				weights.calib=weights.calib,
				add_zero = add_zero, 
				calib.mod.epsilon=calib.mod.epsilon,
				seperation.case_calibrated.threshold.default=seperation.case_calibrated.threshold.default,
				uncalibrated.preds.sub.epsilon.case_default.occurence.thresh_raw=uncalibrated.preds.sub.epsilon.case_default.occurence.thresh_raw,
				do.plot=plot.calibration, plot.calibration.name=plot.calibration.name,
					...)
			occurence_threshold_calibrated <- calibration_run$occurence_threshold_calibrated
			occurence_threshold_raw <- calibration_run$occurence_threshold_raw
			calib.mod <- calibration_run$model_calibration
			# if(is.numeric(occurence_threshold_calibrated)){
			# 	print(paste("occurence threshold, on the calibrated scale: ", round(occurence_threshold_calibrated, 2)))
			# }
			# -------------------------------------------------------------
    			# Model counts for those which cross the threshold into count-land
    			#--------------------------------------------------------------
    			# Sometimes insufficient data to fit this model. So a try-catch
    			# sequence attempts a different set of params, if this is the case.
    			# It tries smaller trees with smaller end-nodes. If this doesn't work
			# the median value for all of the observations is used - effectively
			# a tree with no branches.

			# After subsetting to keep only count data, need to re-remove any covariates with all NA
			D.count <- D[which(pred.zi>occurence_threshold_raw),]  # do this comparisson on the raw scale, not calibrated
			#why? because might break if switching to calibrated scale, still not entirely sure of the control structure here
			D.count <- D.count[ , apply(is.na(D.count), 2, sum) != nrow(D.count)]
			med.count <- median(D.count$count)
			if(nrow(D.count)==0) med.count <- 0
			suff.counts <- TRUE
			if(sum(D$zi)<min.positive.obs) suff.counts <- FALSE
			if(nrow(D.count)<min.positive.obs) suff.counts <- FALSE
			# if(nrow(D)>ss.cutoff & (sum(D$zi)/nrow(D))<percent.cutoff) suff.counts <- FALSE
			if(suff.counts){
     				if(model.formula[[2]]==0){
     					model.formula[[2]] <- get_model_formula_abund(D.count, response_var="count")
     				} 
				d.gbm.count <- tryCatch(gbm(formula=as.formula(model.formula[[2]]),
								data=D.count,
								distribution=count.dist,
								weights = weights.count,
								interaction.depth = int.depth_count,
								n.minobsinnode = min.gbmnode_count,
								...),
								error=function(e) e
								)

				## THis error catch is not catching all relevant problem models.
				## Need to expand it to get other errors which aren't triggering this
				## flag.
				refit <- ifelse(inherits(d.gbm.count, "error"), TRUE, FALSE)

				if(refit){
					d.gbm.count <- tryCatch(gbm(formula=as.formula(model.formula[[2]]),
									data=D.count,
									distribution=count.dist,
									weights = weights.count,
									interaction.depth = ceiling(int.depth_count/2),
									n.minobsinnode = ceiling(min.gbmnode_count/2),
									...),
									error=function(e) e
									)

					use.med <- ifelse(inherits(d.gbm.count, "error"), TRUE, FALSE)

					if(use.med){
						d.gbm.count <- "median"
					} # close if(use.med)
				} # close if(refit)
			} # close if(suff.counts)
			if(!suff.counts){
				use.med <- TRUE
				d.gbm.count <- "median"
			}
		} # close if(!inherits(d.gbm.zi, "error"))
		if(inherits(d.gbm.zi, "error")) {
			d.gbm.zi <- list(NA)
			# print("ZI model ran with errors. Likely data not sufficient range")
		}
	} # close if(suff.data)
	return(list("model_ZI"=d.gbm.zi, "model_count"=d.gbm.count,
		"model_calibration"=calib.mod, "calibration_threshold"=occurence_threshold_calibrated, "median_count"=med.count))
}

##now, helpers for the threshold selection

	abund_ZI_calibration_kappa.max <- function(pred.zi, obs.zi,
						weights.calib,
						add_zero = FALSE, calib.mod.epsilon=.01,
						seperation.case_calibrated.threshold.default=0.5, do.plot=FALSE, plot.calibration.name=NULL,
						uncalibrated.preds.sub.epsilon.case_default.occurence.thresh_raw=1.1,THRESH_SELECTION.METRIC="kappa",
						...){

		calib.mod <- occurence_threshold_calibrated <- occurence_threshold_raw <- list(NA)
		# print(paste("if seperation, calibrated threshold: ", seperation.case_calibrated.threshold.default))

		# Do calibration. Uses gaussian GAM at the moment
		# scam library uses a spline forced to be monotonically increasing
		pred.zi.seq <- seq(0, 1, by=0.01)
		if(add_zero){
			orig.pred.zi <- pred.zi; orig.obs.zi <- obs.zi; orig.weights.calib <- weights.calib
			pred.zi <- c(pred.zi, seq(0.05, 0.95, by=0.1))
			obs.zi <- c(obs.zi, rep(0, 10))
			weights.calib <- c(weights.calib, rep(weights.calib[obs.zi==0][1], 10))
		}
		calib.mod.data <- data.frame(obs.zi, pred.zi, weights.calib)
		#calibration model can't run if all preds are zero
		# if(all(pred.zi==0)){
		#here, there's an issue with not quite all zero situations,
			#the belief on 11/18/14: we saw 1 positive with 49 points, but all had identical very small predictions
									#deemed ok for calibration model, but calibration model failed.
									#believed that the model could fit if there was a different prediction for that 1 point, or other points
									#fix: instead of blocking calibration if all are equal to 0, block calibration if all are epsilon close to zero
		calib.pred <- rep(NA,length(pred.zi))
		calib.train.pred <- rep(NA,length(pred.zi))

		if(all(pred.zi<calib.mod.epsilon)){
			# occurence_threshold_raw <- seperation.case_calibrated.threshold.default ##TODO: note with the name issue, how this is a clear and obvious bad fix. Must be excised and fixed.  
			occurence_threshold_raw <- uncalibrated.preds.sub.epsilon.case_default.occurence.thresh_raw
		}else{
			calib.mod <- tryCatch(scam(obs.zi ~ s(pred.zi, k=6, bs="mpi"), weights=weights.calib, gamma=2, data=calib.mod.data),
				error=function(e) e
				)
			if(!inherits(calib.mod, "error")){
				calib.pred<-predict(calib.mod,data.frame(pred.zi=pred.zi.seq))
				calib.train.pred<-predict(calib.mod,data.frame(pred.zi=pred.zi))
			}
			# If the model didn't run, try with reduced degrees of freedom:
			if(inherits(calib.mod, "error")) {
				calib.mod <- tryCatch(scam(obs.zi ~ s(pred.zi, k=4, bs="mpi"), weights=weights.calib, gamma=2, data=calib.mod.data),
					error=function(e) e
					)
				# If the reduced model ran, force a monotonically increasing spline:
				if(!inherits(calib.mod, "error")){
					calib.pred<-predict(calib.mod,data.frame(pred.zi=pred.zi.seq))
					calib.train.pred<-predict(calib.mod,data.frame(pred.zi=pred.zi))
				}
				if(inherits(calib.mod, "error")) {
					g <- ceiling(runif(1, 1, 1000))
					print(paste("random no:", g))
					write.table(calib.mod.data, paste("calib.failed.data.",g, ".txt", sep=""), row.names=F)
					print(paste("Calibration model failed. Dataset written to file ", getwd(), "/calib.failed.data.txt", sep=""))
				}
			}

	    		# Cap predictions at 0 and 1, required for gausssian model:
				###issue! need to calibrate all predictions, not just the the 100 we supply!


	    		calib.pred <- apply(cbind(calib.pred, rep(1, length(calib.pred))), 1, min)
	    		calib.pred <- apply(cbind(calib.pred, rep(0, length(calib.pred))), 1, max)

	    		calib.train.pred <- apply(cbind(calib.train.pred, 
	    			rep(1, length(calib.train.pred))), 1, min)
	    		calib.train.pred <- apply(cbind(calib.train.pred, 
	    			rep(0, length(calib.train.pred))), 1, max)

	    	####///////////// start the splice!
	    		#todo:
	    			# add arguements to function, then all the way up into the parameter file
	    			# 	- selection.metric
	    			# 	-
	    		#check for seperation
				calibration_hash <- ceiling(runif(1, 1, 1000))	
				# cat("calibration_hash is: ", calibration_hash,"\n")
				infinite_index <- is.infinite(calib.train.pred)
				non_infinite_train.preds <- calib.train.pred[!infinite_index]
	    		obs.zi_index <- as.logical(obs.zi[!infinite_index])
	    		if(min(non_infinite_train.preds[obs.zi_index],na.rm=T) > max(non_infinite_train.preds[!obs.zi_index],na.rm=T)){
	    			# print("seperation exists. Use default .5")
	    			# cat("min pos is: ",min(non_infinite_train.preds[obs.zi_index],na.rm=T),"\n")
					# cat("max neg is: ",max(non_infinite_train.preds[!obs.zi_index],na.rm=T),"\n")
					# print(summary(calib.train.pred))
	    			occurence_threshold_calibrated <- seperation.case_calibrated.threshold.default
	    		}else{
		    		occurence_threshold_calibrated <- select_threshold(thresh_cands=pred.zi.seq,calib.train.pred,obs.zi,THRESH_SELECTION.METRIC,calibration_hash)	
	    		}
	    		# occurence_threshold_calibrated <- .5


	    	######\\\\\\\\\\\\
	    		# Interpolate between the two values of calib.pred which are closest to 0.5
	    		# First check for values exactly equal to 0.5.

			if(any(calib.pred==occurence_threshold_calibrated)) {
				occurence_threshold_raw <- pred.zi.seq[which(calib.pred==occurence_threshold_raw)]
			} else {
				if(all(calib.pred<occurence_threshold_calibrated)) {
					occurence_threshold_raw <- 1
				} else {
					if(all(calib.pred>occurence_threshold_calibrated)) {
						occurence_threshold_raw <- 0
					} else {
						pre_0.5 <- max(which(calib.pred < occurence_threshold_calibrated))
	    					occurence_threshold_raw <- pred.zi.seq[pre_0.5] + (pred.zi.seq[pre_0.5+1] - pred.zi.seq[pre_0.5])*(occurence_threshold_calibrated-calib.pred[pre_0.5])/(calib.pred[pre_0.5+1] - calib.pred[pre_0.5])
	    				}
	    			}
			}

			## Collate calibrated and raw predictions for ZI probabilities:
			raw_occ_prob <- pred.zi
			# calibrated_occ_prob <- calib.pred

			# Produce plot of calibration, if specified:
		} # close else{, on wheter you can find a "occurence_threshold_calibrated, which used to get set back to a raw and returned,
		##pay careful attention to the current early escape case-- when all raw preds are below the testing epsilon
			## for now, leave the occurence_threshold_calibrated as NA, set the raw as **.011**
			##this is introducting NA's into the model as a value for the stixel-threshold, that the stored threshold is now held as NA, 
				#so logic in the prediction function is recquired to handle the NA's, as they were creating issues in building matrices
		return(list("model_calibration"=calib.mod, "occurence_threshold_raw"=occurence_threshold_raw, "occurence_threshold_calibrated"=occurence_threshold_calibrated,"calibrated_pred_zi"=calib.pred))
		##I believe i will add another component to this list-- "calibration_threshold_raw", "calibration_threshold_"
		##actually-- doesn't make sense here. calibration_threshold? not what it really is. 
		### occurence_threshold_raw   occurence_threshold_calibrated
		##so-- add a 4th value, 
			# and change the names!
	}


	##prob should ouput why it was selected!
	##values on x, kappa on right! similar to last time this happened.
	select_threshold <- function(thresh_cands,cur_preds,cur_obs,selection_metric,calib_hash=NA){
		pa.df <- data.frame(
			"space.holder.tag",
			obs=cur_obs,
			ppp=cur_preds
			)

		pa.metrics <- presence.absence.accuracy(pa.df,threshold=thresh_cands,na.rm=T,st.dev=F) #na.rm=T is very important, without, was returning a complete vector of NA's...
		if(is.na(pa.metrics)){
			return(NA)	
		}else{
			optimal_thresh_position <- which.max(pa.metrics$Kappa)
			# zi_thresh.selection_diagnostic_plot(thresh_cands,pa.metrics$Kappa,optimal_thresh_position,selection_metric,calib_hash)
			# #plot kappas, show the selection
			return(thresh_cands[optimal_thresh_position])
			
		}
		# tss <- pa.metrics$sensitivity  + pa.metrics$specificity - 1
		# cur_stat.vec <- c(cur_stat.vec,tss)
		# cur_rel.vec <- log.scale_verify.call(
		# 	obs=pa.df$obs,
		# 	pred=pa.df$ppp,
		# 	bin.n= 15
		# 	)
		# cur_stat.vec <- c(cur_stat.vec,cur_rel.vec$bs.uncert)
		# cur_stat.vec <- c(cur_stat.vec,cur_rel.vec$bs.resol )
		# cur_stat.vec <- c(cur_stat.vec,(1 - cur_rel.vec$bs.resol/cur_rel.vec$bs.uncert ))
		# #relative reliability now
		# cur_stat.vec <- c(cur_stat.vec, cur_rel.vec$bs.reliability)
		# cur_stat.vec <- c(cur_stat.vec, cur_rel.vec$bs.reliability/cur_rel.vec$bs.uncert)

		# n
		# optimal_thresh
	}


	get_model_formula_abund <- function(D,response_var="y"){
		non.na.pred.names <- names(D)[ names(D) != "obs" & names(D) != "row_id"& names(D) != 
			"bs.sample" & names(D) != "data.type" & names(D) != "YEAR" & names(D) != "DAY" & 
			names(D) != "UMD2011_LANDCOVER" & names(D) != "lon" & names(D) != "lat" & names(D) != 
			".id" & names(D) != "window_vec" & names(D) != "count" & substr(names(D), 1, 8) != 
			"weights." & substr(names(D), 1, 8) != "count" & substr(names(D), 1, 8) != "zi"]

		model.formula <- paste(non.na.pred.names, collapse="+")
		model.formula <- paste(response_var,"~", model.formula)
		return(model.formula)
	}


	create.weight.vector <- function(count_vector, weight_type){
		if(weight_type=="equal") weight_vector <- rep(1, length(count_vector))
		if(weight_type=="sqrt") weight_vector <- sqrt(count_vector + 1)
		if(weight_type=="log") weight_vector <- log(count_vector + 2)
		if(weight_type=="ten") weight_vector <- ifelse(count_vector==0, 1, 10)
		return(weight_vector)
	}




#-------------------------------------------------------------------------------
# prediction
#-------------------------------------------------------------------------------
# -----------------------------------------------
# Its important for the user to pass n.trees!!!!

# # ---------------------------------
# # OUTPUTS FROM abund_ZI_GBM_predict:
# # 4 column matrix with columns of
# # 1. prediction from ZI model (0-1)
# # 2. calibration threshold for ZI model. All elements the same. (0-1)
# # 3. prediction from count model (0-Inf)
# # 4. overall prediction calculated from columns 1, 2 and 3. (0-Inf)

abund_ZI_GBM_predict <- function(D, model.obj, verbose=T, ...){
	require(gbm)
	require(mgcv)
	n.points <- nrow(D)
	chunk_preds <- matrix(rep(NA,n.points*5), nrow=n.points, ncol=5,dimnames=
			list(NULL,
			c("pred_zi", "cal_thresh", "pred_count", "trunc_pred","smooth_pred")
				))
	if(length(model.obj[[1]])==1) return(chunk_preds)
	if(class(model.obj[[3]])!="scam") return(chunk_preds)
	if(class(model.obj[[1]])!="gbm") return(chunk_preds) 
	# Each element of model.obj is an output from abund_ZI_GBM_train, which is
	# a list with five elements:
	# 1. "model_ZI": GBM ZI model object -- NA if insufficient data in stixel
	# 2. "model_count": GBM count model object -- NA if insufficient data in stixel for
	#    ZI model. If sufficient data for ZI model but insufficient for count model,
	#    this is a character element =="median";
	# 3. "model_calibration": Calibration model - currently GAM object. -- NA if
	# 	 insufficient data in stixel.
	# 4. "calibration_threshold": Calibration cutoff - numeric object with the
	#    threshold for ZI predictions
	# 5. "median_count": median count value from the training dataset
	# if(class(model.obj[[1]])!="gbm") model.na <- TRUE
		# -- so the verbose printing no longer reads for the occurence-model-failures 
	##new condition in place, should be looked at again...
	###debug-- print out situations!
		## can you have non-scam but a count or median object?
	#//\\//\\	
	if(model.obj[[1]]$cv.folds>0) best.iter.zi <- gbm.perf(model.obj[[1]], method="cv", plot.it=F)
	if(model.obj[[1]]$cv.folds==0) best.iter.zi <- 1000
	raw_zi_preds <- predict(model.obj[[1]], D, best.iter.zi, type="response")
	calib_zi_preds <- predict(model.obj[[3]],data.frame(pred.zi=raw_zi_preds))
	calib_zi_preds[calib_zi_preds < 0] <- 0
	calib_zi_preds[calib_zi_preds > 1] <- 1
	chunk_preds[,1] <- calib_zi_preds
	if(is.na(model.obj[[4]])){
		##model.obj[[4]] is the calibrated occurence threshold
		##
		chunk_preds[,2] <- rep(NA, n.points)	
	}else{
		chunk_preds[,2] <- rep(model.obj[[4]], n.points)	
	}
	if(class(model.obj[[2]])=="gbm"){
			if(model.obj[[2]]$cv.folds>0) best.iter.count <- gbm.perf(model.obj[[2]], method="cv", plot.it=F)
			if(model.obj[[2]]$cv.folds==0) best.iter.count <- 1000
			chunk_preds[,3] <- predict(model.obj[[2]], D, best.iter.count, type="response")
	} else {
			chunk_preds[,3] <- rep(model.obj[[5]], n.points)
	}
	chunk_preds[,4] <- ifelse(chunk_preds[,1]>chunk_preds[,2], 1, 0)*chunk_preds[,3]
	chunk_preds[,5] <- chunk_preds[,1] * chunk_preds[,3]
	#\\//\\//
	return(chunk_preds)
}




