
# rm(list=ls(all=TRUE)) # REMOVE ALL PREVIOUS OBJECTS
options(warn=-1);
pkgs.needed <- c("survival", "tidyverse", "survminer", "coxme", "ClusterBootstrap", "xtable")
if (any(!is.element(pkgs.needed, installed.packages()[,1]))) install.packages(pkgs.needed)
library(survival)
library(tidyverse)
options(tidyverse.quiet=TRUE)
library(survminer)
library(coxme)
library(ClusterBootstrap) # ALTERNATIVE WAY TO RE-SAMPLE WITH CLUSTERS
library(xtable)

# ===============================================
# FUNCTION rdat() GENERATES RECURRENT TIMES
# ================================================

# n0 = NUMBER OF SUBJECTS/PATIENTS
# K0 = MAXIMUM NUMBER OF RECURRENT EVENTS FOR EACH INDIVIDUALS
# d = # COVARIATES
# c.censor = RANGE OF CENSORING TIME UNIFORM(0, c.censor)
# digits = SIGNIFICANT DIGITS OF ROUNDING, WHICH HELPS REDUCE THE NUMBER OF DISTINCT VALUES OF EACH X
# rdat.exponential() GENREATES DATA WITH EXPONENTIAL GAP TIMES
rdat.exponential <- function(n0=200, beta=c(1, 1, 1, 2), 
                             c1.cut=0.5,  c2.cut=0.5,
                             frailty.dist="lognormal", theta=1, 
                             r0=0.25, p0=1/2, 
                             K0=100, d=4, digits=1, c.censor=2){
	# GENERATE COVARIATES, TRT
	X <- matrix(round(runif(n0*d, 0, 1), digits=digits), nrow=n0, ncol=d) # ROUNDING HELPS REDUCE THE NUMBER OF DISTINCT VALUES OF EACH X
	X <- as.data.frame(X); names(X) <- paste("X", 1:d, sep="") 
	(runif(n0)<p0) %>% as.numeric() -> trt

	# FRAILTY
	if (frailty.dist=="lognormal") {
	  Mu0 <- 0; C.0 <-  theta^2/exp(2*Mu0)
	  sigma <- sqrt(log((1+sqrt(1+4*C.0))/2))    # TO HAVE SCALE (SD) theta
	  eta <- rlnorm(n=n0, meanlog = Mu0, sdlog = sigma)
	} else eta <- rgamma(n0, shape=1/theta, scale= theta)
	# MODEL 
  	linear.pred <- beta[1]*trt + beta[2]*sign(X[,1] <= c1.cut) + beta[3]*sign(X[,2] <= c2.cut)  + 
  	  beta[4]*sign(X[,1] <= c1.cut)*sign(X[,2] <= c2.cut)*trt
#				beta[4]*sign(X[,1] <= c.cut)*trt + beta[5]*sign(X[,2] <= c.cut)*trt

	lambda <- r0*(exp(linear.pred)*eta)

	# GENERATE SURVIVAL/CENSORING TIME
	censor.time <- runif(n0, 0, c.censor) 
	
	INCREM <- matrix(0, n0, K0+1)
	for (j in 1:K0)  INCREM[,j+1] <- rexp(n0, rate=lambda)
	INCREM <- as.data.frame(INCREM); names(INCREM) <- paste("incre", 0:K0, sep="") 
	dat0 <- as_tibble(cbind(IDNUM=1:n0, trt=trt, X, INCREM, censor.time)) 	
	dat0 %>%  pivot_longer(
   			cols = starts_with("incre"),
   			names_to = "rep",
   			names_prefix = "incre",
   			values_to = "time") %>% 
		group_by(IDNUM) %>% 
		mutate(rep = lead(rep),
			start.time=accumulate(time, `+`), 
			time.stop=lead(start.time), 
			status=as.numeric(time.stop <= censor.time), 
			stop.time=pmin(censor.time, time.stop)) %>%
		na.omit() %>%
		group_by(IDNUM) %>% 
		mutate(first= as.numeric(row_number() == min(row_number()[status==0]))) %>% 
		filter((first + status)!=0) %>%
		select(-censor.time, -time, -time.stop, -first) %>%
		select(IDNUM, rep, start.time, stop.time, status, everything()) %>%
		rename(starttime=start.time, stoptime=stop.time, event=status) %>%
	  as.data.frame() -> dat 
	return(dat)
}


# IN FUNCTION rdat.weibull(), WEIBULL BASELINE HAZARDS
rdat.weibull <- function(n0=200, beta=c(-1, 1, 1, 2), 
                         c1.cut=0.5,  c2.cut=0.5, 
                         frailty.dist="lognormal", theta=1,   # SCALE (SD)
                         r0=0.25, p0=1/2, 
                         K0=100, d=4, digits=1, c.censor=2){
  # GENERATE COVARIATES, TRT
  X <- matrix(round(runif(n0*d, 0, 1), digits=digits), nrow=n0, ncol=d) # ROUNDING HELPS REDUCE THE NUMBER OF DISTINCT VALUES OF EACH X
  X <- as.data.frame(X); names(X) <- paste("X", 1:d, sep="") 
  (runif(n0)<p0) %>% as.numeric() -> trt
  
  # FRAILTY
  if (frailty.dist=="lognormal") {
    Mu0 <- 0; C.0 <-  theta^2/exp(2*Mu0)
    sigma <- sqrt(log((1+sqrt(1+4*C.0))/2))    # TO HAVE SCALE (SD) theta
    eta <- rlnorm(n=n0, meanlog = Mu0, sdlog = sigma)
  } else eta <- rgamma(n0, shape=1/theta, scale= theta)
  # MODEL
  linear.pred <- beta[1]*trt + beta[2]*sign(X[,1] <= c1.cut) + beta[3]*sign(X[,2] <= c2.cut)  + 
    beta[4]*sign(X[,1] <= c1.cut)*sign(X[,2] <= c2.cut)*trt
#				 beta[4]*sign(X[,1] <= c.cut)*trt + beta[5]*sign(X[,2] <= c.cut)*trt
 		           
  lambda <- r0*(exp(linear.pred)*eta)^(1/2)
  
  # GENERATE SURVIVAL/CENSORING TIME
  censor.time <- runif(n0, 0, c.censor) 
  
  INCREM <- matrix(0, n0, K0+1)
  for (j in 1:K0)  INCREM[,j+1] <- -log(runif(n0))/lambda^2
  INCREM <- as.data.frame(INCREM); names(INCREM) <- paste("incre", 0:K0, sep="") 
  dat0 <- as_tibble(cbind(IDNUM=1:n0, trt=trt, X, INCREM, censor.time)) 	
  dat0 %>%  pivot_longer(
    cols = starts_with("incre"),
    names_to = "rep",
    names_prefix = "incre",
    values_to = "time") %>%
    group_by(IDNUM) %>% 
    mutate(rep = lead(rep),
           start.time=accumulate(time, .f=function(x, dx)sqrt(x^2 + dx)),   ######
           time.stop=lead(start.time), 
           status=as.numeric(time.stop <= censor.time), 
           stop.time=pmin(censor.time, time.stop)) %>%
    na.omit() %>%
    group_by(IDNUM) %>% 
    mutate(first= as.numeric(row_number() == min(row_number()[status==0]))) %>% 
    filter((first + status)!=0) %>%
    select(-censor.time, -time, -time.stop, -first) %>%
    select(IDNUM, rep, start.time, stop.time, status, everything()) %>%
    rename(starttime=start.time, stoptime=stop.time, event=status) %>%
    as.data.frame() -> dat 
  return(dat)
}








# ==========================================================================================  #
#  TREE-RELATED FUNCTION - THE MAIN PART
# ==========================================================================================  #

control.coxph.0 <- coxph.control(eps =1e-03, toler.chol = 1.818989e-012, 
                        iter.max = 100, toler.inf = 0.001, outer.max = 100) 

control.0 <- coxme.control(eps = 1e-08, toler.chol = 1.818989e-012,iter.max=5, refine.method="control")
                                                                                              
# control.0<-coxme.control(eps = 1e-04, toler.chol = 1.818989e-012,iter.max=20, refine.method="control")                                                                                                                                                                
# control.0<-coxme.control(eps = 1e-08, toler.chol = .Machine$double.eps^0.75,iter.max=20, 
           # inner.iter = Quote(max(4, fit0$iter+1)),sparse.calc = NULL,optpar = list(method = "BFGS", 
            # control=list(reltol = 1e-5)),refine.df=4, refine.detail=FALSE, 
           #refine.method="control",sparse=c(50, .02),varinit=c(.02, .1, .4, .8)^2, corinit = c(0, .3))


# --------------------------
# SPLITTING STATISTIC 
# --------------------------

split.stat <- function(dat, z, n0=2, method=NULL, strata=FALSE)
# THREE CHOICES FOR method: c("lrt", "score", "wald"), WITH DEFAULT lrt
{  
 # write.csv(dat,file="data_split_stat_coxph.csv")
  #write.csv(z,"z1.csv")
  
    IDNUM<-dat$IDNUM;starttime <- dat$starttime;stoptime <- dat$stoptime;event <- dat$event; trt <- dat$trt; test.stat <- NA; n <- nrow(dat)
    
    if (sum(is.na(z))>0) {print(cbind(dat, z)); stop("Something is wrong! There should no missing value in z.")}
    if (length(starttime)!=length(z)) stop("the vector z must have the same length as data.") 
    if (length(stoptime)!=length(z)) stop("the vector z must have the same length as data.") 
      
    
    ## select the subset of data with at least 1 recurrent event for each id
    dat_new <- data.frame(cbind(dat,z))
    # print(dat_new)
    
    if (FALSE) {
    		n11 <-length(unique(as.data.frame(dat_new%>% 
                                        group_by(IDNUM) %>%
                                        filter(trt==0 & z==1) %>%
                                        ungroup)$IDNUM)) 
    		n12 <- length(unique(as.data.frame(dat_new%>% 
                                               group_by(IDNUM) %>%
                                               filter(trt==0 & z==0) %>%
                                               ungroup)$IDNUM)) 
    		n21 <- length(unique(as.data.frame(dat_new%>% 
                                               group_by(IDNUM) %>%
                                               filter(trt==1 & z==1) %>%
                                               ungroup)$IDNUM))  
    		n22 <- length(unique(as.data.frame(dat_new%>% 
                                               group_by(IDNUM) %>%
                                               filter(trt==1 & z==0) %>%
                                               ungroup)$IDNUM)) 
    
    		tmp1<-aggregate(event ~ IDNUM, data=dat_new, FUN=function(x) with(rle(x), sum(lengths>=1 & values==1)))
    		id_sel<-tmp1[tmp1$event==1,]$IDNUM
    		dat_sub<-dat_new[dat_new$IDNUM %in% id_sel,]
    
    		M00<-length(unique(as.data.frame(dat_sub %>% 
                                       group_by(IDNUM) %>%
                                       filter(z==0 & trt==0) %>%
                                       ungroup)$IDNUM))
    
    		M01<-length(unique(as.data.frame(dat_sub %>% 
                                       group_by(IDNUM) %>%
                                       filter(z==0 & trt==1) %>%
                                       ungroup)$IDNUM))
    
    		M10<-length(unique(as.data.frame(dat_sub %>% 
                                       group_by(IDNUM) %>%
                                       filter(z==1 & trt==0) %>%
                                       ungroup)$IDNUM))
    
    		M11<-length(unique(as.data.frame(dat_sub %>% 
                                       group_by(IDNUM) %>%
                                       filter(z==1 & trt==1) %>%
                                       ungroup)$IDNUM))
    }

	suppressMessages(dat_new %>% filter(z==1) %>%
		group_by(trt) %>%
		count(IDNUM, event) %>%
		summarize(n.L=sum((event*n)==1)) -> tbl.L)

   suppressMessages(dat_new %>% filter(z==0) %>%
		group_by(trt) %>%
		count(IDNUM, event) %>%
		summarize(n.R=sum((event*n)==1)) -> tbl.R)
	M0 <- ifelse((length(tbl.L$n.L) + length(tbl.R$n.R) < 4),  0, 
	             min(tbl.L$n.L, tbl.R$n.R))
	

  if (is.null(method)) method <- "lrt"# DEFAULT
    test.stat <- NA 
    # if (min(n11, n12, n21, n22) >=n0 && min(M00, M01, M10,M11) >= 5 ) { 
    if (M0 >= n0){  
		if (method=="lrt"){	
			if (strata) {
        			fit1 <- coxme(Surv(starttime,stoptime,event)~trt+(1|IDNUM)+ z:trt + strata(z), control=control.0)
        			fit0 <- coxme(Surv(starttime,stoptime,event)~trt+(1|IDNUM)+ strata(z), control=control.0)
			} else {
			  fit1 <- coxme(Surv(starttime,stoptime,event)~trt+z+ z:trt+ (1|IDNUM),control=control.0)
			  fit0 <- coxme(Surv(starttime,stoptime,event)~trt+z+ (1|IDNUM), control=control.0)
        			#fit1<-coxph(Surv(starttime,stoptime,event)~trt+z+ z:trt+ frailty.gamma(IDNUM),control=control.coxph.0)
			# fit0 <- coxph(Surv(starttime,stoptime,event)~trt+z+ frailty.gamma(IDNUM),control=control.coxph.0)
			}
			test.stat <- 2*(fit1$loglik[2] - fit0$loglik[2])
		} else if (method=="score") {
			if (strata) {
				fit0 <- coxme(Surv(starttime,stoptime,event)~trt + strata(z)+(1|IDNUM), control=control.0)
        			                    fit1 <- coxme(Surv(starttime,stoptime,event)~trt + z:trt + strata(z)+(1|IDNUM),init=c(fit0$coef,0), eps= 1e-20, iter.max = 1)
			} else {
				fit0 <- coxme(Surv(starttime,stoptime,event)~trt+z+(1|IDNUM), control=control.0)
                                                                                fit1 <- coxme(Surv(starttime,stoptime,event)~trt+z + z:trt+(1|IDNUM), eps= 1e-20, iter.max = 1)
			}
        		                    test.stat <- fit1$score
		} else if (method == "wald") {
			if (strata) {
				fit1 <- coxme(Surv(starttime,stoptime,event)~trt +z:trt + strata(z)+(1|IDNUM),control=control.0)
				test.stat <- (coef(fit1)[2])^2/fit1$var[2,2] 
			} else {
				# fit1 <- coxme(Surv(starttime,stoptime,event)~trt+ z + z:trt+(1|IDNUM)) # ,control=control.0 )
        error <- tryCatch(expr={coxph(Surv(starttime,stoptime,event)~trt+ z + z:trt+frailty.gamma(IDNUM)); FALSE}, # control=control.coxph.0),
                         error = function(error) TRUE, 
                         finally = {
                           # message("There could have been an error here with coxph().")
                           NULL
                           # print(c(tbl.L$n.L, tbl.R$n.R));   
                           # print(dat_new) ####################### CHECK POINT ###########
                          }) 
        if (error) test.stat <- NA
        else {
          fit1 <- coxph(Surv(starttime,stoptime,event)~trt+ z + z:trt+frailty.gamma(IDNUM))
				  # test.stat <- as.numeric((coef(fit1)[3])^2/vcov(fit1)[3,3])
          test.stat <- (coef(fit1)[3])^2/fit1$var[3,3]
          # print(c(tbl.L$n.L, tbl.R$n.R));
          # print(test.stat)    ############################# CHECK POINT ############
        }
			}
		} else stop("Wrong choice for the method= argument!")
    }
   if (!is.na(test.stat) && test.stat < 0) {
      # print("Hmmm, you got a negative LRT/score/Wald test statistic!"); 
      test.stat <- 0
     }
    return(as.numeric(test.stat))
}



# ONE SINGLE SPLIT USING THE T TEST
partition.INT <- function(dat, test, name="1", min.ndsz=10, n0=2, split.var, max.depth=10, mtry=length(split.var), 
	split.method="lrt", split.strata=FALSE)
{
 
    call <- match.call(); out <- match.call(expand = F)
    out$info <- out$name.l <- out$name.r <- out$left <- out$right <- out$... <- NULL
    out$name.l <- paste(name, 1, sep=""); out$name.r <- paste(name, 2, sep="")
    n <- length(unique(dat$IDNUM)); n.test <- length(unique(test$IDNUM)); xcol <- NA; cut <- NA; max.score <- -1e20; miss <- NA; score.test <- NA   
    
    # cox0 <- coxme(Surv(starttime,stoptime,event)~trt+(1|IDNUM),data=dat) 
    cox0 <- coxph(Surv(starttime,stoptime,event)~trt+frailty.gamma(IDNUM),data=dat) 
    trt.effect <- as.numeric(exp(coef(cox0))); # trt.effect is the hazard ratio;   
    # trt.effect <- (data.frame((summary(cox0))$coefficients))[1,2]; # trt.effect is the hazard ratio;   
    
    depth <- nchar(name) # CONTROL THE MAX TREE DEPTH
    vnames <- colnames(dat)
    if (depth < max.depth && n >= min.ndsz) {
        for(i in sample(split.var, size=mtry, replace=F)) {
          x <- dat[,i]; temp <- sort(unique(x));  
            if(length(temp) > 1)         { 
                zcut <- temp[-length(temp)]
                  for(j in zcut) {
                    score <- NA;miss_<-NA
			  if (sum(is.na(x)) > 0) {        ##if there are missing values in x
                    	grp.L <- sign(x <= j|is.na(x)) 
			  	score.L <- split.stat(dat, z=grp.L, n0=n0, method=split.method, strata=split.strata)
                    	grp.R <- sign(x > j|is.na(x))
			  	score.R <- split.stat(dat, z=grp.R, n0=n0, method=split.method, strata=split.strata)
				score <- max(score.L, score.R)
				miss_ <- ifelse(score.L>=score.R, "L", "R")
			   } else {
			      grp <- sign(x <= j)
            score <- split.stat(dat, z=grp, n0=n0, method=split.method, strata=split.strata)
                    	# print(score)
                    }
                    # print(cbind(xcol=i, cut=j, score=score))  ####################
            if (!is.na(score) && score >= max.score) {max.score <- score; xcol <- i; cut <- j; miss<-miss_;}
} }} }
    
    if (!(is.na(xcol))) {
	  if (is.na(miss)) {
		grp.test <- sign(test[, xcol] <= cut)
		if (sum(is.na(grp.test))>0) {  # RANDOMLY ASSIGN 0 &1 IN THIS SPECIAL CASE
			grp.test[is.na(grp.test)] <- sample(x=0:1, size=length(sum(is.na(grp.test))), replace = TRUE)
			print("Randomly assign 0 and 1 to missing values in this special case")
		}
	  }
	  else if (miss=="L") grp.test <- sign(test[, xcol] <= cut|is.na(test[, xcol]))
	  else if (miss=="R") grp.test <- sign(test[, xcol] > cut|is.na(test[, xcol]))
	  else print("Nothing else. You should not see this. Something is wrong!")
        score.test <- split.stat(test, z=grp.test, n0=n0/2, method=split.method, strata=split.strata)
    }
    if (!is.na(score.test)){
	  if (is.na(miss)) {
	                            left.dat <- which(dat[,xcol]<= cut);  
                               out$left  <- dat[left.dat,]; out$right <- dat[-left.dat,]
                                test_xcol<-test[,xcol]
                                  if(sum(is.na(test_xcol))>0)
                                         {
                                      test_xcol[is.na(test_xcol)] <- sample(x=test_xcol[!is.na(test_xcol)], size=length(sum(is.na(test_xcol))), replace = TRUE)  
                                             }
          		left.test <- which(test_xcol<= cut); 
		out$left.test <- test[left.test,]; out$right.test <- test[-left.test,]
	   } else if (miss=="L") {
		left.dat <- which(dat[,xcol]<= cut|is.na(dat[, xcol]));
		out$left  <- dat[left.dat,]; out$right <- dat[-left.dat,]
		left.test <- which(test[,xcol]<= cut|is.na(test[, xcol]));  
		out$left.test <- test[left.test,]; out$right.test <- test[-left.test,]
	   } else if (miss=="R") {
		right.dat <- which(dat[,xcol] > cut|is.na(dat[, xcol]));
		out$left  <- dat[-right.dat,]; out$right <- dat[right.dat,]
		right.test <- which(test[,xcol] > cut|is.na(test[, xcol]));  
		out$left.test <- test[-right.test,]; out$right.test <- test[right.test,]
	   } else print("Nothing else. You should not see this. Something is wrong!")
	} else {xcol <- NA; cut <- NA; max.score <- NA; miss <- NA}  
    	vname <- ifelse(is.na(xcol), NA, vnames[xcol])	
    	# score=ifelse(max.score==-1e20|is.na(max.score), NA, max.score)
    	# print(score);print(score.test);print(n.test);print(trt.effect)
    	out$info <- data.frame(node=name, size = n, xcol = xcol,vname=vname, cut= cut,miss=miss,
            score=ifelse(max.score==-1e20|is.na(max.score), NA, max.score), score.test=score.test, size.test=n.test, trt.effect=trt.effect)
    	out 
}
# EXAMPLE
if (FALSE){
 split <- partition.INT(dat, dat, name="1", min.ndsz=20, n0=5, split.var=7:10, max.depth=4, mtry=2, 
	split.method="lrt", split.strata=TRUE)
 print(split$info)
}



grow.INT <- function(data, test, min.ndsz=10, n0=2, split.var, max.depth=10, mtry=length(split.var), 
	split.method="lrt", split.strata=FALSE)
{
    out <- list.nd <- list.test <- temp.list <- temp.test <- temp.name <- NULL
    list.nd <- list(data); list.test <- list(test)
    name <- 1
    while (length(list.nd)!=0) {    
      for (i in 1:length(list.nd)){
        if (!is.null(dim(list.nd[[i]])) && nrow(list.nd[[i]])>1){ 
        split <- partition.INT(list.nd[[i]], list.test[[i]], name[i], min.ndsz=min.ndsz, n0=n0, 
			split.var=split.var, max.depth=max.depth, mtry=mtry, 
			split.method=split.method, split.strata=split.strata)
        # print(split$info)
        out <- rbind(out, split$info)
        # out <- rbind(out, c(as.numeric(paste(rep(depth.tre,depth.tre), collapse="")), n.split, split$info))
        if (!is.null(split$left) && !is.null(split$left.test)){
            temp.list <- c(temp.list, list(split$left, split$right))
              temp.test <- c(temp.test, list(split$left.test, split$right.test))
            temp.name <- c(temp.name, split$name.l, split$name.r)
        }}}
        list.nd <- temp.list; list.test <- temp.test; name <- temp.name
        temp.list <- temp.test <- temp.name <- NULL
    }   
    out$node <- as.character(out$node)
    out <- out[order(out$node), ]
    out
}
# EXAMPLE
if (FALSE){
 tree0 <- grow.INT(dat, dat, min.ndsz=20, n0=10, split.var=7:10, max.depth=4, mtry=2, 
	split.method="lrt", split.strata=TRUE)
 tree0 
}



# ====================================================================
# Pruning and Size Selection Based on LeBlanc and Crowley (JASA, 1992)
# ====================================================================

prune.size <- function(tre)
{
     if(is.null(dim(tre))) stop("No Need to Prune Further.")
     result <- NULL; n.tmnl <- sum(is.na(tre[,4])); subtree <- 1            
     while (n.tmnl > 1 ) {
            # if (n.tmnl==6) {btre <- tre; print(btre)}
            internal <- tre$node[!is.na(tre$cut)]; l <- length(internal); 
            r.value <- 1:l
            for(i in 1:l) {
                branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
                score <- as.numeric(as.vector(branch$score))
                r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
            }
            alpha <- min(r.value)
            nod.rm <- internal[r.value == alpha]; 
            # if (length(nod.rm)>1) print("Multiple Nodes will be pruned. Check!")
            G <- sum(as.numeric(as.vector(tre$score)), na.rm=T);
            G.test <- sum(as.numeric(as.vector(tre$score.test)), na.rm=T)
            result <- rbind(result, cbind(subtree=subtree, node.rm=nod.rm, size.tree=nrow(tre), 
                    size.tmnl=nrow(tre)-l, alpha=alpha, G=G, G.test=G.test))
            tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
            tre[match(nod.rm, tre$node), 3:8] <- NA          ############## REVISED 6/2021 ########
            n.tmnl <- sum(is.na(tre$cut))
            subtree <- subtree + 1          
      }
      # HANDLE THE NULL TREE WITH THE ROOT NODE ONLY
      result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre), 
                    size.tmnl=1, alpha=9999, G=0, G.test=0))    
      result <- as.data.frame(result)
      result
}




de <- function(x, tree)
{
    if(length(x) != 1) stop("The length of x in function de must be 1.")    
    y <- tree$node;  de <- NA
    if(sum(match(x, y), na.rm = T) != 0) {
        temp <- 1:length(y)
        start <- match(x, y) + 1    
        end <- length(y)
        if(start <= length(y) & nchar(y[start]) > nchar(x)) {
            temp1 <- temp[temp >= start & nchar(y) <= nchar(x)][1] - 1
            if(!is.na(temp1)) end <- temp1
            de <- y[start:end]
    }}
    de
}


# ==========================================================================================  #
#  FUNCTIONS RELATED TO THE BOOTSTRAP FUNCTION
# ==========================================================================================  #


boot.cluster <- function(x,id){
  boot.id <- sample(unique(id), replace=T)
  out <- lapply(1:length(boot.id), function(newid){cbind(x[id%in%boot.id[newid],],newid)})
  return(do.call("rbind",out) )
}


bootstrap.grow.prune <- function(B=30, data, N0=10, n0=2, split.var, max.depth=10, 
	split.method="lrt", split.strata=FALSE)  
{
    call <- match.call(); out <- match.call(expand = F)
    out$boot.tree <- out$boot.prune <- out$... <- NULL
    time.start <- date()
    
    tree0 <- grow.INT(data=data, test=data, min.ndsz=N0, n0=n0, split.var=split.var, 
                      max.depth=max.depth,mtry=length(split.var), 
		split.method=split.method, split.strata=split.strata);  
    print(tree0);  
    prune0 <- prune.size(tree0); 
    boot.tree <- list(tree0); boot.prune <- list(prune0) 
    
	# set.seed(1234)
    # fit.fake <- clusbootglm(trt~1, B=B, data=data, clusterid=IDNUM)
    for (b in 1:B) {
      print(paste("###################### b = ", b, " ###########################", sep=""))
      data.all=boot.cluster(data,data$IDNUM)
      data.all$IDNUM=data.all$newid
      data.all$newid <- NULL
      dat=data.all
      # write.csv(dat,"bootsrap_file.csv")
      # print(head(dat))
      # samp <- sample(1:nrow(data), size=nrow(data),  replace=T) 
      # dat <- data[samp, ];
      # U<-unique(data$IDNUM)
      # samp <- sample(U, size=length(U),replace=T)
      # print(unique(samp))
      # dat<-data[data$IDNUM %in% samp, ]
      # print(dat)            # CHECK POINT #

    	# dat <- clusbootsample(fit.fake, samplenr=b)   # b-TH BOOTSTRAP SAMPLE DATA
    	tre <- grow.INT(data=dat, test=data, min.ndsz=N0, n0=n0, split.var=split.var, max.depth=max.depth,mtry=length(split.var), 
		split.method=split.method, split.strata=split.strata); 
    	# print(tre)   #################################
      boot.tree <- c(boot.tree, list(tre)); 
      prune <- prune.size(tre); # print(prune)
      boot.prune <- c(boot.prune, list(prune));
    }
    time.end <- date(); 
    print(paste("The Start and End time for ", B, "bootstrap runs is:"))
    print(rbind(time.start, time.end))
    out$boot.tree <- boot.tree
    out$boot.prune <- boot.prune
    out
}   


bootstrap.size <- function(boot.prune, penalty=c(2:4))
{   
    #  COMPUTE THE ALPHA PRIME'S
    prune0 <- boot.prune[[1]] 
    n.subtree <- nrow(prune0); 
    print(n.subtree)
    alpha <- as.numeric(as.vector(prune0$alpha));
    alpha[alpha <0] <- 0
    # temp <- c(alpha[-1], alpha[length(alpha)])    # REVISED 6/2021
    temp <- c(alpha[1], alpha[-length(alpha)])
    alpha.prime <- sqrt(alpha*temp)  
    # print(cbind(alpha, alpha.prime=alpha.prime))
    b <- length(boot.prune)
    G <- as.numeric(as.vector(prune0$G)); 
    size.tmnl <- as.numeric(as.vector(prune0$size.tmnl)); 
    subtree <- as.numeric(as.vector(prune0$subtree)); 
    G.a <- matrix(0, n.subtree, length(penalty))
    bi <- G.honest <- rep(0, n.subtree)
    
    for (i in 1:n.subtree) {
        a <- alpha.prime[i]
        bias <- 0
        for (j in 2:b){
            prune.bs <- boot.prune[[j]]
            alpha.bs <- as.numeric(as.vector(prune.bs$alpha)); 
            g <- as.numeric(as.vector(prune.bs$G)); 
            g.test <- as.numeric(as.vector(prune.bs$G.test)); 
            indx <- 1
            if (sum(alpha.bs <=a)>0) {          
                temp1 <- which.max(which(alpha.bs<=a))
                indx <- ifelse(is.null(temp1), 1, temp1)
            }
            # print(cbind(g, g.test))
            temp2 <- (g-g.test)[indx]
            bias <- bias + temp2 
            # print(cbind(i, a, j, bias, indx, temp2))
        }
        bi[i] <- bias/(b-1)
        G.honest[i] <- G[i] - bi[i]
        G.a[i,] <- G.honest[i] - penalty* (size.tmnl[i]-1)
    }
    G.honest[n.subtree] <- 0      # REVISED ON 6/2021
    cbind(size.tmnl, G.a, G=G, bias=bi, G.honest=G.honest)
}



# THIE FUNCTION HELPS OBTAIN THE BEST FINAL TREE STRUCTURE
obtain.btree <- function(tre, bsize=c(4,5))
{
     if(is.null(dim(tre))) stop("No Need to Prune Further.")
     result <- NULL; n.tmnl <- sum(is.na(tre[,4])); subtree <- 1 
     lenth <- length(bsize) 
     btree <-as.list(1:lenth)            
     while (n.tmnl >= (min(bsize)) ) {
            if (is.element(n.tmnl,bsize)) {
                index <- (1:lenth)[bsize==n.tmnl]
                for (j in index) btree[[j]] <- tre
            }
            if (n.tmnl ==1) n.tmnl <- 0
            else {
            internal <- tre$node[!is.na(tre$cut)]; l <- length(internal); 
            r.value <- 1:l
            for(i in 1:l) {
                branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
                score <- as.numeric(as.vector(branch$score))
                r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
            }
            alpha <- min(r.value)
            nod.rm <- internal[r.value == alpha]; 
            G <- sum(as.numeric(as.vector(tre$score)), na.rm=T);
            G.test <- sum(as.numeric(as.vector(tre$score.test)), na.rm=T)
            result <- rbind(result, cbind(subtree=subtree, node.rm=nod.rm, size.tree=nrow(tre), 
                    size.tmnl=nrow(tre)-l, alpha=alpha, G=G, G.test=G.test))
            tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
            tre[match(nod.rm, tre$node), 3:8] <- NA
            n.tmnl <- sum(is.na(tre$cut))
            subtree <- subtree + 1}         
      }
      btree
}



# =========================================== END ===============================================  #
# =========================================================================================
# FUNCTION plot.tree.latex() GENERATES LATEX CODES FOR PLOTTING TREE IN PACKAGE pstricks
# =========================================================================================

plot.tree.latex <- function(tree, file="btree.tex", 
                            digits=4,cols.nominal=NULL, landscape=FALSE)
{
  n.node <- nrow(tree)
  sink(file=file)
  #  SET UP THE LATEX FILE
  cat("\\documentclass[12pt]{article} \n")
  cat("\\usepackage{pstricks,pst-node,pst-tree,lscape} \n")
  cat("\\usepackage[landscape]{geometry} \n")
  cat("\\pagenumbering{gobble} \n\n")   # SUPPRESS PAGE NUMBERING
  cat("\\begin{document} \n \n")
  
  # DEFINE SYMBOLS
  if (landscape) cat("\\begin{landscape} \n")
  cat("\\begin{figure} \n")
  cat("\\begin{landscape} \n")
  cat("\\centering \n")
  cat("\\newcommand{\\Tgcircle}{\\Tcircle[linecolor=black, fillcolor=gray]} \n")
  cat("\\newcommand{\\Tgoval}{\\Toval[linecolor=black, fillcolor=gray]} \n")
  cat("\\newcommand{\\Tgframe}[1]{\\Tr{\\psframebox[linecolor=black, fillstyle=solid, fillcolor=orange]{#1}}} \n")
  # OPTION
  cat("\\psset{nodesep=0.7pt, linecolor=black, treesep=1.2cm, levelsep=1.8cm} \n")
  I0 <- i0 <- NULL
  for (i in 1:n.node) {
    node.i <- as.character(tree[i, 1])
    miss.i<-as.character(tree[i, 6])
    
  if (is.na(miss.i)) s.i<-""
    else if (miss.i=="L")  s.i<-"*"
    else s.i<-"**"
    de.i <- de(node.i, tree) #vector with all node numbers except ith node
    blanks <- paste(rep(" ", (nchar(node.i)-1)*8), sep="")  # 8 SPACES IN ONE TAB		
    n.i <- tree$size[i]
    trt.effect <- round(as.numeric(tree$trt.effect[i]), digits=digits)
  lable.i <- ifelse(!is.null(tree$trt.effect), trt.effect, " "); 	##### THIS MEASURE MAY BE DIFFERENT FOR OTHER TYPES OF TREES
        if (!is.na(de.i[1])) {	# Check the INTERNAL NODE
      if (nchar(node.i)==1 ||  substr(node.i, nchar(node.i), nchar(node.i))=="2") 
        cat(blanks, "\\pstree{\\Tgcircle{~~}} \n",sep = "")
                else cat(blanks, "\\pstree{\\Tgcircle{~~} \\tlput{\\color{blue}",rule.i, "\\hspace{-.6in}}} \n", sep = "")				
      cat(blanks, "{ \n", sep = "") 
      I0 <- c(I0, i)
      i0 <- c(i0, i + length(de.i))
      # UPDATE THE SPLITTING RULE
      vname.i <- tree$vname[i]; col.var <- tree$xcol[i]; 
      cutpoint <- tree$cut[i]
      if (is.element(col.var, cols.nominal)){
        cutpoint <- unlist(strsplit(as.character(cutpoint), split=" "))
        cutpoint <- paste("\\{", paste(cutpoint, collapse=","), "\\}", sep="")
        rule.i <- paste("\\texttt{", vname.i, "}","\\texttt{", s.i, "}","$\\,\\in\\,$", cutpoint, sep="") 
      } else {
        # cutpoint <- strtrim(as.character(cutpoint), width=digits); 				
        cutpoint <- round(as.numeric(cutpoint), digits=digits)
        rule.i <- paste("\\texttt{", vname.i, "}","\\texttt{", s.i, "}", "$\\,\\leq\\,", cutpoint, "$", sep="")
      }
    } else if (substr(node.i, nchar(node.i), nchar(node.i))=="1") { # TERMINAL NODE
      cat(blanks, "\\Tgframe{",  lable.i, "} \\nput{d}{\\pssucc}{\\color{black} \\textbf{", n.i, "}}",
          "\\tlput{\\color{blue} ", rule.i, " \\hspace{-.3in}} \n", sep = "")
    } else cat(blanks, "\\Tgframe{",  lable.i, "} \\nput{d}{\\pssucc}{\\color{black} \\textbf{",  n.i, "}} \n", sep = "")
    if (is.element(i, i0)) {
      rep0 <- rep("}", sum(i0==i))
      node.i0 <- as.character(tree[I0[i0==i][1] , 1])
      blanks0 <- paste(rep(" ", (nchar(node.i0)-1)*8), sep="")  		
      cat(blanks0, rep0, "\n", sep = "") 
    }
  }
  cat("\\end{landscape} \n")
  cat("\\end{figure} \n\n")
  if (landscape) cat("\\end{landscape} \n\n")
  cat("\\end{document} \n")
  sink()  
}
# plot.tree.latex(tree0, file="tree-code.tex", digits=5, cols.nominal=7)






# ===================================================================
# PLOTTING IT TREE STRUCTURE, MODIFIED FROM PETER CALHOUN'S CODES
# ===================================================================

plot.tree <- function(tree, cols.nominal=NULL, 
                      textDepth=3, lines="rectangle", digits=4)
{
  depth<-max(nchar(tree[,1]))
  par(xaxs='i')
  par(mar=c(1,1,1,1))
  par(xpd=TRUE)
  plot(1, type="n", xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), axes=FALSE,xaxs="i",yaxs="i")
  nodes <- tree$node
  nodesBin <- gsub("1", "0", nodes)
  nodesBin <- gsub("2", "1", nodesBin)
  lastObs<-nchar(nodesBin)
  nodesBin<-substr(nodesBin,2,lastObs)
  vname <- tree$vname 
  col.var <- tree$xcol 
  cut <- as.character(tree$cut) 
  size <- tree$n  
  effect <- tree$trt.effect
  
  for(i in 1:length(nodesBin)){
    nChar<-nchar(nodesBin[i])
    if(!is.na(vname[i])){
      if(lines=="rectangle"){
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar)/(depth+1)))
        
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
        
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
        
      } else if(lines=="triangle"){
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
        
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
      }         
      
      if(nChar <= textDepth){ 
        if (is.element(col.var[i], cols.nominal)){
          cutpoint <- unlist(strsplit(as.character(cut[i]), split=" "))
          cutpoint <- paste("{", paste(cutpoint, collapse=","), "}", sep="")
          text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2), (depth-nChar)/(depth+1)+1/(depth+20), 
               bquote(.(as.character(vname[i]))%in%.(cutpoint)),cex=1, col="blue")
        } else {
          cutpoint <- round(as.numeric(cut[i]), digits=digits)
          text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2), (depth-nChar)/(depth+1)+1/(depth+20), 
               bquote(.(as.character(vname[i]))<=.(cutpoint)),cex=1, col="blue")
        }
      }
    } else {
      if(nChar <= textDepth){
        effect.i <- round(as.numeric(effect[i]), digits=digits)
        text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2),(depth-nChar)/(depth+1), 
             paste("n=", size[i], sep=""),cex=1, offset=1, col="red")
        text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2), (depth-nChar)/(depth+1)-0.025, 
             paste("d=", effect.i, sep=""),cex=1, offset=1, col="red")
      }
    }
  }
}

# EXAMPLE
# plot.tree(tree0, cols.nominal=7, textDepth=5, lines="rectangle")




# =========================================== END ===============================================  #

# ===================================================================================
# MODIFIED FUNCTIONS partition AND grow FOR VARIABLE IMPORTANCE VIA BOOTSTRAP METHOD 
# METHOD I: OUT-OF-BAG ESTIMATE. 
# ===================================================================================

# dimnames(sample.oob)
# x <- senddown(sample.oob, tre.b, n0=4/2)
# x <- senddown(permute, tre.b, n0=4/2)
# x$tree

senddown <- function(data, tre, n0, split.method="lrt", split.strata=TRUE)
{
    call <- match.call(); out <- match.call(expand = F)
    out$tre <- out$data  <- out$... <- NULL
    tre <- cbind(tre, chi=NA) # This is a way of checking the codes. 
    dat <- cbind(data, node=1)  #add a new column named node to the data set to put tre$node[1] value
    for (i in 1:nrow(tre)){
        if (!is.na(tre$xcol[i])){  
            var.split <-  dat[,tre$xcol[i]];     
            in.node <- (dat$node)==(tre$node[i]);# Assign true if dat$node= tre splitting node
            IDNUM <- dat$IDNUM[in.node]
            starttime <- dat$starttime[in.node]
            stoptime<- dat$stoptime[in.node]
            trt <- dat$trt[in.node]
            event <- dat$event[in.node]
          dat.0 <- data.frame(IDNUM=IDNUM,starttime=starttime,stoptime=stoptime, event=event, trt=trt) # create a data set with the samples in the splitting node
		# HANDLE MISSING VALUES FOR TREE
            
           		if (is.na(tre$miss[i])) z <- sign(var.split[in.node] <= tre$cut[i])
		else if (tre$miss[i]=="L") z <- sign(var.split[in.node] <= tre$cut[i] | is.na(var.split[in.node]))
		else if (tre$miss[i]=="R") z <- sign(var.split[in.node] > tre$cut[i] | is.na(var.split[in.node]))
		else print("Hmmm. You should not be there. Something is wrong!")
            tre$chi[i] <- split.stat(dat=dat.0, z=z, n0=n0, method=split.method, strata=split.strata)
            #print(cbind(i, xcol=tre$xcol[i], cut=tre$cut[i], tre$score[i], tre$score.test[i], tre$chi[i], sum(in.node)))

		# MISSING IN DATA
           
		if (is.na(tre$miss[i])){
		  if(sum(is.na(var.split[in.node]))>0)
		  {
		    var.split[is.na(var.split[in.node])] <- sample(x=min(var.split[!is.na(var.split[in.node])]):max(var.split[!is.na(var.split[in.node])]), size=length(sum(is.na(var.split[in.node]))), replace = TRUE)
		    print("Randomly assign values to missing values in this special case")
		  }            	      	
		          l.nd <- dat$node[in.node & (var.split <= tre$cut[i])] 
            	dat$node[in.node & (var.split <= tre$cut[i])] <- paste(l.nd, 1, sep="")  
            	r.nd <- dat$node[in.node & (var.split > (tre$cut)[i])] 
            	dat$node[in.node & (var.split >  tre$cut[i])] <- paste(r.nd, 2, sep="") 
		} else if (tre$miss[i]=="L") {
			l.nd <- dat$node[in.node & (var.split <= tre$cut[i] | is.na(var.split))]
			dat$node[in.node & (var.split <= tre$cut[i] | is.na(var.split))] <- paste(l.nd, 1, sep="")
			r.nd <- dat$node[in.node & var.split>tre$cut[i] & !is.na(var.split)]
			dat$node[in.node & var.split>tre$cut[i] & !is.na(var.split)] <- paste(r.nd, 2, sep="")
		} else if (tre$miss[i]=="R") {
		  l.nd <- dat$node[in.node & var.split<= tre$cut[i] & !is.na(var.split)]
		  dat$node[in.node & var.split<= tre$cut[i] & !is.na(var.split)] <- paste(l.nd, 1, sep="")
		  r.nd <- dat$node[in.node & (var.split > tre$cut[i] | is.na(var.split) )]
		  dat$node[in.node & (var.split > tre$cut[i] | is.na(var.split))] <- paste(r.nd, 2, sep="")
		  		  		  		  		  					} else print("Hmmm. You should not be there. Something is wrong!")
     	   }
    }
    out$data <- dat; out$tree <- tre
    out 
}


