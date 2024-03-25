
rm(list=ls(all=TRUE)) # REMOVE ALL PREVIOUS OBJECTS
options(warn=-1);
# options(error = recover)
source("Functions-IT-Recurrent.r")

# ------------------
# EXPONENTIAL MODEL 
# ------------------
set.seed(123)
n0 <- 300; 
beta <- c(-1, 1, 1, 3);
c1.cut <- 0.5;  c2.cut <- 0.5
theta <- 1;  r0 <- 0.25; p0 <- 1/2 
K0 <- 100; d <- 4; digits <- 1; 
c.censor <- 5  # IMPORTANT IN DETERMINING THE SAMPLE SIZE
# SIMULATION RUNS
nrun <- 100; B <- 20; 
split.var <- 7:10; 
max.depth <- 5
BOOTS <- BTREE <- as.list(1:nrun)
SIZE <- matrix(0, nrun, 6)
for (i in 1:nrun){
  	print(paste("#################### run = ",i, " ######################", sep="")) 
    dat <- rdat.exponential(n0=n0, beta=beta, c1.cut=c1.cut,  c2.cut=c2.cut, 
                            frailty.dist="gamma", theta=theta, r0=r0, p0=p0, 
                            K0=K0, d=d, digits=digits, c.censor=c.censor) 
  	print(dim(dat)); print(sort(table(dat$IDNUM), decreasing =TRUE))
  	penalty <- c(1:4, log(n0), log(NROW(dat)))

  	boot.result <- bootstrap.grow.prune(B=B, data=dat, N0=20, n0=3, split.var=split.var, 
		        max.depth=max.depth, split.method="wald", split.strata=FALSE) 
  	BOOTS[[i]] <- boot.result
  	# THE LARGE INITIAL TREE
  	tree0 <- boot.result$boot.tree[[1]]; 
  	# THE PRUNING RESULTS
	  boot.prune <- boot.result$boot.prune; boot.prune
	  # OBTAIN THE BEST TREE STRUCTURE WITH DIFFERENT COMPLEXITY PENALTY a
	  G.a <- bootstrap.size(boot.prune, penalty=penalty); print(G.a)
	  sz <- c(G.a[which.max(G.a[,2]),1], G.a[which.max(G.a[,3]),1], 
       	 G.a[which.max(G.a[,4]), 1],  G.a[which.max(G.a[,5]),1],
       	 G.a[which.max(G.a[,6]), 1],  G.a[which.max(G.a[,7]),1])
	  print(sz) 
	  btree <- obtain.btree(tree0, bsize=sz)
	  SIZE[i, ] <- sz; 
	  BTREE[[i]] <- btree
} 
SIZE
save(BOOTS, BTREE, SIZE, file="result3-GammaFrailty-Strong-exp.Rdat") 




# ------------------
# WEIBULL MODEL 
# ------------------

rm(list=ls(all=TRUE)) # REMOVE ALL PREVIOUS OBJECTS
options(warn=-1);
source("Functions-IT-Recurrent.r")

set.seed(123)
n0 <- 300; 
beta <- c(-1, 1, 1, 3);
c1.cut <- 0.5;  c2.cut <- 0.5
theta <- 1;  r0 <- 0.25; p0 <-1/2 
K0 <- 100; d <- 4; digits <- 1; 
c.censor <- 5  # IMPORTANT IN DETERMINING THE SAMPLE SIZE
# SIMULATION RUNS
nrun <- 100
B <- 20; split.var <- 7:10; 
max.depth <- 5
BOOTS <- BTREE <- as.list(1:nrun)
SIZE <- matrix(0, nrun, 6)
for (i in 1:nrun){
  print(paste("#################### run = ",i, " ######################", sep="")) 
  dat <- rdat.weibull(n0=n0, beta=beta, c1.cut=c1.cut,  c2.cut=c2.cut, 
                      frailty.dist="gamma", theta=theta, r0=r0, p0=p0, 
      K0=K0, d=d, digits=digits, c.censor=c.censor)
  print(dim(dat)); print(sort(table(dat$IDNUM), decreasing =TRUE))
  penalty <- c(1:4, log(n0), log(NROW(dat)))
  
  boot.result <- bootstrap.grow.prune(B=B, data=dat, N0=20, n0=3, split.var=split.var, 
                          max.depth=max.depth, split.method="wald", split.strata=FALSE) 
  BOOTS[[i]] <- boot.result
  # THE LARGE INITIAL TREE
  tree0 <- boot.result$boot.tree[[1]]; 
  # THE PRUNING RESULTS
  boot.prune <- boot.result$boot.prune; boot.prune
  # OBTAIN THE BEST TREE STRUCTURE WITH DIFFERENT COMPLEXITY PENALTY a
  G.a <- bootstrap.size(boot.prune, penalty=penalty); print(G.a)
  sz <- c(G.a[which.max(G.a[,2]),1], G.a[which.max(G.a[,3]),1], 
          G.a[which.max(G.a[,4]), 1],  G.a[which.max(G.a[,5]),1],
          G.a[which.max(G.a[,6]), 1],  G.a[which.max(G.a[,7]),1])
  print(sz)
  btree <- obtain.btree(tree0, bsize=sz)
  SIZE[i, ] <- sz; 
  BTREE[[i]] <- btree
} 
SIZE
save(BOOTS, BTREE, SIZE, file="result3-GammaFrailty-Strong-Weibull.Rdat") 

# load("result3-GammaFrailty-Strong-Weibull.Rdat")
# SIZE


# =========================
# SUMMARIZE RESULTS
# =========================

library(tidyverse)
library(knitr)

FILES <- c("result3-GammaFrailty-Strong-exp.Rdat", "result3-GammaFrailty-Strong-Weibull.Rdat")
Models <- c("Exp.Tree", "Weibull.Tree")
Complexity <- c(1:4, "log.n", "log.N")
size.c <- function(x, c=3) sum(x==c, na.rm = TRUE)
nullmodel <- FALSE
n.complexity <- length(Complexity)
TBL.OUT <- NULL
for (k in 1:length(FILES)) {
  filename <- FILES[k]
  load(file=filename); 
  Model <- Models[k] 
  # TREE SIZE
  # size.crct.perc <- apply(SIZE, 2, size.c, c=c0)/nrun
  size.mean <- apply(SIZE, 2, mean)
  size.sd <- apply(SIZE, 2, sd)
  n.NullTree <- apply(SIZE, 2, size.c, c=1)
  
  for (j in 1:n.complexity) {
    index <- j    # 4 FOR BIC
    complexity <- Complexity[j]
    n.x1 <- n.x2 <- n.x1.Or.x2 <- n.x1.x2 <- 0
    for (i in 1:length(BTREE)) {
      print(cbind(file=k, tree=i, index=j))  
      Btree <- BTREE[[i]]
      btree <- Btree[[index]]     
      if (NROW(btree) > 1) {
        split.variables <- na.omit(btree$vname)
        if (is.element("X1", split.variables)) n.x1 <- n.x1 +1
        if (is.element("X2", split.variables)) n.x2 <- n.x2 +1  
        if (is.element("X1", split.variables) && is.element("X2", split.variables))  
          n.x1.x2 <- n.x1.x2 +1 
        if (is.element("X1", split.variables) || is.element("X2", split.variables))  
          n.x1.Or.x2 <- n.x1.Or.x2 +1  
      }
    }
    out <- c(Model=Model, complexity=complexity, n.x1=n.x1, 
             n.x2=n.x2, n.x1.x2=n.x1.x2, n.x1.Or.x2=n.x1.Or.x2, 
             n.NullTree[index], size.mean[index], size.sd[index])
    TBL.OUT <- rbind(TBL.OUT, out)
  }
}
TBL.OUT <- as.data.frame(TBL.OUT)
row.names(TBL.OUT) <- NULL
names(TBL.OUT) <- c("Model", "complexity", "X1", "X2", "X1&X2", "X1orX2", 
                    "nNullTree", "size.mean", "size.sd")
TBL.OUT


  
  
  













#
