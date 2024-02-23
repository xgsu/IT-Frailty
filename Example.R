
rm(list=ls(all=TRUE)) # REMOVE ALL PREVIOUS OBJECTS
options(warn=-1);
# options(error = recover)
source("Functions-IT-Recurrent.r")

set.seed(123)
n0 <- 300
beta <- c(-1, 1, 1, 3);
dat <- rdat.exponential(n0=n0, beta=beta, c1.cut=0.5,  c2.cut=0.5, 
                        frailty.dist="lognormal", theta=1, r0=0.25, p0=1/2, 
                        K0=100, d=4, digits=1, c.censor=5) 
print(dim(dat)); print(sort(table(dat$IDNUM), decreasing =TRUE))
penalty <- c(2:4, log(n0), log(NROW(dat)))

split.var <- 7:10; 
boot.result <- bootstrap.grow.prune(B=30, data=dat, N0=20, n0=5, split.var=split.var, 
                                    max.depth=5, split.method="lrt", split.strata=FALSE) 
# THE LARGE INITIAL TREE
tree0 <- boot.result$boot.tree[[1]]; 
# THE PRUNING RESULTS
boot.prune <- boot.result$boot.prune; boot.prune
# OBTAIN THE BEST TREE STRUCTURE WITH DIFFERENT COMPLEXITY PENALTY a
G.a <- bootstrap.size(boot.prune, penalty=penalty); print(G.a)
sz <- c(G.a[which.max(G.a[,2]),1], G.a[which.max(G.a[,3]),1], 
        G.a[which.max(G.a[,4]), 1],  G.a[which.max(G.a[,5]),1],
        G.a[which.max(G.a[,6]), 1])
print(sz) 
btree <- obtain.btree(tree0, bsize=sz)
plot.tree.latex(btree[[4]], file="btree-BIC.tex", 
                            digits=4,cols.nominal=NULL, landscape=FALSE)


















#
