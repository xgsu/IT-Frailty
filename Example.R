
# ============================
# EXAMPLE WITH SIMULATED DATA
# ============================

rm(list=ls(all=TRUE)) 
source("Functions-IT-Recurrent.r")

set.seed(2024)
n0 <- 300
beta <- c(-1, 1, 1, 3);
dat <- rdat.exponential(n0=n0, beta=beta, frailty.dist="lognormal")
dim(dat); head(dat)
penalty <- c(2:4, log(n0), log(NROW(dat)))
split.var <- 7:10; 
boot.result <- bootstrap.grow.prune(B=25, data=dat, N0=10, n0=3, 
    split.var=split.var, max.depth=10, split.method="wald", 
    split.strata=FALSE) 
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
btree.BIC <- btree[[4]]


# PLOT THE FINAL TREE (SELECTED BY BIC)
plot.tree(btree.BIC, cols.nominal=NULL, 
          textDepth=5, lines="rectangle", digits=3)

# OBTAIN LATEX PLOT OF THE FINAL TREE 
plot.tree.latex(btree.BIC, file="btree-BIC.tex", 
          digits=4, cols.nominal=NULL, landscape=FALSE)
# TO COMPILTE THE LATEX FILE WITH THE pstricks PACKAGE, 
# GET dvi USING LATEX; THEN dvips; FINALLY ps2pdf












#
