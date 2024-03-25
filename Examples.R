
rm(list=ls(all=TRUE)) 
source("Functions-IT-Recurrent.r")

# ============================
# EXAMPLE 1: SIMULATED DATA
# ============================

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



# ===========================
# EXAMPLE 2: THE CGD DATA
# ===========================

dat <- read_excel("dat-CGD.xlsx") %>% as.data.frame
dim(dat); head(dat)

# OVERALL TREATMENT EFFECT
fit <- coxme(Surv(starttime,stoptime,event)~trt +(1|IDNUM), 
            data=dat, control = control.0)
fit0 <- coxme(Surv(starttime,stoptime,event)~(1|IDNUM), 
            data=dat, control = control.0)
anova(fit0,fit,test="LRT")


# INTERACTION TREE
# -------------------
penalty <- c(2:4, log(length(unique(dat$IDNUM))), log(NROW(dat)))
B <- 30; 
split.var <- c(3:8,13) 
max.depth <- 6
# BOOTSTRAP-BASED TREE SIZE SELECTION
set.seed(668)
boot.result <- bootstrap.grow.prune(B=B, data=dat, N0=10, n0=2, split.var=split.var, 
                                    max.depth=6, split.method="lrt", split.strata=FALSE) 
# THE LARGE INITIAL TREE
tree0 <- boot.result$boot.tree[[1]]
# THE PRUNING RESULTS
boot.prune <- boot.result$boot.prune
# OBTAIN THE BEST TREE STRUCTURE WITH DIFFERENT COMPLEXITY PENALTY TERM
G.a <- bootstrap.size(boot.prune, penalty=penalty)
sz <- c(G.a[which.max(G.a[,2]),1], G.a[which.max(G.a[,3]),1], 
        G.a[which.max(G.a[,4]), 1],  G.a[which.max(G.a[,5]),1],
        G.a[which.max(G.a[,6]), 1])
btree <- obtain.btree(tree0, bsize=sz)
tre.s <- btree[[4]]  
tre.s

# PLOT IN R
plot.tree(tre.s, textDepth=3, lines="rectangle", digits = 3)

# GENERATE LATEX TREE FIGURE FILE
plot.tree.latex(tre.s, file="btree-CGD.tex", digits=5)










#
