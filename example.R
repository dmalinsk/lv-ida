## requires packages: pcalg, igraph, RGBL

source("lvida.R")

### sometimes the algorithms in pcalg can return a cyclic graph (for some alpha setting)
### LV-IDA (and IDA) won't work in such cases so this is a script to alert you
source("iscyclic.R")
###

#### here is a very basic simulation: create a DAG, generate some data, search using FCI
#### then estimate some causal effects using LV-IDA

# first, choose number of variables and sample size
pvar <- 8
sample <- 1000
x <- 2
y <- pvar

rDAG <- randomDAG(n = pvar, prob= 0.5, lB = 0.5, uB = 1.5)
data <- rmvDAG(sample,rDAG,errDist="normal")
suffStat <- list(C=cor(data),n=nrow(data))
rules <- c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE) # don't want any possible --- edges
fci.est <- fci(suffStat, indepTest = gaussCItest, p = ncol(data), alpha=0.01, rules=rules)

plot(fci.est)

if(is.cyclic(fci.est@amat)){
  cat("#### FOUND CYCLIC GRAPH #### \n LV-IDA won't work here! \n try again! \n")
}

lv.ida.est <- lv.ida(x,y,cov(data),fci.est@amat,method="local")
lv.ida.est # this is the multiset of causal effects of x = 2 on y = pvar
