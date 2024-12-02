library(pcalg)
library(igraph)

lm.cov <- function (C, y, x) {
  solve(C[x, x], C[x, y, drop = FALSE])[1, ]
}

find.sink2 <- function(gm) {
  ## Purpose: Find sink of an adj matrix; return numeric(0) if there is none;
  ## a sink may have incident undirected edges, but no directed ones
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gm: Adjacency matrix (gm_i_j is edge from j to i)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006;  speedup: Martin Maechler, Dec.2013
  ## New speedup: DMalinsky, Feb.2017
  
  uncon <- which(colSums(gm) == 0) # added 2.27.2017 to speed things up
  
  ## treat undirected edges
  gm[gm == t(gm) & gm == 1] <- 0
  ## treat directed edges
  setdiff(which(colSums(gm) == 0),uncon)
}

allDags.fast <- function(gm,a,tmp, verbose=FALSE)
{
  ## Purpose: Find all DAGs for a given PDAG
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gm: Adjacency matrix of initial PDAG; only 0-1 entries
  ##   i -> j iff gm(j,i)=1
  ## - a: copy of gm
  ## - tmp: NULL
  ## ----------------------------------------------------------------------
  ## Value:
  ## - one 0/1 adj.matrix per row
  ## Reversion to graph: as(matrix(res[i,],p,p),"graphNEL")
  ## Reversion to wgtMatrix (i->j iff a[j,i]=1): t(matrix(res[i,],p,p))
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date:  7 Apr 2008, 14:08
  if (sum(a) == 0) {
    if (verbose) {
      cat("Last Call - Final Graph: \n")
      print(gm)
      cat("#################### \n")
    }
    tmp2 <- rbind(tmp,c(t(gm)))
    if (all(!duplicated(tmp2))) tmp <- tmp2
  } else {
    sinks <- find.sink2(a)
    if (verbose) {
      cat("Main Call: ################## \n")
      print(gm)
      print(a)
      cat("Sinks: ",sinks,"\n")
    }
    for(x in sinks) {
      if (verbose) cat("Try removing", x," in a.\n")
      gm2 <- gm
      a2 <- a
      if (adj.check(a,x)) {
        inc.to.x <- a[, x] == 1 & a[x, ] == 1
        if (any(inc.to.x)) {
          real.inc.to.x <- as.numeric(row.names(a)[inc.to.x])
          real.x <- as.numeric(row.names(a)[x])
          gm2[real.x, real.inc.to.x] <- 1
          gm2[real.inc.to.x, real.x] <- 0
        }
        a2 <- a[-x,-x]
        if (verbose) {
          cat("Removed sink",as.numeric(row.names(a)[x]),
              "in g (", x,"in a).\n")
          cat("New graphs: \n")
          print(gm2)
          print(a)
        }
        tmp <- allDags.fast(gm2,a2,tmp, verbose)
      }
    }
  }
  tmp
}

adj.check <- function(gm,x) {
  ## Purpose:  Return "TRUE", if:
  ## For every vertex y, adj to x, with (x,y) undirected, y is adjacent to
  ## all the other vertices which are adjacent to x
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gm: adjacency matrix of graph
  ## - x: node number (number)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006;
  ## several smart speedups: Martin Maechler, Dec.2013
  
  gm.1 <- (gm == 1)
  xr <- gm.1[x,]
  xc <- gm.1[,x]
  nx <- which(xr | xc)
  ## undirected neighbors of x
  un <- which(xr & xc)
  for(y in un) {
    adj.x <- setdiff(nx, y)
    adj.y <- setdiff(which(gm.1[y,] | gm.1[,y]), x)
    if(!all(adj.x %in% adj.y))
      return(FALSE)
  }
  TRUE
}

#removed memoisation on 02/04/16 because it doesn't seem to be working. version compatability issue?
#library(memoise) 

#######
#cond1 is a function which checks if condition t1 from Lemma 4.4.2 in Zhang (2006) is passed
#i.e., if A-->B in MAG G is transformed into A<->B, check if there is a directed path from A to B
#in G, other than A-->B.
#######
#f <- memoize(is.path)
foundpath <- c()
is.path <- function(a, b, g, internal = FALSE){
	ind1 <- which(g==3, arr.ind=TRUE, useNames=FALSE)
	ind1 <- subset(ind1,ind1[,2]==a) # array of indices for tails out of A
	if(!internal){
	ind1 <- subset(ind1,!(ind1[,1]==b)) # minus A-->B
	}
	if(nrow(ind1)==0) return(FALSE) # addition: if there are no tails at A aside from the A-->B edge
	for(x in 1:nrow(ind1)){ # loop through tails out of A
		if(g[ind1[x,2],ind1[x,1]]==2){ # if there is an arrowhead at the other end of the x-th tail (call this C)
			if(ind1[x,1]==b){
				foundpath <- append(foundpath,TRUE)
				break
				}
			if(any(g[,ind1[x,1]]==3)){ # if there are any tails out of C, i.e., A-->C--*
				a_old <- a
				a2 <- ind1[x,1]
				if(a2==a_old) next
				foundpath <- append(foundpath,is.path(a2,b,g,internal=TRUE))
				if(any(foundpath)==TRUE) break
			}
		} # if there isn't an arrowhead at C - !(A-->C) - don't return anything
	} # for x in 1:nrow(ind1)
	if(any(foundpath)==TRUE) return(TRUE)
	else return(FALSE)
} # end function

cond1 <- function(tmp_old,pag,q){
	ind3 <- which(pag==1)
	ind3.a <- which(pag==1,arr.ind=TRUE,useNames=FALSE)
	xy3 <- ind3.a[which(ind3==q),]
	a <- xy3[2]
	b <- xy3[1]
	if(is.path(a,b,tmp_old)) return(FALSE)
	else return(TRUE)
} # end function


#######
#cond2 is a function which checks if condition t2 from Lemma 4.4.2 in Zhang (2006) is passed
#i.e., if A-->B in MAG G is transformed into A<->B, check the following:
#if there is a C such that C-->A in G, then C-->B is also in G
#if there is a C such that C<->A in G, then either C-->B or C<->B is also in G
#######
cond2 <- function(tmp_old,pag,q){
	ind <- which(pag==1)
	ind.a <- which(pag==1,arr.ind=TRUE,useNames=FALSE)
	xy <- ind.a[which(ind==q),] # a two-element list which is the i,j coordinate of an element in ind
	# xy[1] is the i coordinate
	# xy[2] is the j coordinate
	# can do things like pag[xy[1],xy[2]]
	if (any(tmp_old[,xy[2]]==2)){ # there is at least one C*->A
		cd.ind <- which(tmp_old[,xy[2]]==2) + (nrow(tmp_old)*(xy[2]-1)) # list of positions of variables with arrowheads into A
		cd.ind2 <- which(tmp_old==2, arr.ind=TRUE, useNames=FALSE)
		cd.ind.a <- subset(cd.ind2,cd.ind2[,2]==xy[2])
		#cd.ind.a <- which(tmp_old[,xy[2]]==2, arr.ind=TRUE, useNames=FALSE) # same list as above but in <row,col> form
		for(f in cd.ind){
			cd.xy <- cd.ind.a[which(cd.ind==f),] # cd.xy[1] is i coordinate, cd.xy[2] is j
			if(tmp_old[cd.xy[2],cd.xy[1]]==3){ # if C-->A
				#cat("THERE IS A C-->A *****************", "\n")
				if(!(tmp_old[cd.xy[1],xy[1]]==2 && tmp_old[xy[1],cd.xy[1]]==3)){
					#cat("BUT NO C-->B !!!! *****************", "\n")
					return(FALSE)
				} #return(FALSE) # if not C-->B
			    else next ##### addition 4/29
			} # if(tmp_old[cd.xy[2],cd.xy[1]]==3)
			if(tmp_old[cd.xy[2],cd.xy[1]]==2){ # if C<->A
				#cat("THERE IS A C<->A **********************", "\n")
				if(!(tmp_old[cd.xy[1],xy[1]]==2 && (tmp_old[xy[1],cd.xy[1]]==3 || tmp_old[xy[1],cd.xy[1]]==2))){
					#cat("BUT NO C<->B OR C-->B !!!! **********************", "\n")
					return(FALSE)
				} #return(FALSE) # if neither C-->B or C<->B
				else next ##### addition 4/29 
			} # if(tmp_old[cd.xy[2],cd.xy[1]]==3)
			else return(FALSE) # if neither of these if-statements are entered... (something is wrong!)
		}
		return(TRUE) # for(f in cd.ind) ##### also added 4/29
	} # if any(tmp_old[,j]==2)
	else return(TRUE)
} # end of function

#######
#cond3 is a function which checks if condition t3 rom Lemma 4.4.2 in Zhang (2006) is passed
#i.e., if A-->B in MAG G is transformed into A<->B, check the following:
#there is no discriminating path for A on which B is the endpoint adjacent to A. 
#######
updateList <- function(path, set, old.list)
{
  ## Purpose: update the list of all paths in the iterative functions
  ## minDiscrPath, minUncovCircPath and minUncovPdPath
  ## ----------------------------------------------------------------------
  ## Arguments: - path: the path under investigation
  ##            - set: (integer) index set of variables to be added to path
  ##            - old.list: the list to update
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 21 Oct 2011; Without for() by Martin Maechler
  c(old.list, lapply(set, function(s) c(path,s)))
}
minDiscrPath.TF <- function(path, pag) ## (pag, a,b,c, verbose = FALSE)
{
  ## Purpose: find a minimal discriminating path for a,b,c.
  ## If a path exists this is the output, otherwise NA
  ## ----------------------------------------------------------------------
  ## Arguments: - pag: adjacency matrix
  ##            - a,b,c: node positions under interest
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 25 Jan 2011; speedup: Martin Maechler
  ## Modified by DMalinsky to return TRUE when a path is found and FALSE otherwise
  ## path ends with ... a *-> b *-*c, b is the edge being discriminated and c is the endpoint
  ## a must be adjacent to b and parent of c
  
  c <- path[3]
  b <- path[2]
  a <- path[1]
  p <- as.numeric(dim(pag)[1])
  visited <- rep(FALSE, p)
  visited[c(a,b,c)] <- TRUE # {a,b,c} "visited"
  ## find all neighbours of a  not visited yet
  indD <- which(pag[a,] != 0 & pag[,a] == 2 & !visited) ## d *-> a
  if (length(indD) > 0) {
    path.list <- updateList(a, indD, NULL)
    while (length(path.list) > 0) {
      ## next element in the queue
      mpath <- path.list[[1]]
      m <- length(mpath)
      d <- mpath[m]
      if (pag[c,d] == 0 & pag[d,c] == 0)
        ## minimal discriminating path found :
        return(TRUE) ## return( c(rev(mpath), b,c) ) ## (previously returned the path)
      
      ## else :
      pred <- mpath[m-1]
      path.list[[1]] <- NULL
      
      
      ## d is connected to c -----> search iteratively
      if (pag[d,c] == 2 && pag[c,d] == 3 && pag[pred,d] == 2) {
        visited[d] <- TRUE
        ## find all neighbours of d not visited yet
        indR <- which(pag[d,] != 0 & pag[,d] == 2 & !visited) ## r *-> d
        if (length(indR) > 0)
          ## update the queues
          path.list <- updateList(mpath[-1], indR, path.list)
      }
    } ## {while}
  }
  ## nothing found:  return
  return(FALSE) ## NA ## (previously returned NA if no path found)
} ## {minDiscrPath}

# DEPRECATED 
# ff <- memoize(is.discr.path) 
# founddpath <- c()
# is.discr.path <- function (path, pag)
# {
#   stopifnot((n <- length(path)) >= 3)
#   if (n > nrow(pag)) return(FALSE)
# 
#   pag <- pag
#   c <- path[1]
#   b <- path[2]
#   a <- path[3]
#   first.pos <- path[n]
#   del.pos <- path[n - 1]
#   indD <- which(pag[first.pos, ] != 0 &
#                 pag[, first.pos] == 2)
#   indD <- setdiff(indD, del.pos)
#   for (d in indD)  if(all(d != path)) {
#     if (pag[c, d] == 0 && pag[d, c] == 0) {
# 	### found discr path There is a discriminating path between d and c for b
# 		cat("There is a discriminating path between:",
#               d, "and", c, "for", b, "\n")
# 		founddpath <- append(founddpath,TRUE)
# 		break
#     }
#     else {
#       if (pag[first.pos, d] == 2 && pag[d, c] == 2 && pag[c, d] == 3) {
#         founddpath <- append(founddpath, is.discr.path(path = c(path, d), pag = pag))
#       	if(any(founddpath)==TRUE) break
#       } ## else : keep 'pag'
#     }
#   } ## for( d )
# 
#   if(any(founddpath)==TRUE) return(TRUE)
#   else return(FALSE)
# }## {discr.path}

cond3 <- function(tmp_old,pag,q){
#need to specify b, c in col coordinates
ind4 <- which(pag==1)
ind4.a <- which(pag==1,arr.ind=TRUE,useNames=FALSE)
xy4 <- ind4.a[which(ind4==q),]
b <- xy4[2]
c <- xy4[1]
pag <- tmp_old ### !!! since is.discr.path is written with 'pag', replace that with the MAG under condisideration
indA <- which((pag[b, ] == 2 & pag[, b] != 0) & (pag[c, ] == 3 & pag[, c] == 2))
if(length(indA)==0) return(TRUE) #if indA is empty return false
for(a in indA){
	# founddpath <- c()
	# pathexists <- c()
	# pathexists <- append(pathexists, is.discr.path(path=c(c,b,a),pag=pag))
	# if(any(pathexists)==TRUE) return(FALSE)
  if(minDiscrPath.TF(c(a,b,c),pag)) return(FALSE) ## modified by DMalinsky 08.2024
}
return(TRUE)
}

####### END DEFINITIONS OF COND1, COND2, COND3

############### BEGIN FUNCTION: listMags


listMags <- function(pag, nMags = 500, method=method, Z_i = NULL, opag = pag){
cat("Entering listMags... \n")
tag <- pag
ccomp <- matrix(0,nrow(pag),ncol(pag)) #empty copy of pag
p <- ncol(pag)
cat("original pag = ", "\n")
print(pag)

ind <- which(pag==1)
lind <- length(ind)
if(lind==0) return(list(pag)) ### if there are no circles...

### CREATE THE TAIL AUGEMENTED GRAPH (TAG) ###

for(i in 1:nrow(tag)) {
	for(j in 1:ncol(tag)){
		#double_circle <- 999
		tag[i,j] <- ifelse(tag[i,j]==1 && tag[j,i]!=1,3,tag[i,j]) # circles on o-> and o-- become tails 
		if(tag[i,j]==1 && tag[j,i]==1) ccomp[i,j]<-ccomp[j,i]<-1
	}
}
#cat("tag (tail augmented graph) = ", "\n")
#print(tag)
#cat("circle component of tag = ", "\n")
#print(ccomp)

### END TAG ###

# NEED TO MAKE SURE CCOMP IS A WGTMATRIX NOT JUST DEFAULT MATRIX
ccomp_g <- as(ccomp,"graphNEL") # make the circle component of the pag into a graph
ccomp_m <- wgtMatrix(ccomp_g) # get the wgtMatrix of the circle component graph (a pattern/CPDAG)
cat("Entering allDags method... \n")
#cat("print circle component: \n")
#print(ccomp_m)
listDags <- allDags.fast(ccomp_m,ccomp_m,NULL) # get all possible DAGs from that pattern
cat("Finished allDags. \n")
if(is.null(listDags)){
	cat("ERROR! Circle compenent of input PAG cannot be oriented into a DAG! \n")
	return(NULL)
}

if(FALSE){
### this stuff is new as of 12.14.2014 (below)###

if(method=="local"){

##function which requires igraph package
connected <- function(ccomp){
	if(!require("igraph")) stop("Package 'igraph' must be installed!")
	ccomp_g2 <- graph.adjacency(ccomp,mode="undirected") # turn ccomp into an igraph object
	return(is.connected(ccomp_g2)) # if the circle component is connected, return TRUE; otherwise, FALSE.
}

## if the circle component is NOT connected...
if(!connected(ccomp)){
	cat("@@@@@@@ not-connected circle component @@@@@@ \n")
  #opag is the original pag
	ccomp_opag <- matrix(0,nrow(opag),ncol(opag))
	for(i in 1:nrow(opag)){
		for(j in 1:ncol(opag)){
			if(opag[i,j]==1 && opag[j,i]==1) ccomp_opag[i,j] <- ccomp_opag[j,i] <- 1
		}
	}
	# now ccomp_opag is the circle component of original pag
	rem <- c()
	for(k in 1:nrow(listDags)){
		ccomp_opag_tmp <- ccomp_opag
		dag <- matrix(listDags[k,],p,p) # current dag in list
		for(i in 1:p){
			for(j in 1:p){
				if(dag[i,j]==1){
					ccomp_opag_tmp[Z_i[[i]],Z_i[[j]]] <- 1
					ccomp_opag_tmp[Z_i[[j]],Z_i[[i]]] <- 0
				}
			}
		}
		# now ccomp_opag_tmp has the dag as a subgraph
		# check if ccomp_opag_tmp is extendable to a dag; if not, then
		# add k to rem (list of dags to remove from allDags)		
		check <- pdag2dag(as(ccomp_opag_tmp,"graphNEL"))	
		if(!check$success) rem <- c(rem,k)
		if(!check$success) cat("@@@@@@@ dag thrown out! @@@@@@ \n")	
	} # for k in listDags
	if(!is.null(rem)) listDags <- listDags[-rem,]
	
} # if !connected

} # if method=="local"

### this stuff is new as of 12.14.2014 (above)###
} # if FALSE temporary as of 11.27.2016

mags <- list() # make a list to fill with MAGs
for(k in 1:nrow(listDags)){
	mags[[k]] <- tag
	for(i in 1:p){
		for(j in 1:p){ # loop thru the list of DAGs
			if(matrix(listDags[k,],p,p)[i,j]==1) # where there are edges in the DAG...
			{
				mags[[k]][i,j] <-2
				mags[[k]][j,i] <-3 # put those edges in a copy of the original TAG
			}
		}
	}
}

###### Now we have a list of MAGs, each of which is a TAG with the circle component oriented as one of the possible DAGs ######


### The following code works roughly like this: 
### the outermost loop (with index k) iterates through the mags that have only 
### invariant double-headed arrows (these are stored in a list called mags). 
### Then for n=1 it looks at all the possible combinations of length n of circle marks
### in the original pag. Then it loops through the matrix entries which correspond to 
### circles in the original pag. For each of these, test if the current mag under study 
### has a tail at the right matrix location. If it does, change it to an arrowhead. 
### This is a single mark change from the mag in mags. Then look at all the possible 
### combinations of circles, length n=2. Try changing both marks. That would be a two-mark 
### different mag. Then try n=3 three mark changes, then 4... Only save those mags that 
### haven't been saved before. ###

bigmaglist <- vector(mode="list", nMags) # default value: the list can hold max 500 mags!!!
for(k in 1:length(mags)){bigmaglist[[k]]<-mags[[k]]}
i <- k+1 # counter to fill the list called bigmaglist
mag_t <- tmp <- tmp_old <- matrix(0,nrow(pag),ncol(pag))
for(k in 1:length(mags)){ # loop over elements in the list of mags
	mag_t <- tmp <- mags[[k]]
	for(n in 1:lind){
	#cat("Entering combn for circle marks... \n")
	comb <- combn(ind,m=n,FUN=NULL,simplify=FALSE)
	#cat("Finished combn. \n")
	### BUG FIX
	### if input to combn(x) is an integer (i.e., ind is of length 1), the it returns seq(1:x)
	### we don't want that! so just let comb <- ind in that case
	if(lind==1) comb <- ind
	###
	for(y in 1:length(comb)){
			tmp <- mag_t
			for(q in comb[[y]]){
				tmp_old <- tmp
				if(mag_t[q]==3){ # if the mark is tail
					tmp[q] <- 2 # change it to an arrowhead
				} # if mag_t[q]==3
			} # for q
			#n <- n+1		
			#cat("Checking for duplicates... \n")
			if(any(sapply(bigmaglist,identical,tmp))) next
			if(any(sapply(bigmaglist,identical,tmp_old))){
				if(i > nMags){
				  cat("MORE MAGS THAN SPECIFIED BY PARAMETER 'nMags' !!!!! ", "\n")
				  bigmaglist <- bigmaglist[!sapply(bigmaglist, is.null)]
				  cat("Finishing listMags. \n")
				  return(bigmaglist)
				  #break
				}
				foundpath <- c()
				#cat("Checking conditions 1, 2, and 3 for transformation... \n")
				if(cond1(tmp_old,pag,q) && cond2(tmp_old,pag,q) && cond3(tmp_old,pag,q)) { # conditions for transformational equivalence....
					#cat("Passed all 3 conditions. \n")
				  bigmaglist[[i]] <- tmp
					i <- i+1
					#cat("saving graph with k = ",k," and m = ", m, "and q = ", q, "\n")
				} # conditions are TRUE		
			} # if tmp_old is in bigmaglist
		} # for y	
	} # for n
} # for k

bigmaglist <- bigmaglist[!sapply(bigmaglist, is.null)]
cat("Finishing listMags. \n")
return(bigmaglist)

} # end function listMags

# DEPRECATED
# is.visible <- function(a,b,g){
# 	#checking edge from a to b (need to have previously checked that there is a directed edge from a to b)
# 	#if there is a vertex c not adjacent to b such that c has an arrowhead into a, then return true
# 	if (any(g[,a]==2)){ # there is at least one C*->A
# 		cd.ind <- which(g[,a]==2) + (nrow(g)*(a-1)) # list of positions of variables with arrowheads into A
# 		cd.ind2 <- which(g==2, arr.ind=TRUE, useNames=FALSE)
# 		cd.ind.a <- subset(cd.ind2,cd.ind2[,2]==a)
# 		for(f in cd.ind){
# 			cd.xy <- cd.ind.a[which(cd.ind==f),] # cd.xy[1] is i coordinate, cd.xy[2] is j
# 			#cat("THERE IS A C*->A *****************", "\n")
# 			if(!(g[cd.xy[1],b]!=0 && g[b,cd.xy[1]]!=0)){
# 				#cat("BUT NO C*-*B !!!! *****************", "\n")
# 				return(TRUE)
# 			} # if there is no C*-*B, then the edge is visible
# 			else{ 
# 				#if there is vertex c such that there is a collider path between c and a that is into a and every non-endpoint vertex on the path is a parent of b, then return true
# 				if(is.discr.path(c(b,a,cd.xy[1]),g)) return(TRUE)
# 				else next
# 			}
# 		} 
# 		return(FALSE) # for(f in cd.ind)
# 	} # if any(g[,a]==2)
# 	else return(FALSE)	
# 	#else return false
# }


remove.visible.edges <- function(x,g){
	#removes visible edges out of node x, returns updated graph
	indX3 <- which(g==3, arr.ind=TRUE, useNames=FALSE) # location of all tails in graph
	indX4 <- subset(indX3,indX3[,2]==x) # only tails out of node x
	indX5 <- subset(indX4,g[x,indX4[,1]]==2) # only directed edges out of x

	for(i in indX5[,1]){
		if(visibleEdge(g,x,i)){ ## if(is.visible(x,i,g)){ ## modified to use visibleEdge from pcalg 08.2024
			g[x,i] <- g[i,x] <- 0 # if directed edge out of x is visible, remove it
		}
		else next
	}
	return(g) # return updated graph
}

### is a an ancestor of b in graph g?
#fff <- memoize(is.ancestor)
is.ancestor <- function(a, b, g){
	if(a==b) return(TRUE) ## fix 7/26... if x=y, x is an ancestor of y
	foundpath <- c()
	ind1 <- which(g==3, arr.ind=TRUE, useNames=FALSE)
	ind1 <- subset(ind1,ind1[,2]==a) # array of indices for tails out of A
	if(nrow(ind1)==0) return(FALSE) # addition: if there are no tails at A
	for(x in 1:nrow(ind1)){ # loop through tails out of A
		if(g[ind1[x,2],ind1[x,1]]==2){ # if there is an arrowhead at the other end of the x-th tail (call this C)
			if(ind1[x,1]==b){
				foundpath <- append(foundpath,TRUE)
				break
				}
			if(any(g[,ind1[x,1]]==3)){ # if there are any tails out of C, i.e., A-->C--*
				a_old <- a
				a2 <- ind1[x,1]
				if(a2==a_old) next
				foundpath <- append(foundpath,is.ancestor(a2,b,g))
				if(any(foundpath)==TRUE) break
			}
		} # if there isn't an arrowhead at C - !(A-->C) - don't return anything
	} # for x in 1:nrow(ind1)
	if(any(foundpath)==TRUE) return(TRUE)
	else return(FALSE)
} # end function

### is a an ancestor of b in graph g?
#fff.2 <- memoize(is.poss.ancestor)
is.poss.ancestor <- function(a, b, g,visited=NULL){
	if(a==b) return(TRUE)
	foundpath <- c()
	ind1 <- which(g==3, arr.ind=TRUE, useNames=FALSE) #tails
	ind11 <- which(g==1, arr.ind=TRUE, useNames=FALSE) #circles
	ind1 <- rbind(ind1,ind11) ## tails and circles
	ind1 <- subset(ind1,ind1[,2]==a) # array of indices for tails and circles out of A
	if(nrow(ind1)==0) return(FALSE) # addition: if there are no tails or circles at A
	for(x in 1:nrow(ind1)){ # loop through tails and circles out of A
		if(ind1[x,1] %in% visited) next
		if(g[ind1[x,2],ind1[x,1]]==2 || g[ind1[x,2],ind1[x,1]]==1){ # if there is an arrowhead or circle at the other end of the x-th tail (call this C)
			if(ind1[x,1]==b){
				foundpath <- append(foundpath,TRUE)
				break
				}
			if(any(g[,ind1[x,1]]==3 | g[,ind1[x,1]]==1)){ # if there are any tails or circles out of C
				a_old <- a
				a2 <- ind1[x,1]
				if(a2==a_old) next
				foundpath <- append(foundpath,is.poss.ancestor(a2,b,g,visited=c(visited,a_old)))
				if(any(foundpath)==TRUE) break
			}
		} # if there isn't an arrowhead at C - !(A-->C) - don't return anything
	} # for x in 1:nrow(ind1)
	if(any(foundpath)==TRUE) return(TRUE)
	else return(FALSE)
} # end function

## Function that computes the set D-SEP(X,Y), modified by DMalinsky from version done by Spirtes
dsepset.reach <- function(a,b,c,adjacency)
{
  ## reachable      set of vertices;
  ## edgeslist      array[1..maxvertex] of list of edges
  ## numvertex      integer
  ## labeled        array (by depth) of list of edges that have been labeled

  makeedge <- function(x,y) list(list(x,y))
  
    legal.dsep <- function(r,s) {
    ## Modifying global 'edgeslist'
   if ((adjacency[r[[1]],r[[2]]] == 2 &&
         adjacency[s,     r[[2]]] == 2 && r[[1]] != s) &&  (is.ancestor(s,a,adjacency) || is.ancestor(s,b,adjacency)) && (is.ancestor(r[[2]],a,adjacency) || is.ancestor(r[[2]],b,adjacency))) {
      edgeslist[[r[[2]]]] <<- setdiff(edgeslist[[r[[2]]]],s)
      makeedge(r[[2]],s)
    }
  }


  initialize.pdsep <- function(x,y) mapply(makeedge, x=x, y=y)

  labeled <- list()
  numvertex <- dim(adjacency)[1]
  edgeslist <- list()
  for (i in 1:numvertex)
    edgeslist <- c(edgeslist,list(which(adjacency[,i] != 0)))
  labeled[[1]] <- initialize.pdsep(a, edgeslist[[a]])
  edgeslist[[a]] <- list()
  depth <- 2
  repeat {
    labeled[[depth]] <- list()
    for (i in seq_along(labeled[[depth-1]])) {
      lab.i <- labeled[[depth-1]][[i]]
      edgestemp <- edgeslist[[lab.i[[2]]]]
      if (length(edgestemp) == 0) break
      for (j in seq_along(edgestemp))
        labeled[[depth]] <- union(legal.dsep(lab.i, edgestemp[[j]]),
                                  labeled[[depth]])
    }
    if (length(labeled[[depth]]) == 0)
      break
    ## else :
    depth <- depth  + 1
  }
  ### BUG FIX 7/26
  ### the function previously included variables adjacent to a automatically
  ### because of the step labeled[[1]] <- initialize.pdsep(a, edgeslist[[a]])
  ### this is fine for possible-d-sep but not for d-sep
  dsep <- unique(unlist(labeled))
  for(k in dsep){
  	if(adjacency[a,k]!=0){
  		if(is.ancestor(k,a,adjacency)||is.ancestor(k,b,adjacency)) dsep <- dsep
  		else dsep <- setdiff(dsep,k)
  	} else dsep <- dsep
  }
  dsep
} # end function

## Function that computes the set Possible-D-SEP(X,Y), modified by DMalinsky from version done by Spirtes
pdsepset.reach <- function(a,b,c,adjacency)
{
  ## reachable      set of vertices;
  ## edgeslist      array[1..maxvertex] of list of edges
  ## numvertex      integer
  ## labeled        array (by depth) of list of edges that have been labeled

  makeedge <- function(x,y) list(list(x,y))
  
    legal.dsep <- function(r,s) {
    ## Modifying global 'edgeslist'
   if (((adjacency[r[[1]],r[[2]]] == 2 &&
         adjacency[s,     r[[2]]] == 2 && r[[1]] != s) || ((adjacency[r[[1]],s] != 0 && r[[1]] != s))) &&  (is.poss.ancestor(s,a,adjacency) || is.poss.ancestor(s,b,adjacency))    && (is.poss.ancestor(r[[2]],a,adjacency) || is.poss.ancestor(r[[2]],b,adjacency))) {
      edgeslist[[r[[2]]]] <<- setdiff(edgeslist[[r[[2]]]],s)
      makeedge(r[[2]],s)
    }
  }


  initialize.pdsep <- function(x,y) mapply(makeedge, x=x, y=y)

  labeled <- list()
  numvertex <- dim(adjacency)[1]
  edgeslist <- list()

  for (i in 1:numvertex)
    edgeslist <- c(edgeslist,list(which(adjacency[,i] != 0)))
  labeled[[1]] <- initialize.pdsep(a, edgeslist[[a]])
  edgeslist[[a]] <- list()
  depth <- 2
  repeat {
    labeled[[depth]] <- list()
    for (i in seq_along(labeled[[depth-1]])) {
      lab.i <- labeled[[depth-1]][[i]]
      edgestemp <- edgeslist[[lab.i[[2]]]]
      if (length(edgestemp) == 0) break
      for (j in seq_along(edgestemp))
        labeled[[depth]] <- union(legal.dsep(lab.i, edgestemp[[j]]),
                                  labeled[[depth]])
    }
    if (length(labeled[[depth]]) == 0)
      break
    ## else :
    depth <- depth  + 1
  }
  dsep <- unique(unlist(labeled))
  dsep
} # end function

reach <- function(a,b,c,adjacency)
{
  ## reachable      set of vertices;
  ## edgeslist      array[1..maxvertex] of list of edges
  ## numvertex      integer
  ## labeled        array (by depth) of list of edges that have been labeled

  makeedge <- function(x,y) list(list(x,y))

  legal.pdsep <- function(r,s) {
    ## Modifying global 'edgeslist'
    if ((adjacency[r[[1]],r[[2]]] == 2 &&
         adjacency[s,     r[[2]]] == 2 && r[[1]] != s) ||
        (adjacency[r[[1]],s] != 0 && r[[1]] != s)) {
      edgeslist[[r[[2]]]] <<- setdiff(edgeslist[[r[[2]]]],s)
      makeedge(r[[2]],s)
    }
  }

  initialize.pdsep <- function(x,y) mapply(makeedge, x=x, y=y)

  labeled <- list()
  numvertex <- dim(adjacency)[1]
  edgeslist <- list()
  for (i in 1:numvertex)
    edgeslist <- c(edgeslist,list(which(adjacency[,i] != 0)))
  labeled[[1]] <- initialize.pdsep(a, edgeslist[[a]])
  edgeslist[[a]] <- list()
  depth <- 2
  repeat {
    labeled[[depth]] <- list()
    for (i in seq_along(labeled[[depth-1]])) {
      lab.i <- labeled[[depth-1]][[i]]
      edgestemp <- edgeslist[[lab.i[[2]]]]
      if (length(edgestemp) == 0) break
      for (j in seq_along(edgestemp))
        labeled[[depth]] <- union(legal.pdsep(lab.i, edgestemp[[j]]),
                                  labeled[[depth]])
    }
    if (length(labeled[[depth]]) == 0)
      break
    ## else :
    depth <- depth  + 1
  }
  unique(unlist(labeled))
}

### taken from pcalg package 11/27/2016 ###
possibleDe <- function(amat,x)
{
  ## Purpose: in a DAG, CPDAG, MAG, or PAG determine which nodes are
  ##          possible descendants of x on definite status paths
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - amat: matrix corresponding to the DAG, CPDAG, MAG, or PAG
  ## - x: node of interest
  ## ----------------------------------------------------------------------
  ## Value:
  ## - de.list: array containing the possible descendants of x
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 26 Apr 2012, 16:58
  
  stopifnot(is.matrix(amat))
  p <- nrow(amat)
  is.de <- rep.int(FALSE, p) ##
  ## 1. case: x is a possible child of itself
  is.de[x] <- TRUE
  ## 2. case: find all the possible children of x
  indD <- which(amat[x,] != 0  & amat[,x] != 2 & !is.de) ## x (o,-)-* d
  i.pr <- rep(x,length(indD))
  while (length(indD) > 0) {
    ##next element in the queue
    d <- indD[1]
    indD <- indD[-1]
    pred <- i.pr[1]
    i.pr <- i.pr[-1]
    is.de[d] <- TRUE
    a.d <- amat[,d]
    a.d.p <- a.d[pred]
    ## find all possible children of d not visited yet
    indR <- which(amat[d,] != 0 & a.d != 2 & !is.de) ## d (o,-)-* r
    for(j in seq_along(indR)) {
      ## check that the triple <pred,d,r> is of a definite status
      ## 1. d is a collider on this subpath; this is impossible
      ##    because the edge between d and r cannot be into d
      ## 2. d is a definite non-collider
      r <- indR[j]
      if (a.d.p == 3 || a.d[r] == 3 ||
          (a.d.p == 1 && a.d[r] == 1 && amat[pred,r] == 0)) {
        ## update the queues
        indD <- c(indD, r)
        i.pr <- c(i.pr, d)
      }
    }
  }
  ## return 'de.list' :
  which(is.de)
  
} ## {possibleDe}


###################################################################
###################################################################
###################################################################
###################################################################
lv.ida <- function(x.pos,y.pos,mcov,pag,method="global",nMags=500, localcap=NULL,
                possdsep="small", verbose=FALSE, mags.local=FALSE,Z_i=NULL,opag=pag,bugwatch=FALSE)
{
  ## Purpose: Estimate the causal effect of x on y; the graphEst and correlation
  ## matrix have to be precomputed; all MAGs can be precomputed
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - x.pos, y.pos: Column of x and y in d.mat
  ## - mcov: Covariance matrix that was used to estimate graphEst
  ## - pag : PAG matrix found by FCI **this is necessary for the time being because there isn't a way to get from the graphical object to the matrix**
  ## - graphEst: Fit of FCI Algorithm (a PAG)
  ## - method: "local" - local
  ##           "global" - all MAGs
  ## - verbose: if TRUE, details on regressions that were used
  ## - all.mags: All MAGs in the format of function listMAGs; if this is
  ##   available, no new function call allMags is done
  ## ----------------------------------------------------------------------
  ## Value: causal values
  ## ----------------------------------------------------------------------
  if (verbose) cat("Starting...")
  pag <- matrix(pag,nrow(pag),ncol(pag),dimnames = NULL)
  #amat <- pag
  #x.pos <- as.character(x.pos)
  #y.pos <- as.character(y.pos)
  if (verbose) cat("Loaded the adjacency matrix", "\n")
  if (method=="local") {
##############################
    ## local method
    ## Main Input: mcov, graphEst
##############################
	
	beta.hat <- c()
	
	if(all(pag[,x.pos]==0) || all(pag[,y.pos]==0)) return(0)
	
	# set of adjacencies of x
	adj <- which(pag[,x.pos]!=0)
	if (verbose) cat("Calculated adjacencies of x", "\n")
	
	# set of possible descendents of x
  pdes <- c()
    #for(k in 1:nrow(pag)){
    #  if (verbose) cat("Ancestor search row = ", k, "\n")
    # 	if(is.poss.ancestor(x.pos,k,pag)){
    #  		pdes <- c(pdes,k)
    #  	}
    #  }
    #
  pdes <- possibleDe(pag,x.pos) ## trying this out 11/27/2016
	if (verbose) cat("Constructed set of possible descendents", "\n")
  if (!(y.pos %in% pdes)) return(0)
	# possible d-sep
	if(possdsep=="small"){
		pdsep <- union(pdsepset.reach(x.pos,y.pos,-1,pag),y.pos) ##added y.pos to this set 8/29
	}
	if(possdsep=="big"){
		pdsep <- union(reach(x.pos,y.pos,-1,pag),y.pos) ##added y.pos to this set 8/29
	}
	if (verbose) cat("Constructed possible-d-sep set", "\n")
  
  ### %%%%%%% MAJOR BUG FIX 08.2024 %%%%%%%% ###
  pdsep.possDe <- c()
  for(i in pdsep){ pdsep.possDe <- c(pdsep.possDe,possibleDe(pag,i))}

  
	# Z_i = union of adj, pdes, and pdsep ## AND possDe(pdsep)
  Z_i <- sort(unique(c(adj,pdes,pdsep,pdsep.possDe))) ## Z_i <- sort(unique(c(adj,pdes,pdsep)))
  ### END BUG FIX ###	

	if (!is.null(localcap) && length(Z_i)>localcap){
	  cat("WARNING: cannot localize when calculating the effect of ", x.pos, " on ", y.pos, 
	      ". Z_i is too big. Size of Z_i = ", length(Z_i), ". Returning NA", "\n")
	  return(NA)
	}
	subpag <- pag[c(Z_i),c(Z_i)]
	mcov <- mcov[c(Z_i),c(Z_i)]
	x.pos <- which(Z_i==x.pos)
	y.pos <- which(Z_i==y.pos)
	if(bugwatch==TRUE){
	  cat("$$ adj = ", adj, "\n")
	  cat("$$ pdes = ", pdes, "\n")
	  cat("$$ pdsep = ", pdsep, "\n")
	  cat("$$ Z_i = ", Z_i, "\n \n \n")
  }
	beta.hat <- lv.ida(x.pos,y.pos,mcov,subpag,method="global",nMags,
                verbose=FALSE, mags.local=TRUE, Z_i = c(Z_i), opag = pag)

  } else {
##############################
    ## global method
    ## Main Input: mcov, pag, graphEst
##############################
    p <- nrow(pag)
    am.pag <- pag

    ## find all MAGs if not provided externally
    am <- if(mags.local) listMags(am.pag,nMags = nMags, method="local",Z_i=Z_i,opag=opag) else listMags(am.pag,nMags = nMags,method="global")
    #cat("mag list = ", "\n")
    #print(am)
    n.mags <- length(am)
    cat("# of MAGs to enumerate = ", n.mags, "\n")
    #if(n.mags == nMags) cat("Number of MAGs listed is maxed out w.r.t. nMags; there might be MAGs in the equivalence class which are missing. Try increasing nMags.")
    beta.hat <- rep(NA,n.mags)
    for (i in 1:n.mags) {
      ## compute effect for every MAG
      
      gMag <- am[[i]]

      ### if x is not an ancestor of y, the causal effect is 0
      if(!is.ancestor(x.pos,y.pos,gMag)){
      	beta.hat[i] <- 0
      	next
      }
      ###
      
      gMag_x <- remove.visible.edges(x.pos,gMag) ### This is the MAG when remove invisible edges out of x

      dsepset <- setdiff(dsepset.reach(x.pos,y.pos,-1,gMag_x),x.pos)
      
      ### the set of descendents of x
      des <- c()
      for(k in 1:nrow(gMag)){
      	if(is.ancestor(x.pos,k,gMag)){
      		des <- c(des,k)
      	}
      }
      if(gMag_x[x.pos,y.pos] != 0) { ### if y is adjacent to x in gMag_x
        beta.hat[i] <- NA
      } else if(length(intersect(dsepset,des)) != 0){ ### if the intersection of dsepset and the descendents of x is non-empty
      	       beta.hat[i] <- NA
      		 } else {
                    beta.hat[i] <- lm.cov(mcov,y.pos,c(x.pos,dsepset))
                    #cat("printing Mag with visible edges removed: ", gMag_x, "\n")
                    cat("dsepset = ", dsepset, "\n")
                    cat("effect estimate = ", beta.hat[i], "\n")
               }
    } ## for ( i  n.dags)
  } ## else : method = "global"
  beta.hat
}
