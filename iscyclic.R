### input is a square matrix representing a graph
### returns TRUE if the graph is cyclic, returns FALSE otherwise
is.cyclic <- function(mat){
	amat <- matrix(0, nrow(mat),ncol(mat))
	for(i in 1:nrow(mat)){
		for(j in 1:nrow(mat)){
			if(mat[i,j]==2 && mat[j,i]==3){
				amat[i,j]<-1
				amat[j,i]<-0
			}
		} # for j
	} #for i	
	g2 <- as(amat,"graphNEL")
	g.tc <- RBGL::transitive.closure(g2)
	amat.tc <- as(g.tc,"matrix")
	if(any(diag(amat.tc)==1)) return(TRUE)
	else return(FALSE)
} # end function