library(r2r)

n = 10
edges = as.matrix(read.csv("edges.csv"))
nodes = 1:n

makeMat= function(edgeMat, n){
  mat = matrix(0, n, n)
  edges_num = dim(edgeMat)[1]
  for (i in 1:edges_num){
    mat[edgeMat[i,1],edgeMat[i,2]] = edgeMat[i,3]
    mat[edgeMat[i,2],edgeMat[i,1]] = edgeMat[i,3]
  }
  return(mat)
}

isConnected = function(edgeMat, n){
  mat = makeMat(edgeMat, n)
  currEle = 1
  nextExplore = c(1)
  i = 1
  while(i <= length(nextExplore)){
    currEle = nextExplore[i]
    adj = (1:n)[mat[currEle,] != 0]
    nextExplore = c(nextExplore, adj[!(adj %in% nextExplore)])
    i = i+1
  }
  if (i == n+1){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

isCyclic = function(mat,n){
  currEle = 1
  nextExplore = c(1)
  i = 1
  while(i <= length(nextExplore)){
    currEle = nextExplore[i]
    adj = (1:n)[mat[currEle,] != 0]
    if (!is.null(adj %in% nextExplore)){
      return(TRUE)
    }
    nextExplore = c(nextExplore, adj)
    i = i+1
  }
  return(FALSE)
}

prim = function(edgeMat, n){
  if(!isConnected(edges, n)){
    return(FALSE)
  }
  
}

isCyclic(makeMat(edges,n), n)