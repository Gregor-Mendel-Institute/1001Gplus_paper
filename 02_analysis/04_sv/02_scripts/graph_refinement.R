#' Create and visualize a network graph from an igraph object
#'
#' @param g An igraph object representing the graph.
#' @param label Logical value for node labeling (default is TRUE).
#'
#' @return A ggplot2 object of the visualized network graph.
#'
ggigraph <- function(g, label = T){
  edges.linear = get.edgelist(g)
  g.part <- network(edges.linear, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
  set.seed(20)
  p <- ggnet2(g.part, label = label, edge.color = "black",
              size = 5, color = '#468B97',
              arrow.gap = 0.04, arrow.size = 5,
              # mode = "kamadakawai"
  )
  return(p)
}

#' Create a network graph and visualize it using ggnet2
#'
#' @param edges.linear A list of edges in the form of a matrix (n x 2).
#' @param label Logical value indicating whether to label the nodes (default is TRUE).
#'
#' @return A ggplot2 object representing the visualized network graph.
ggedges <- function(edges.linear, label = T){
  g.part <- network(edges.linear, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
  set.seed(20)
  p <- ggnet2(g.part, label = label, edge.color = "black",
              size = 5, color = '#468B97',
              arrow.gap = 0.04, arrow.size = 5,
              # mode = "kamadakawai"
  )
  return(p)
}


#' Remove direct edges from a graph. WORKING VERSION!!
#'
#' This function takes a list of edges and removes direct edges if a longer path exists.
#' Example list of edges:
#' (1) -> (2)
#' (2) -> (3)
#' (1) -> (3)
#' 
#' The edge (1) -> (3) should de removed, because path (1) -> (2) -> (3) exists.
#'
#' @param b.graph A list of edges in the form of a matrix (n x 2).
#' @param echo Logical value for printing progress (default is TRUE).
#'
#' @return A list of edges without direct edges
#'
refineDirectEdges <- function(edges.compact, echo = T){
  g = igraph::make_graph(t(edges.compact), directed = TRUE)
  E(g)$name <- paste(edges.compact[,1], edges.compact[,2], sep = '-')
  
  nodes.keep = c()
  n.nodes.keep = -1
  while(length(nodes.keep) != n.nodes.keep){
    n.nodes.keep = length(nodes.keep)
    # print(n.nodes.keep)
    print(vcount(g))
    
    deg.in <- igraph::degree(g, mode = "in")
    deg.out <- igraph::degree(g, mode = "out")
    # tails = names(deg.in)[((deg.in + deg.out) == 1) & ((deg.in * deg.out) == 0) | ((deg.in * deg.out) == 1)]
    tails = names(deg.in)[((deg.in + deg.out) == 1) & ((deg.in * deg.out) == 0) ]
    
    # if( 'R4396' %in% tails) stop()
    
    nodes.keep = c(nodes.keep, tails)
    g <- delete_vertices(g, tails)
  }
  
  
  # Remove components with size 1
  g.comp <- igraph::components(g)
  id.single = which(g.comp$csize == 1)
  nodes.single = names(g.comp$membership)[g.comp$membership %in% id.single]
  nodes.keep = c(nodes.keep, nodes.single)
  g <- delete_vertices(g, nodes.single)
  
  
  edges.polised = c()
  g.comp <- igraph::components(g)
  for(i.comp in 1:g.comp$no){
    # message(i.comp)
    names.comp = names(g.comp$membership)[g.comp$membership == i.comp]
    g.sub <- induced_subgraph(g, names.comp)
    
    # print(c(vcount(g.sub), ecount(g.sub)))
    # ggigraph(g.sub)
    
    deg.out <- igraph::degree(g.sub, mode = "out")
    targets = names(deg.out)[deg.out != 0]  
    edges.sub =  get.edgelist(g.sub)
    
    idx.target = edges.sub[,1] %in% targets
    out1 = edges.sub[idx.target, 2]
    out1.comb = paste(edges.sub[idx.target, 1], edges.sub[idx.target, 2], sep = '_')
    
    idx.out1 = which(edges.sub[,1] %in% unique(out1))
    out2 = tapply(edges.sub[idx.out1,2], edges.sub[idx.out1,1], function(x) list(x))
    
    out12 = out2[out1]
    names(out12) = paste(edges.sub[idx.target, 1], '_', out1, '_', sep = '')
    
    out12 = unlist(out12)
    out12.pref = sub("_.*", "", names(out12))
    out2.comb = paste(out12.pref, out12, sep = '_')
    
    out.common = intersect(out1.comb, out2.comb)
    if(length(out.common) != 0){
      edges.delete <- as.matrix(data.frame(do.call(rbind, strsplit(out.common, "_"))))
      
      id.edges.delete <- sapply(1:nrow(edges.delete), function(i) get.edge.ids(g.sub, edges.delete[i,]))
      g.sub <- delete_edges(g.sub, id.edges.delete)
    }
    # print(c(vcount(g.sub), ecount(g.sub)))
    
    # 
    # OLD CODE, not efficient, but still here to remember the fail state.
    # i = 1
    # for(target.vertex in targets){
    #   i = i + 1
    #   if(round(i/100) == i/100) print(i)
    #   out1 <- V(g.sub)[neighbors(g.sub, target.vertex, mode = "out")]$name
    #   out2 <- unlist(sapply(out1, function(s) V(g.sub)[neighbors(g.sub, s, mode = "out")]$name))
    #   v.remove = intersect(out1, out2)
    #   if(length(v.remove) == 0) next
    #   
    #   edge_to_delete <- sapply(v.remove, function(s) get.edge.ids(g.sub, c(target.vertex, s)))
    #   
    #   g.sub <- delete_edges(g.sub, edge_to_delete)  
    # }
    # ggigraph(g.sub, label = F)
    # ggigraph(g.sub)
    
    edges.sub = get.edgelist(g.sub)
    
    edges.polised = rbind(edges.polised, edges.sub)
  }
  
  edges.add = edges.compact[(edges.compact[,1] %in% nodes.keep) | (edges.compact[,2] %in% nodes.keep), ]
  edges.final = rbind(edges.polised, edges.add)
  return(edges.final)
}


#' Remove direct edges from a graph. OLD VERSION!!
#'
#' This function takes a list of edges and removes direct edges if a longer path exists.
#' Example list of edges:
#' (1) -> (2)
#' (2) -> (3)
#' (1) -> (3)
#' 
#' The edge (1) -> (3) should de removed, because path (1) -> (2) -> (3) exists.
#'
#' @param b.graph A list of edges in the form of a matrix (n x 2).
#' @param echo Logical value for printing progress (default is TRUE).
#'
#' @return A list of edges without direct edges
#'
refineDirectEdges_old <- function(b.graph, echo = T){
  # reduce indirect arrows
  idx.remove = c()
  for(i.edge in 1:nrow(b.graph)){
    
    if(echo){
      if(i.edge %% 1000 == 0) print(i.edge)
    }
    
    tmp.to = b.graph[b.graph[,1] == b.graph[i.edge,1],2]
    tmp.from = b.graph[b.graph[,2] == b.graph[i.edge,2],1]
    if(length(intersect(tmp.to, tmp.from)) > 0) idx.remove = c(idx.remove, i.edge)
  }
  idx.remove = unique(idx.remove)
  b.graph = b.graph[-idx.remove,]
  
  return(b.graph)
}