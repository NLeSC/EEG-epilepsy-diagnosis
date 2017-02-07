getmatrixfromfile <- function( file ) { # Get weighted sparse matrix.
  return( as.matrix( read.csv( file, row.names = 1 ) ) )
}	
getmatrixfromfile2 <- function( file ) { # Get weighted sparse matrix.
  X <- as.matrix( read.csv( file, row.names = 1 ) )
  X <- as.matrix( Matrix::forceSymmetric( Matrix::Matrix( X ) ) )
  diag( X ) <- 0
  X[ is.na( X ) ] <- 0
  return( X )
}
matrix2graph <- function( m )  { # Convert sparse Matrix to igraph object
  return( igraph::graph.adjacency( m, mode = "undirected", weighted = TRUE, diag = FALSE ) )
}
graph2matrix <- function( m ) { # Convert igraph object to sparse Matrix
  return( igraph::get.adjacency( m, attr = 'weight' ) )
}
getdensity <- function( m ) { # Return matrix density.
  #	return( round( length( na.omit( as.vector( m ) ) ) / ( nrow( m ) * ncol( m ) ), 4 ) )
  N <- ncol( m )
  denom <- ( N * ( N - 1 ) ) / 2	
  edges <- as.vector( m[ upper.tri( m ) ] )
  total.edges <- length( stats::na.omit( edges ) )
  density <- total.edges / denom
  return( round( density, 4 ) )
}
###
# Weighted clustering coefficient undirected graphs, Onnela 2005.
# The weighted clustering coefficient is the average "intensity" of triangles around a node.
# W: input is undirected weighted matrix.
# Example:
#   W <- rbind( c( 0,2,3 ), c( 2,0,2 ), c( 3,2,0 ) )
#   clustering_onnela( W )
##
clustering_onnela <- function( g ){
  if( is.null( igraph::E( g )$weight ) ) 
    stop( "no weights in graph!" )
  degree <- igraph::degree( g ) 
  w <- apply( igraph::get.adjacency( g, attr = 'weight' ), 1, as.numeric ) # if edges are 'chars' 
  w[ is.na( w ) ] <- 0 # replace NA with 0	
  t <- w^( 1 / 3 ) 
  cyc3 <- diag( t %*% t %*% t ) 
  degree[ cyc3 == 0 ] <- Inf              # if no 3-cycles exist, make C=0 (via K = Inf)
  cc <- cyc3 / ( degree * ( degree - 1 ) ) # clustering coefficient
  out <- cbind( clustering = cc )
  rownames( out ) <- igraph::V( g )$name  
  return( out )
}
strength <- function( g )  { # Return strength of weighted undirected graph.
  if( is.null( igraph::E( g )$weight ) ) 
    stop( "no weights in graph!" )
  w <- apply( igraph::get.adjacency( g, attr = 'weight' ), 1, as.numeric ) # if edges are 'chars' 
  strength <- rowSums( w, na.rm = TRUE )
  out <- cbind( strength = strength )
  rownames( out ) <- igraph::V( g )$name
  return( out )
}
pathlength <- function( g, retmat = FALSE ) { # Shortest path length (harmonic mean)
  if( is.null( igraph::E( g )$weight ) )
    stop( "no weights in graph!" )
  # invert weights
  igraph::E( g )$weight <- 1 / ( 1 / as.numeric( igraph::E( g )$weight ) )
  paths <- igraph::shortest.paths( g, mode = "out", weights = igraph::E( g )$weight )
  diag( paths ) <- NA
  out <- cbind( paths = 1 / rowMeans( paths, na.rm = TRUE ) )
  rownames( out ) <- igraph::V( g )$name
  if( retmat )
    return( paths )
  else
    return( out )
}

###
# Remain degree distribution of randomized network.
# d => weighted or binary, undirected network (as matrix, missing edges as NA or 0.0). 
# niter => randomization iterations per existing edge.
# technical details: http://www.ncbi.nlm.nih.gov/pubmed/11988575
##
randomize <- function( d, niter = 3 ) {
  ###
  # Internal function: make symmetric matrix.
  ##
  make_symmetric <- function( m ) {
    m[ is.na( m ) ] <- 0
    m <- m + t( m )
    m[ m == 0 ] <- NA
    return( m )
  }
  
  # check input size
  if( nrow( d ) != ncol( d ) )
    stop( "*** ERROR ***: input is not symmetrical. Input should be a undirected matrix" )
  # remove lower corner and diagonal (->NA)
  d[ is.na( d ) ] <- 0
  d[ lower.tri( d ) ] <- NA
  diag( d ) <- NA
  # check if matrix is not sparse.
  print(d)
  print(any(d[ upper.tri( d ) ] ==0))
  if( ! any( d[ upper.tri( d ) ] == 0 ) )
    stop( "*** ERROR ***: input network is not sparse. Edge shuffle is not possible." )
  # determine edges	
  K <- nrow( data.frame( r = row( d )[ which( !d == 0 ) ], c = col( d )[ which( !d == 0 ) ] ) ) 
  # maximal number of rewiring attempts per 'iter'
  niter <- K * niter
  n <- nrow( d )
  maxAttempts <- round( n * K / ( n * ( n - 1 ) ) )
  suc <- 0
  for( i in 1:niter )  {
    # indices of edges with and without weights
    idx.full <- data.frame( r = row( d )[ which( !d == 0 ) ], c = col( d )[ which( !d == 0 ) ] )
    e <- idx.full[ sample.int( nrow( idx.full ), 2 ), ] 
    ab <- d[ e[1, 'r'], e[1, 'c'] ]
    cd <- d[ e[2, 'r'], e[2, 'c'] ]
    ad <- d[ e[2, 'r'], e[1, 'c'] ]
    cb <- d[ e[1, 'r'], e[2, 'c'] ]
    if( ( is.na( ad ) || is.na( cb ) ) || ( ( ad > 0 ) || ( cb > 0 ) ) )  {
      #print( paste( "skip: ", ad, cb ) )
    } else {
      # set ab, cd to 0
      d[ e[1, 'r'], e[1, 'c'] ] <- 0
      d[ e[2, 'r'], e[2, 'c'] ] <- 0
      # set ad => ab, cb => cd
      d[ e[2, 'r'], e[1, 'c'] ] <- ab
      d[ e[1, 'r'], e[2, 'c'] ] <- cd
      suc <- suc + 1
    }
  }
  print( paste( "Acceptance: ", ( suc / niter ) * 100, "%", sep = '' ) )
  # set back to full matrix
  return( make_symmetric( d ) )
}
make.symmetric <- function( mat ) { # average two triangles of matrix to convert to symmetrical version.
  # in two gamma matrices very small differences between lower and upper triangle, so average out.
  if( mean( mat[ upper.tri(mat) ], na.rm = T ) != mean( mat[ lower.tri(mat) ], na.rm = T ) ) mat <- ( mat + t( mat ) ) / 2 
  return( mat )
}
graphproperties <- function( g, mat ) { # Return list with graph properties
  out <- list()
  out[[ 'd' ]] <- mat.degree( mat )
  out[[ 's' ]] <- mat.strength( mat )
  out[[ 'c' ]] <- clustering_onnela( g )
  out[[ 'l' ]] <- pathlength( g )
  return( as.data.frame( out ) )
}
mat.strength <- function( m ) { # Return strength of Matrix
  return( ( rowSums( m, na.rm = TRUE ) ) / 2 )
}
mat.degree <- function( m ) { # Return degree of Matrix
  return( ( rowSums( binarize( m ), na.rm = TRUE ) ) / 2 )
}
binarize <- function(m) { # Binarize Matrix
  m[] <- ifelse( m > 0, 1, 0 ) 
  return( m )
}
realproperties <- function( mat ) { # Return real properties
  realg <- matrix2graph( make.symmetric( mat ) )
  return( graphproperties( realg, mat ) )
}
randomproperties <- function( mat, nrandom = 10 ) {
  rout <- list()
  symmat = make.symmetric(mat)
  for( i in 1:nrandom ) {
    rsymmat = randomize(symmat)
    rg <- matrix2graph(rsymmat)
    rout[[ i ]] <- graphproperties( rg, mat )
  }
  return( rout )
}
getgraph <- function( mat ) { # matrix -> graph
  g <- igraph::graph.adjacency( mat, mode = c( "undirected" ), weighted = TRUE, diag = FALSE )
  return( g )
}
getmst <- function( g ) { # mst graph
  # get minimum spanning tree using Prim's algorithm (based on weights)	
  mst.g <- igraph::minimum.spanning.tree( g, algorithm = 'prim' )
  # binarize minimum spanning tree
  igraph::E( mst.g )$weight <- 1
  return( mst.g )
}