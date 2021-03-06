phaselagindex = function(EEGdata,frequency=128) {
  inputw2 <- function(wave, f, channel=1)  { # function by w.m.otte@umcutrecht.nl (Wim Otte) 11-02-2015
    #renamed from inputw to inputw2 because there already is a inputw function in seewave
    if(is.data.frame(wave))   {f<-f ; wave <- as.matrix(wave[,channel])}
    if(is.vector(wave))       {f<-f ; wave <- as.matrix(wave)}
    # mts objects are matrix by default, there is then a conflict between is.matrix and is.mts
    if(is.matrix(wave) && !stats::is.mts(wave)) {f<-f ; wave <- wave[,channel,drop=FALSE]}  
    if(is.ts(wave))           {f<-frequency(wave) ; wave <- as.matrix(wave)} 
    if(stats::is.mts(wave))          {f<-frequency(wave) ; wave <- as.matrix(wave[, channel])} 
    if(class(wave)=="Sample") {f<-wave$rate ; wave <- as.matrix(wave$sound[channel, ])}
    if(class(wave)=="audioSample"){f<-wave$rate ; wave <- as.matrix(wave)}
    if(class(wave)=="Wave") {    
      f <- wave@samp.rate
      if(channel==1) {wave <- as.matrix(wave@left)}   
      if(channel==2) {wave <- as.matrix(wave@right)}     
    }    
    return(list(w=wave,f=f))
  } 
  hilbert <- function( wave, f )  { # function by w.m.otte@umcutrecht.nl (Wim Otte) 11-02-2015
    wave<- inputw2(wave=wave,f=f)$w
    n<- nrow(wave)
    ff<-fft(wave)
    h<-rep(0,n)
    if(n>0 & 2*floor(n/2)==n){h[c(1, n/2+1)]<-1; h[2:(n/2)]<-2}
    else{if(n>0){h[1]<-1; h[2:((n+1)/2)]<-2}}
    ht<-fft(ff*h,inverse=TRUE)/length(ff)
    return(ht)
  } 
  pli <- function( chan1, chan2, f ) { # function by w.m.otte@umcutrecht.nl (Wim Otte) 11-02-2015
    return( abs( mean( sign( Arg( hilbert( chan1, f = f ) ) - Arg( hilbert( chan2, f = f ) ) ) ) ) )
  }
  nchan = ncol(EEGdata)
  plimatrix= matrix(0,nchan,nchan)  
  for (i in 1:nchan) {
    for (j in 1:nchan) {
      plimatrix[i,j] = pli(EEGdata[,i], EEGdata[,j], f =frequency)
    }
  }
  
  output = list(meanpli=mean(plimatrix),sdpli=sd(plimatrix),ninetyppli=as.numeric(quantile(abs(plimatrix),probs=0.9)))
  #==========================
  # extra pli features
  # check
  if( mean( plimatrix, na.rm = T ) != 0 )  {
    plimatrix2 = plimatrix
    diag(plimatrix2) = 0
    plimatrix2[is.na(plimatrix2)] = 0
    
    g <- getgraph(plimatrix2) # get weighted graph
    mst.g <- getmst( g ) # get binary mst
    max.bc <- max( igraph::betweenness( mst.g ) )
    max.degree <- max( igraph::degree( mst.g ) ) # get max degree
    max.strength <- max( igraph::strength( mst.g ) ) # get max strength
    leafnumber <- sum( igraph::degree( mst.g ) == 1 ) # get leaf-number (e.g. nodes with degree 1 )
    m <- length( igraph::V( mst.g ) ) - 1 # get m (links)
    diameter <- m - leafnumber + 2 # get diameter
    ecc <- mean( igraph::eccentricity( mst.g ) ) # mean eccentricity
    radius <- igraph::radius( mst.g ) # smallest eccentricity
    Th <- leafnumber / ( 2 * m * max.bc ) # tree hierarchy
    kappa <- mean( igraph::degree( mst.g )^2 ) / mean( igraph::degree( mst.g ) ) # kappa
    density <- getdensity( plimatrix )  # get density
    d <- data.frame(densitypli = density, maxdegreepli = max.degree, 
                    maxstrengthpli = max.strength, mpli = m, diameterpli = diameter, leafnumberpli = leafnumber, maxbcpli = max.bc,
                    eccpli = ecc, radiuspli = radius, Thpli = Th, kappapli = kappa )    
    # final <- rbind( final, d )
    dlist <- as.list(d)
    output = append(output,dlist)
    invisible(output)
  } else {
    print( paste( "*** PROBLEM with:  phaselagindex calculation" ) )
    invisible(output)
  }
  #==========================
  # invisible(list(meanpli=mean(plimatrix),sdpli=sd(plimatrix),ninetyppli=as.numeric(quantile(abs(plimatrix),probs=0.9))))
}