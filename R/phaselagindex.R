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
  invisible(list(meanpli=mean(plimatrix),sdpli=sd(plimatrix),ninetyppli=as.numeric(quantile(abs(plimatrix),probs=0.9))))
}