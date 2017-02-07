# # Read edf file and convert it to csv
# rm(list=ls())
# graphics.off()
# library( 'edfReader' )
# 
# edf2header <- function( fname ) {
#   # Get header as csv from edf.
# 	S <- edfReader::readEdfSignals( H <- readEdfHeader( fname ), fragments = TRUE )
# 	samplerate <- S$AF3$sRate
# 	header <- data.frame( 
# 		fname = fname,	
# 		continuousRecording = S$AF3$isContinuous,
# 		recordedPeriod =  S$AF3$recordedPeriod,
# 		startTime = S$AF3$startTime,
# 		from = S$AF3$from,
# 		start = S$AF3$start,
# 		totalPeriod = S$AF3$totalPeriod,
# 		till = S$AF3$till,
# 		transducerType = S$AF3$transducerType,
# 		sampleRate = S$AF3$sRate,
# 		preFilter = S$AF3$preFilter,
# 		bitsPerSample = S$AF3$sampleBits,
# 		range = S$AF3$range )
# 	return( header )
# }
# edf2signal <- function( fname ) {
#   # Get signals as csv from edf.
# 	#fname = "data/emotivexample/example-data-control--26.05.2016.22.37.50.edf"
# 	S <-  edfReader::readEdfSignals( H <- readEdfHeader( fname ), fragments = TRUE )
# 	store <- data.frame( 
# 				AF3 = S$AF3$signal, AF4 = S$AF4$signal, F3 = S$F3$signal, F4 = S$F4$signal,
# 				F7 = S$F7$signal, F8 = S$F8$signal, FC5 = S$FC5$signal, FC6 = S$FC6$signal,
# 				O1 = S$O1$signal, O2 = S$O2$signal,	P7 = S$P7$signal, P8 = S$P8$signal,
# 				T7 = S$T7$signal, T8 = S$T8$signal, COUNTER = S$COUNTER$signal, INTERPOLATED = S$INTERPOLATED$signal,
# 				GYROX = S$GYROX$signal, GYROY = S$GYROY$signal, RAW_CQ = S$RAW_CQ$signal, CQ_CMS = S$CQ_CMS$signal,
# 				CQ_F7 = S$CQ_F7$signal, CQ_T7 = S$CQ_T7$signal, CQ_O2 = S$CQ_O2$signal,	 CQ_FC6 = S$CQ_FC6$signal,
# 				CQ_AF4 = S$CQ_AF4$signal, CQ_F3 = S$CQ_F3$signal, CQ_P7 = S$CQ_P7$signal, CQ_P8 = S$CQ_P8$signal,
# 				CQ_F4 = S$CQ_F4$signal, CQ_AF3 = S$CQ_AF3$signal, CQ_FC5 = S$CQ_FC5$signal, CQ_O1 = S$CQ_O1$signal,
# 				CQ_T8 = S$CQ_T8$signal, CQ_F8 = S$CQ_F8$signal, CQ_DRL = S$CQ_DRL$signal )
# 	return( store )
# }
# 
# indir <- 'data'
# outdir <- 'converted'
# dir.create( outdir, showWarnings = FALSE )
# 
# # select .edf files
# files <- dir( indir )
# files <- files[ grep( ".edf", files ) ]
# file.ids <- sapply( strsplit( as.character( files ), "-" ), "[", 1 )
# session.ids <- sapply( strsplit( as.character( files ), "-" ), "[", 2 )
# id <- paste( file.ids, '-', session.ids, sep = '' )
# 
# # merge in data.frame
# d <- data.frame( file = files, id = file.ids, session = session.ids )
# 
# headers <- NULL
# for( i in 1:nrow( d ) ) {
# 	file <- d[ i, 'file' ]
# 	id <- d[ i, 'id' ]
# 	session <- d[ i, 'session' ]
# 	fname <- paste( indir, '/', file, sep = '' )
# 	# write output if file does not exists
# 	outfile <- paste( outdir, "/signal-", id, "-", session, ".csv.gz", sep = '' )
# 	print( paste( file, session, id, outfile ) )
# 	if( ! file.exists( outfile ) )	{
# 		gz.file <- gzfile( outfile, "w" )
# 		header <- edf2header( fname )
# 		header$subject.id <- id
# 		header$session.id <- session
# 		headers <- rbind( headers, header )
# 		signal <- edf2signal( fname )
# 		write.csv( signal, gz.file )
# 		close( gz.file )
# 	}
# }
# dim( headers )
# write.csv( headers, file = paste( outdir, "/meta-info.csv", sep = '' ) )