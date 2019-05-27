
# ===========================================================
#  FUNCTIONS FOR RAIN GAUGE CALIBRATION AND PROCESSING
#  Author: Martin Fencl, CVUT (2019)

# Functions for processing rain gauge data designed specifically for
# Rain gauges and weather stations manufactured by  FIEDLER AMS s.r.o. 
# (https://www.fiedler.company/):
# Includes:
# - evaluation of laboratory (dynamic) calibraiton of rain gauges
# - dynamic calibration of rain gauge data
# - reading of logged data (in different formates)
# - data cleaning and processing
# - analysis of rain gauge data



read.scale <- function(path, fil.n, m.date){
  ## function to read and reformat data from laboratory scale
  ##
  ## inputs:  path   - path to the file
  ##          fil.n  - file name
  ##          m.date - date of measurement ("YYYY-mm-dd")
  ## outputs: data.frame with time [POSIXct] and weight [g]

  dat <- read.table(paste(path,"/" ,fil.n, sep=""))
  tim <- as.POSIXct(paste(m.date, dat[ ,1]), tz="UTC", format="%Y-%m-%d %H:%M:%S")
  return(data.frame("time" = tim, "weight" = as.numeric(dat[,2])))
}

###################


read_RG_formatV2 <- function(path, fil.n, m.date){
  ## function to read and reformat data from Fiedler RGs (from MOST, format v2)
  ##
  ## inputs:  path   - path to the file
  ##          fil.n  - file name
  ##          m.date - date of measurement ("YYYY-mm-dd")
  ##
  ## outputs: data.frame with time stamps of tips and number of tips

  dat <- read.table(paste(path,"/" ,fil.n, sep=""), skip=7)
  id.filt <- which(dat$V3 =="???")
  tim <- paste(m.date, dat[-id.filt, 2], sep=" ")
  tim <- as.POSIXct(tim, tz="UTC", format="%Y-%m-%d %H:%M:%S")
  tip <- as.numeric(as.character(dat[-id.filt, 3]))
  return(data.frame("time"=tim, "tip"=tip))
  
}


###################


read_fiedler_dta <- function(file_name) {
  # Read dta fiels from Fiedler weather station and format and convert them
  # to be in POSIXct time with UTC time zone.
  # Arguments:
  # file_name - file name (and path)
  # Value:
  # output - list with data (data.frame) and metadata
    
  dat <- readLines(file_name)
  rowid_station <- which(dat == '<station>')
  station_nams <- strsplit(dat[rowid_station + 1], '\t')[[1]]
  station_nams <- unlist(strsplit(station_nams, '\"'))
  station_nams <- station_nams[nchar(station_nams) > 0]
  station_vals <- strsplit(dat[rowid_station + 2], '\t')[[1]]
  station_vals <- as.data.frame(matrix(station_vals, 1, length(station_nams),
                                       byrow = T, dimnames = list(NULL, station_nams)),
                                stringsAsFactors = F)
  
  
  rowid_channels <- which(dat == '<channels>')
  channels_nams <- strsplit(dat[rowid_channels + 1], '\t')[[1]]
  channels_nams <- unlist(strsplit(channels_nams, '\"'))
  channels_nams <- channels_nams[nchar(channels_nams) > 0]
  rowid_data <- which(dat == '<data>')
  channels_vals <- strsplit(dat[(rowid_channels + 2):(rowid_data - 2)], '\t')
  channels_vals <- as.data.frame(matrix(unlist(channels_vals), length(channels_vals),
                                        length(channels_nams), byrow = T,
                                        dimnames = list(NULL, channels_nams)),
                                 stringsAsFactors = F)
  
  
  rowid_binaries <- which(dat == '<binaries>')
  data_vals <- strsplit(dat[(rowid_data + 2):(rowid_binaries - 2)], '\t')
  data_vals <- as.data.frame(matrix(unlist(data_vals), length(data_vals),
                                    nrow(channels_vals) + 1, byrow = T,
                                    dimnames = list(NULL, c('Time', channels_vals$Label))),
                             stringsAsFactors = F, dec = ',')
  for(i in 2 : ncol(data_vals)) {
    data_vals[, i] <- suppressWarnings(as.numeric(gsub(',', '.', data_vals[, i])))
  }
  
  data_vals[, 1] <- strptime(data_vals[ , 1], tz = 'UTC', format = "%d.%m.%Y %H:%M:%S") - as.numeric(station_vals$UTC_Offset)
                    
  
  
  output <- list('data' = data_vals, 'metadata' = list('station' = station_vals,
                                                       'data' = channels_vals))
  return(output)

}  
  

####################

calibrate_rg <- function (tim, tab) {
    # function to calibrate RG based on calibration table 
    # inputs:
    # tim - pulse times
    # tab - calibration table (lenght of pulse; volume of tip)
    # Outputs: vol - volume of tip (mm)
    
    l.puls <- tim[-1] - tim[-length(tim)]
    l.puls <- c(tab[1, 1], l.puls)
    l.puls[which(l.puls > tab[1, 1])] <- tab[1, 1]        #pulses longer than the longest interval in calibration table take as having the same length as this interval
    l.puls[which(l.puls < min(tab[, 1]))] <- min(tab[, 1])        #pulses longer than the shortes interval in calibration table take as having the same length as this interval
    
    #apply calibration table
    vol <- round(approx(tab[,1], tab[,2], l.puls, method  ="linear")$y, 3)
    
    return (vol)
}

calibrate_rg_zoo <- function (rg, tab) {
    # function for applying RG calibration on zoo object
    require(zoo)
    vol <- calibrate_rg(index(rg), tab)
    rg[] <- vol
    return (rg)
}


####################

RGtip.aggreg.by <- function(x, tim, step){
    ## function to aggregate time step of a time series of RG tips
    ##
    ## inputs
    ## x    - values of time series
    ## tim  - time indexes
    ## step - time step to which agregate the series (in minutes)
    ##
    ## outputs:
    ## y - aggregated time series
    
    t.num <- as.numeric(tim)
    t.num <- round(t.num/(step*60), 0)*step*60
    
    tim.agr <-  as.POSIXct(t.num, origin="1970-01-01 00:00:00 UTC", tz="UTC") #time indexes for aggregtion
    
    y <- aggregate(x, list(tim.agr), sum, na.rm=T)
    y$x <- y$x/step
    
    return(y)
    
}

RGtip.aggreg.by.zoo <- function (z, step) {
    require(zoo)
    y <- RGtip.aggreg.by(coredata(z), index(z), step)
    z2 <- zoo(y[ ,2], y[ ,1])
    return(z2)
}

#####################

reg.series <- function(dat, step){
    ## function to insert time steps to make time series regular
    ##
    ## inputs
    ## dat  - data.frame of n columns with time stamp in the first column
    ## step - time step which should have the resulting time series (in minutes)
    ##
    ## outputs:
    ## dat2 - data.frame with consistent time series
    tim0 <- dat[,1]
    
    tim <- seq(tim0[1],tim0[length(tim0)], by=step*60)
    id.add <- which(is.element(tim, tim0)==T)
    dat2 <- as.data.frame(matrix(NA, ncol=ncol(dat), nrow=length(tim)))
    dat2[ ,1] <- tim
    dat2[id.add, -1] <- dat[ ,-1]
    
    # name colnames according to original data.frame
    colnames(dat2) <- colnames(dat)
    
    return(dat2)
}


reg.series.zoo <- function (ser, step) {
    require(zoo)
    
    dat <- data.frame(index(ser), coredata(ser))
    print(head(dat))
    dat2 <- reg.series(dat, step)
    
    ser2 <- zoo(dat2[ , -1], dat[ , 1])
    return (ser2)
}
#####################


reg.series2 <- function(dat, step, int){
    ## function to insert time steps to make time series regular
    ##
    ## inputs
    ## dat  - data.frame of n columns with time stamp in the first column
    ## step - time step which should have the resulting time series (in minutes)
    ## int  - time vector with beggining and end of the series
    ##
    ## outputs:
    ## dat2 - data.frame with consistent time series
    tim0 <- dat[ ,1]
    
    
    # if series beggins before defined interval, cut it
    if(min(tim0, na.rm=T) < int[1]){  
        dat <- dat[-which(tim0 < int[1]), ]
        tim0 <- dat[ ,1]
    }
    
    # if series ends after defined interval, cut it
    if(max(tim0, na.rm=T) > int[2]){  
        dat <- dat[-which(tim0 > int[2]), ]
        tim0 <- dat[ ,1]
    }
    
    # add time stamps with NA values to create regular series
    tim <- seq(int[1],int[2], by=step*60)
    id.add <- which(is.element(tim, tim0)==T)
    dat2 <- as.data.frame(matrix(NA, ncol=ncol(dat), nrow=length(tim)))
    dat2[ ,1] <- tim
    dat2[id.add, -1] <- dat[ ,-1]
    
    # name colnames according to original data.frame
    colnames(dat2) <- colnames(dat)
    
    return(dat2)
}

######################


RG_smooth <- function(y, gap){
    ## function to create smooth time series from RG rainfall
    ## y   - RG values
    ## gap - number of time steps before tip when 0 is asumed to be zero rainfall.
    ##
    ##
    ##
    onetip <- min(y[which(y!=0)], na.rm=T)
    
    x <- y
    x[which(x > onetip)] <- onetip
    id <- which(x > 0)
    n.zer <- c(id[1], id[-1] - id[-length(id)])
    stamp <- rep(id, n.zer)
    #data.frame(stamp, x)
    #x.agr <- aggregate(x, by=list(stamp), mean)
    y.mean <- y - x
    for(i in 1:length(id)){
        id2 <- which(stamp == id[i])
        rain.st <- 3
        if(i != length(id)){rain.st <- length(which(stamp == id[i+1]))}
        if(length(id2) <= gap){
            y.mean[id2] <- y.mean[id2] + mean(x[id2]) 
        }else{
            id3 <- id2[length(id2)-c(2,1,0)]
            y.mean[id3] <- y.mean[id3] + mean(x[id3]) 
        }
    }
    
    return(y.mean)
}


###############
zoo_aggreg_by <- function(x_zoo, step, fun, align = 'center',
                          insert.missing = T, ...){
    ## function to aggregate time series to given time series
    ##
    ## inputs:
    ## x_zoo - zoo time series
    ## step - time step to which agregate the series (in minutes)
    ## step - time step to which agregate the series (in minutes)
    ##
    ## outputs:
    ## ag_zoo - aggregated time series
    
    require(zoo)
    tim <- index(x_zoo)
    
    t_num <- as.numeric(tim)
    
    if(align == 'left'){
        t_num <- floor(t_num/(step*60))*step*60    
    }else if(align == 'center'){
        t_num <- round(t_num/(step*60), 0)*step*60
    }else if(align == 'right'){
        t_num <- ceiling(t_num/(step*60))*step*60
    }else{stop('Error in match.arg(align) : \n 
                arg should be one of center, left, right)')}
    
    
    tim_agr <-  as.POSIXct(t_num, origin="1970-01-01 00:00:00 UTC") #time indexes for aggregtion
    
    ag_zoo <- aggregate(x_zoo, list(tim_agr), fun)
    
    if(is.regular(ag_zoo, strict = T) == F) {
        if(insert.missing == T){
            ag_zoo <- insert_missing_records(ag_zoo, step)
            print('missing time steps where inserted as NA values')  
        } else {
            warning("new time series is not strictly regular!")
        }
    }
    
    
    return(ag_zoo)
    
}

#################


insert_missing_records <- function(zoo_ser, step){
    ## function to insert time steps to make time series regular
    ##
    ## inputs
    ## dat  - data.frame of n columns with time stamp in the first column
    ## step - time step which should have the resulting time series (in minutes)
    ##
    ## outputs:
    ## dat2 - data.frame with consistent time series
    
    require(zoo)
    
    dat <- data.frame('time' = index(zoo_ser), 'data'= coredata((zoo_ser)))
    
    dat2 <- reg_series(dat, step)
    zoo2_ser <- zoo(dat2[ , -1], dat2[ , 1])
    colnames(zoo2_ser) <- colnames(zoo_ser)
    
    return(zoo2_ser)
}

###################


reg_series <- function(dat, step){
    ## function to insert time steps to make time series regular
    ##
    ## inputs
    ## dat  - data.frame of n columns with time stamp in the first column
    ## step - time step which should have the resulting time series (in minutes)
    ##
    ## outputs:
    ## dat2 - data.frame with consistent time series
    tim0 <- dat[,1]
    tim <- seq(tim0[1],tim0[length(tim0)], by=step*60)
    id.add <- which(is.element(tim, tim0)==T)
    dat2 <- as.data.frame(matrix(NA, ncol=ncol(dat), nrow=length(tim)))
    dat2[ ,1] <- tim
    dat2[id.add, -1] <- dat[ ,-1]
    
    # name colnames according to original data.frame
    colnames(dat2) <- colnames(dat)
    
    return(dat2)
}



# =============================

# Rain gauge statistics

identify_Revents <- function(R, win.max = 30, min.len=10, NAs = 'na.or.complete') {
  
  ## function to identify rainy periods from RG series of regular time step
  ## Inputs:  
  ##          R -  time series with rain rates from one (vector) or more (matrix)
  ##               rain gauges
  ##          win.max - maximum size of dry window between two wet time steps
  ##                   to assume them to belong to same rain event
  ##          min.len - minimum length of period [minutes] to assume it as event
  ##          NAs - 'pass', 'remove', 'pass.or.remove', should NAs be passed,
  ##                 removed, or removed unless there are only NAs.
  ## Outputs:  data frame with begginings and ends of rainfall periods
  
  #argument test and parameter replication
  
  tim <- index(R)
  mtx <- coredata(R)
  mtx <- as.matrix(mtx) 
  
  if(missing(win.max)==T){win.max <- 30}
  if(missing(min.len)==T){min.len <- 10}

  #drywet vecotr
  nas.count <- apply(apply(mtx, 2, is.na), 1, sum)
  if (NAs %in%  c('pass', 'rm', 'pass.or.rm') == F){
    stop('Argument NAs has to be \'pass\', \'rm\', \'pass.or.rm\'')
  }

  
  if (NAs == 'pass') {
    drywet0 <- apply(mtx!=0, 1, sum) != 0
  }
  
  if (NAs == 'rm') {
    if (sum(nas.count == ncol(mtx)) > 0) {
      warning('one or more rows with no non-NA values were classified as dry')
    }
    drywet0 <- apply(mtx!=0, 1, sum, na.rm = T) != 0
  }
  
  if (NAs == 'pass.or.rm') {
    drywet0 <- apply(mtx!=0, 1, sum, na.rm = T)
    drywet0[nas.count > nrow(mtx)] <- NA
    drywet0 <- drywet0 != 0
  }
  
  #find wet weather time stamps and time difference between them
  wet.tim <- tim[which(drywet0 > 0)]
  wet.int <- as.numeric(wet.tim[-1]) - as.numeric(wet.tim[-length(wet.tim)])
  
  #find begginings and ends of events
  ev.st <- c(wet.tim[1], wet.tim[which(wet.int > win.max*60)+1])
  ev.en <- c(wet.tim[which(wet.int > win.max*60)], wet.tim[length(wet.tim)])
  
  event.times <- data.frame("st" = ev.st, "en"=ev.en)
  event.times <- event.times[-which(difftime(ev.en, ev.st, "mins") < min.len), ] 
  
  return(event.times)
}




summarize_singleRevent <- function (R, na.rm = T) {
  ## function to provide summary statistics of rainfall events 
  ##
  ## Inputs:  R  -  time series with rain rate from a single instrument
  ##
  ## Outputs: res - vector with basic rainfall statistics
  
  res <- c("duration" = NA, "height" = NA, "Rmax" = NA, "Rmax10" = NA)


  if (length(which(!is.na(R))) == 0){
    return(res)
  } else {

    res[1] <- difftime(end(R), start(R), units="mins") #duration of rainfall
    res[2] <- sum(R, na.rm = na.rm)/60    #total height of rainfall [mm]
    res[3] <- max(R, na.rm = na.rm)       #max rain rate [mm/h]
    
    if(length(R) < 10){               #max 10min rain rate [mm/h]
      res[4] <- sum(R, na.rm = na.rm) / 10
    }else{
      res[4] <- max(rollapply(R, 10, mean, na.rm = na.rm), na.rm = na.rm) 
    }
  }
  
  return(res)
}  



summarize_Revents <- function (R, st = NULL, en = NULL, na.rm = T) {
  ## function to provide summary statistics of rainfall events 
  ##
  ## Inputs:  R  -  time series with rain rate from a single instrument
  ##          st  -  vector with times where the interval begins
  ##          en  -  vector with times where the interval ends
  ##
  ## Outputs: tab - table
  
  if (is.null(st)) {st <- start(R)}
  if (is.null(en)) {en <- end(R)}
  if (length(st) != length(en)) {stop('st and en vectors do not have the same length!')}

  tab <- data.frame("start" = st, "end" = en, "duration" = NA, "height" = NA,
                    "Rmax" = NA, "Rmax10" = NA)
  tab[ ,1] <- st
  tab[ ,2] <- en
  
  for (i in 1 : length(st)) {
    tab[i, 3 : 6] <- summarize_singleRevent(window(R, start = st[i], end = en[i]), na.rm = na.rm)
  }
  
  return(tab)
}  






