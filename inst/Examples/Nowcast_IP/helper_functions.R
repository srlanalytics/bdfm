
#9c39f8a4fbcc38c22475bc7f26988367

#A lightweight function to call data from the fred API

get_fred_data <- function(series_id, api_key, series_name = NULL, observation_start = "1776-07-04", observation_end = "9999-12-31", output_type = 1, vintage_dates = NULL, realtime_start = NULL, realtime_end = NULL, frequency = NULL){
  fred_call <- paste0("https://api.stlouisfed.org/fred/series/observations?series_id=",
                      series_id,
                      "&observation_start=",
                      observation_start,
                      "&observation_end=",
                      observation_end,
                      "&observation_type=",
                      output_type)
  if(!is.null(vintage_dates)){
    fred_call <- paste0(fred_call, "&vintage_dates=",vintage_dates)
  }else{
    if(!is.null(realtime_start)){
      fred_call <- paste0(fred_call,"&realtime_start=",realtime_start)
    }
    if(!is.null(realtime_end)){
      fred_call <- paste0(fred_call,"&realtime_end=",realtime_end)
    }
  }#cannot use realtime_start and realtime_end with vintage dates
  if(!is.null(frequency)){
    fred_call <- paste0(fred_call,"&frequency=",frequency)
    #fred_call <- paste0(fred_call,"&aggregation_method=avg")
  }
  fred_call <- paste0(fred_call,
                      "&api_key=",api_key,
                      "&file_type=json")
  X   <- fromJSON(fred_call)$observations
  X[X$value == ".","value"] <- NA 
  X <- data.frame(realtime_start = as.Date(X$realtime_start), realtime_end = as.Date(X$realtime_end), date = as.Date(X$date), value = as.numeric(X$value))
  if(is.null(series_name)){
    series_name = series_id
  }
  X <- cbind.data.frame(rep(series_name,nrow(X)), X)
  colnames(X)[1] <- "series_name"
  return(X)
}

# DATA3 <- get_fred_data(series_id = 'INDPRO', api_key = "9c39f8a4fbcc38c22475bc7f26988367", observation_start = "1980-01-01", realtime_start = "2018-09-30", realtime_end = "2018-09-30", frequency = "q")

#A wrapper for dcast (part of data.table) to put our time series data into matrix format

call_data <- function(..., dt = Data, as_of = NULL, date = "date", series_id = "series_id", value = "value"){
  
  dt <- data.table(dt)
  
  ##### This section makes sure column names are consistent since data.table uses colnames for most  operations ####
  
  
  colnames(dt)[colnames(dt)==date]        <- "date"
  colnames(dt)[colnames(dt)==series_id]   <- "series_id"
  colnames(dt)[colnames(dt)==value]       <- "value"
  
  setkey(dt, series_id, date)
  
  # for convenience, so we can write call_data("dfsd", "dsfsdf")
  l <- list(...)
  is.char <- (sapply(l, class) == "character")
  if (any(!is.char)){
    stop("some elements are not character vectors")
  }
  oi <- do.call("c", l) #of interest
  
  unique_names   <- unique(dt[,series_id])
  missing_series <- oi[!oi%in%unique_names]
  
  if(length(missing_series)>0){
    warning(paste("The following series are missing:", paste(missing_series, collapse = ", ")))
  }
  
  Out <- dcast(dt[series_id%in%oi], date ~ series_id,    value.var = "value")
  
  return(Out)
}

