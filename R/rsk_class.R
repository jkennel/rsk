#' @export
Rsk <- R6Class("Rsk",
               public = list(

  # database name
  file_name     = NA_character_,


  # tables from rsk
  blob         = list(),
  channels     = data.table(),
  errors       = data.table(),
  instruments  = data.table(),
  events       = data.table(),
  coefficients = data.table(),
  continuous   = data.table(),


  # table of measurements
  data         = data.table(),


  time_1        = as.POSIXct(NA_real_),
  time_interval = NA_real_,

  names         = c(),
  units         = NA_character_,

  pres_index    = NA_integer_,
  temp_index    = NA_integer_,

  n_obs         = NA_integer_,
  n_chan        = NA_integer_,

  serial        = NA_character_,
  model         = NA_character_,

  latitude      = NA_real_,
  longitude     = NA_real_,
  elevation     = NA_real_,

  time_start    = NULL,
  time_end      = NULL,
  has_gaps      = NULL,



  initialize = function(file_name, subset = FALSE) {

    self$file_name <- file_name

    # get data base tables
    db <- dbConnect(SQLite(), file_name)
    # dbListTables(db)
    self$blob         <- setDT(dbGetQuery(db, "SELECT * FROM downloads"))
    self$coefficients <- setDT(dbGetQuery(db, "SELECT calibrationID, key, cast(value as REAL) as value FROM coefficients"))
    self$channels     <- setDT(dbGetQuery(db, "SELECT * FROM channels"))
    self$instruments  <- setDT(dbGetQuery(db, "SELECT * FROM instruments"))
    self$errors       <- setDT(dbGetQuery(db, "SELECT * FROM errors"))
    self$events       <- setDT(dbGetQuery(db, "SELECT * FROM events"))
    self$continuous   <- setDT(dbGetQuery(db, "SELECT * FROM continuous"))

    dbDisconnect(db)


    # simple constants
    self$n_chan        <- nrow(self$channels)
    self$time_1        <- self$events$tstamp[[1]] / 1000 - self$rbr_time_reference()
    self$time_interval <- self$continuous[["samplingPeriod"]] / 1000
    self$model         <- self$instruments[["model"]]
    self$serial        <- self$instruments[["serialID"]]


    co_compensation    <- self$comp_coefficients()
    self$channel_names()
    self$channel_units()


    # get raw values
    raw_val <- unlist(lapply(self$blob$data, function(x) {
      readBin(x,
              n = 100000L,
              what = 'integer',
              size = 4L,
              signed = TRUE,
              endian = 'little')
    }))


    # determine non representative values
    h_len  <- max(which(raw_val[1:300] == self$time_1))
    to_rem <- c(1:h_len, (self[['events']][['sampleIndex']][-1]-1) * self$n_chan + h_len + 1:2)


    # calculate mv, dbar, temperature
    self$data <- setDT(
      test_all(
        raw_val,
        to_rem - 1, # c++ base 0 indexing
        self$n_chan,
        as.integer(self[['events']][['tstamp']] * 1e-3),
        self$events[["sampleIndex"]],
        self$time_interval,
        self$base_coefficients(),
        self$is_temperature(),
        co_compensation,
        self$pres_index * 2L,
        self$temp_index * 2L,
        self$names
      ), key = 'datetime')

    self$n_obs <- nrow(self$data)
    self$update_time_ranges()

    invisible(self)

  },
  update_time_ranges = function() {

    self$time_start <- self$data[1]$datetime
    self$time_end   <- self$data[.N]$datetime
    self$has_gaps   <- ((as.numeric(self$time_end)-
      as.numeric(self$time_start)) /
      self$time_interval + 1) != self$n_obs


    invisible(self)
  },

  channel_names = function() {

    base_name <- tolower(self$channels[["longNamePlainText"]])

    if(self$temp_index != 0L & self$pres_index != 0) {

      base_name[self$temp_index] <- paste0(base_name[self$temp_index], '_onboard')
      base_name  <- rep(base_name, each = 2)
      base_name  <- paste0(base_name, rep(c('_raw', ''), times = self$n_chan))
      base_name  <- c(base_name, 'pressure_dbar_comp')
      self$names <- base_name

    } else {

      base_name  <- rep(base_name, each = 2)
      base_name  <- paste0(base_name, rep(c('_raw', ''), times = self$n_chan))
      self$names <- base_name

    }

    invisible(self)

  },

  channel_units = function() {

    units <- rep(NA_character_, length(self$names))
    units[grepl('pressure', self$names)] <- 'dbar'
    units[grepl('temperature', self$names)] <- 'degrees_c'
    units[grepl('_raw', self$names)] <- 'mv'

    self$units <- units

    invisible(self)
  },

  # which channels are temperature
  is_temperature = function() {
    self$channels[["longNamePlainText"]] == 'Temperature'
  },

  # return a matrix (or vector) of the calibration coefficients (i.e. c0,c1,...)
  base_coefficients = function(id = NULL) {

    if(is.null(id)) {

      co <- matrix(
        self$coefficients[grepl('c', key)][['value']],
        ncol = self$n_chan)

      return(co)

    } else {

      return(self$coefficients[grepl('c', key) & calibrationID == id][['value']])
    }
  },

  # return a vector of the compensation coefficients (i.e. x0,x1,...)
  comp_coefficients = function() {
    comp_cols <- self$coefficients[key == 'n0']

    if(nrow(comp_cols) > 0) {

      self$pres_index <- comp_cols$calibrationID
      self$temp_index <- comp_cols$value

      return(self$coefficients[grepl('x', key)]$value)

    }

    self$pres_index <- 0L
    self$temp_index <-0L

    return(NA_real_)

  },

  glance = function() {
    if(self$n_obs < 1e5){
      print(plot_ly(self$data,
                    x = ~datetime,
                    y = ~value,
                    color = ~variable,
                    type = 'scatter',
                    mode = 'lines'))
    } else {

      data_sub <-  self$data[as.integer(seq.int(1, .N, length.out = 50000))]
      data_sub <- melt(data_sub, id.vars = 'datetime')

      print(plot_ly(data_sub,
                    x = ~datetime,
                    y = ~value,
                    color = ~variable,
                    type = 'scatter',
                    mode = 'lines'))

    }

    invisible(self)
  },


  rbr_time_reference = function() {
    as.numeric(as.POSIXct('2000-01-01 00:00:00', tz = 'UTC'))
  })


)









# library(rsk)
# library(R6)
# library(data.table)
# library(RSQLite)
#
# x <- '/home/jonathankennel/Storage/data/rbr/rd45a 081871_20191118_1213.rsk'
# x <- '/home/jonathankennel/Storage/data/rbr/rd130 077623_20191119_1500.rsk'
# x <- '/home/jonathankennel/Storage/data/rbr/rd45a 081871_20191118_1213.rsk'
# system.time(z <- Rsk$new(x))
# system.time(z <- Rsk$new(x)$data)
# system.time(zz <- transducer::read_rbr(x)$data[[1]])
# profvis::profvis(z <- Rsk$new(x))
# zz <- z$data
#
# plot(V2~datetime, zz[40:2000], type='l')
# points(V6~datetime, zz[40:2000], type='l', col = 'red')
# # str(z$base_coefficients(id = 1))
# # str(z$comp_coefficients(id = 2))

