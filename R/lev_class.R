#' #' @export
#' Transducer <- R6Class("Transducer",
#'                public = list(
#'
#'                  # database name
#'                  file_name     = NA_character_,
#'                  file_extension = NA_character_,
#'
#'                  # tables from rsk
#'                  blob         = NULL,
#'                  channels     = NULL,
#'                  errors       = NULL,
#'                  instruments  = NULL,
#'                  events       = NULL,
#'                  coefficients = NULL,
#'                  continuous   = NULL,
#'
#'
#'                  # table of measurements
#'                  data          = NULL,
#'
#'
#'                  time_1        = NULL,
#'                  measurement_interval = NULL,
#'
#'                  names         = NULL,
#'                  units         = NULL,
#'
#'                  pressure_index    = 0,
#'                  temperature_index = 0,
#'                  is_temperature_compensated = FALSE,
#'
#'                  n_obs         = NULL,
#'                  n_chan        = NULL,
#'
#'                  serial        = NULL,
#'                  model         = NULL,
#'
#'                  latitude      = NULL,
#'                  longitude     = NULL,
#'                  elevation     = NULL,
#'
#'                  time_start    = NULL,
#'                  time_end      = NULL,
#'                  has_gaps      = NULL,
#'
#'                  db_tables     = NULL,
#'
#'
#'
#'                  #' @description
#'                  #' Class to hold data from an .rsk or .duckdb file
#'                  #'
#'                  #' @param file_name
#'                  #' @param subset
#'                  #' @param latitude
#'                  #' @param longitude
#'                  #' @param elevation
#'                  #' @param start
#'                  #' @param end
#'                  #' @param by
#'                  #' @param times
#'                  #' @param raw
#'                  #' @param simplify_names
#'                  #'
#'                  #' @return
#'                  #' @export
#'                  #'
#'                  #' @examples
#'                  initialize = function(file_name,
#'                                        subset = FALSE,
#'                                        latitude = NULL,
#'                                        longitude = NULL,
#'                                        elevation = NULL,
#'                                        start = NULL,
#'                                        end = NULL,
#'                                        by = NULL,
#'                                        times = NULL,
#'                                        raw = FALSE,
#'                                        keep_raw = FALSE,
#'                                        simplify_names = FALSE
#'                  ) {
#'
#'                    self$file_name <- file_name
#'                    self$file_extension <- tools::file_ext(file_name)
#'
#'                    if (self$file_extension == "rsk") {
#'
#'                    } else if (self$file_extension == "lev") {
#'
#'                    }
#'
#'                    # get data base tables
#'                    if (tools::file_ext(file_name) == "rsk") {
#'                      db <- dbConnect(SQLite(), file_name)
#'                      self$db_tables    <- RSQLite::dbListTables(db)
#'                      if ("downloads" %in% self$db_tables) {
#'                        self$blob <- (dbGetQuery(db, "SELECT * FROM downloads"))
#'                      } else {
#'                        if ("data" %in% self$db_tables) {
#'
#'                          keep_raw <- FALSE
#'                          raw <- FALSE
#'
#'                        } else {
#'                          dbDisconnect(db)
#'
#'                          return(self)
#'                        }
#'                      }
#'                    } else if (tools::file_ext(file_name) == "duckdb") {
#'                      if (raw) {
#'                        stop("raw cannot be TRUE if file is type .duckdb")
#'                      }
#'                      db <- dbConnect(duckdb(), file_name)
#'                    } else {
#'                      dbDisconnect(db)
#'                      stop("File type not currently supported")
#'                    }
#'                    self$db_tables    <- RSQLite::dbListTables(db)
#'
#'                    if ("coefficients" %in% self$db_tables) {
#'                      self$coefficients <- collapse::qDT(dbGetQuery(db, "SELECT calibrationID, key, cast(value as REAL) as value FROM coefficients"))
#'                    } else if ("calibrations" %in% self$db_tables) {
#'                      tmp <- collapse::qDT(dbGetQuery(db, "SELECT * FROM calibrations"))
#'                      nms <- names(tmp)
#'                      nms <- nms[nms == "calibrationID" | nms %in% paste0("c", 0:20) | nms %in% paste0("x", 0:20) | nms %in% paste0("n", 0:20) ]
#'                      tmp <- tmp[, nms, with = FALSE]
#'                      tmp <- suppressWarnings(melt(tmp, id.vars = "calibrationID"))
#'                      if (ncol(tmp) <= 1) {
#'                        warning("Could not find calibration parameters.  Setting keep_raw and raw to FALSE")
#'                        keep_raw <- FALSE
#'                        raw <- FALSE
#'                      } else {
#'                        setnames(tmp, c("calibrationID", "key", "value"))
#'                      }
#'                      self$coefficients <- na.omit(tmp)
#'                    }
#'                    self$channels     <- collapse::qDT(dbGetQuery(db, "SELECT * FROM channels"))
#'                    self$instruments  <- collapse::qDT(dbGetQuery(db, "SELECT * FROM instruments"))
#'                    self$errors       <- collapse::qDT(dbGetQuery(db, "SELECT * FROM errors"))
#'                    self$events       <- collapse::qDT(dbGetQuery(db, "SELECT * FROM events"))
#'
#'                    if ("continuous" %in% self$db_tables) {
#'                      self$continuous   <- collapse::qDT(dbGetQuery(db, "SELECT * FROM continuous"))
#'                    } else {
#'                      self$continuous <- collapse::qDT(dbGetQuery(db, "SELECT * FROM schedules"))
#'                    }
#'
#'                    # simple constants
#'                    self$n_chan        <- nrow(self$channels)
#'                    self$time_1        <- self$events$tstamp[[1]] / 1000 - self$rbr_time_reference()
#'                    self$measurement_interval <- self$continuous[["samplingPeriod"]] / 1000
#'                    self$model         <- self$instruments[["model"]]
#'                    self$serial        <- self$instruments[["serialID"]]
#'
#'
#'                    co_compensation    <- self$comp_coefficients()
#'                    self$channel_names(keep_raw)
#'                    self$channel_units()
#'
#'
#'                    if (raw) {
#'
#'                      # get raw values
#'                      raw_val <- unlist(lapply(self$blob$data, function(x) {
#'                        readBin(x,
#'                                n = 68000,
#'                                what = "raw")
#'                      }), use.names = FALSE, recursive = FALSE)
#'
#'
#'                      # determine non representative values
#'
#'                      # calculate mv, dbar, temperature
#'                      self$data <- collapse::qDT(
#'                        rsk:::rsk_read_bin(raw_val,
#'                                           self$is_temperature(),
#'                                           t(self$base_coefficients()),
#'                                           co_compensation,
#'                                           self$measurement_interval,
#'                                           self$pressure_index,
#'                                           self$temperature_index,
#'                                           keep_raw)
#'                      )
#'                      setnames(self$data, self$names)
#'
#'                      # subset samples
#'                      if (!is.null(by)) {
#'                        self$data <- self$data[as.numeric(datetime) %% by == 0]
#'                      }
#'
#'                      if (simplify_names) {
#'                        self <- self$rename_data()
#'                        self <- self$simplify_names()
#'                      } else if (!keep_raw) {
#'                        self <- self$rename_data()
#'                      }
#'
#'                    } else {
#'
#'                      sql_suffix <- rsk_generate_sql_times(start, end, by, times)
#'                      if(sql_suffix != "") {
#'                        sql_text <- paste0(rsk_generate_sql(db, start=NULL, end=NULL, by=NULL),
#'                                           sql_suffix)
#'                      } else {
#'                        sql_text <- rsk_generate_sql(db, start, end, by)
#'                      }
#'
#'                      self$data <- rsk_read_data_table(db, sql_text)
#'                      # setnames(self$data, self$channels$shortName, self$names)
#'                    }
#'
#'
#'                    dbDisconnect(db)
#'
#'
#'                    self$n_obs <- nrow(self$data)
#'                    self$update_time_ranges()
#'
#'                    invisible(self)
#'
#'                  },
#'                  rename_data = function() {
#'
#'
#'                    nms <- names(self$data)
#'                    nms <- intersect(nms, rbr_channels$channel_name)
#'                    nms_dt <- data.table::data.table(channel_name = nms)
#'
#'                    nms_new <- rbr_channels[nms_dt, on = "channel_name"][["field_name"]]
#'
#'                    setnames(self$data, nms, nms_new)
#'
#'                    invisible(self)
#'
#'                  },
#'                  simplify_names = function() {
#'                    # remove non-final columns
#'
#'                    nms <- names(self$data)
#'
#'                    if ("pressure_compensated" %in% nms) {
#'                      self$data[, pressure := NULL]
#'                      setnames(self$data, "pressure_compensated", "pressure", skip_absent = TRUE)
#'                    }
#'
#'                    if ("temperature" %in% nms & "temperature_onboard" %in% nms) {
#'                      self$data[, temperature_onboard := NULL]
#'                    }
#'                    setnames(self$data, "temperature_onboard", "temperature", skip_absent = TRUE)
#'
#'                    wh <- grep("raw", nms)
#'                    data.table::set(self$data, , wh, NULL)
#'
#'                    invisible(self)
#'
#'                  },
#'                  update_time_ranges = function() {
#'
#'                    self$time_start <- self$data[1]$datetime
#'                    self$time_end   <- self$data[.N]$datetime
#'                    self$has_gaps   <- ((as.numeric(self$time_end)-
#'                                           as.numeric(self$time_start)) /
#'                                          self$measurement_interval + 1) != self$n_obs
#'
#'                    invisible(self)
#'                  },
#'                  channel_names = function(keep_raw = FALSE) {
#'
#'                    base_name <- tolower(self$channels[["shortName"]])
#'
#'                    if(keep_raw) {
#'
#'                      if(self$is_temperature_compensated) {
#'                        base_name <- c("datetime",
#'                                       paste0(rep(base_name, each = 2), c("_raw", "")),
#'                                       paste0(base_name[self$pressure_index],"_compensated"))
#'                      } else {
#'                        base_name <- c("datetime", paste0(rep(base_name, each = 2), c("_raw", "")))
#'                      }
#'
#'                    } else {
#'                      if(self$is_temperature_compensated) {
#'                        base_name <- c("datetime",
#'                                       paste0(base_name),
#'                                       paste0(base_name[self$pressure_index], "_compensated"))
#'                      } else {
#'                        base_name <- c("datetime", paste0(base_name))
#'                      }
#'                    }
#'
#'                    self$names <- base_name
#'
#'                    invisible(self)
#'
#'                  },
#'                  channel_units = function() {
#'
#'                    units <- rep(NA_character_, length(self$names))
#'                    units[grepl("pressure", self$names)] <- "dbar"
#'                    units[grepl("temperature", self$names)] <- "degrees_c"
#'                    units[grepl("_raw", self$names)] <- "mv"
#'
#'                    self$units <- units
#'
#'                    invisible(self)
#'                  },
#'
#'                  # which channels are temperature
#'                  is_temperature = function() {
#'                    nms <- names(self$channels)
#'
#'                    if("longNamePlainText" %in% nms) {
#'                      field_name <- "longNamePlainText"
#'                    } else if ("longName" %in% nms) {
#'                      field_name <- "longName"
#'                    }
#'
#'                    self$channels[[field_name]] == "Temperature"
#'
#'                  },
#'
#'                  # return a matrix (or vector) of the calibration coefficients (i.e. c0,c1,...)
#'                  base_coefficients = function(id = NULL) {
#'
#'                    if(is.null(id)) {
#'
#'                      co <- matrix(
#'                        self$coefficients[grepl("c", key)][["value"]],
#'                        ncol = self$n_chan)
#'
#'                      return(co)
#'
#'                    } else {
#'
#'                      return(self$coefficients[grepl("c", key) & calibrationID == id][["value"]])
#'                    }
#'                  },
#'
#'                  # return a vector of the compensation coefficients (i.e. x0,x1,...)
#'                  comp_coefficients = function() {
#'
#'                    if("key" %in% names(self$coefficients)){
#'                      comp_cols <- self$coefficients[key == "n0"]
#'
#'                      if(nrow(comp_cols) > 0) {
#'                        if(comp_cols$value != -1) {
#'
#'                          self$pressure_index <- comp_cols$calibrationID
#'                          self$temperature_index <- comp_cols$value
#'                          self$is_temperature_compensated <- TRUE
#'
#'                          return(self$coefficients[grepl("x", key)]$value)
#'                        }
#'                      }
#'                    }
#'
#'                    self$pressure_index <- 0L
#'                    self$temperature_index <- 0L
#'
#'                    return(NA_real_)
#'
#'                  },
#'
#'                  glance = function(n_max = 1e5, type = "calculated", exclude_times = NULL) {
#'
#'                    # do you exclude times
#'                    if(is.null(exclude_times)) {
#'                      do_subset <- FALSE
#'                    } else {
#'                      do_subset <- TRUE
#'                    }
#'
#'                    # subset columns
#'                    if(type == "calculated") {
#'                      subs <- self$names[!grepl("_raw", self$names)]
#'                    } else if(type == "raw") {
#'                      subs <- self$names[grepl("_raw", self$names)]
#'                      subs <- c("datetime", subs)
#'                    } else {
#'                      subs <- self$names
#'                    }
#'                    data_sub <- self$data[, subs, with = FALSE]
#'
#'
#'                    if(self$n_obs < n_max){
#'
#'                      data_sub <- melt(data_sub, id.vars = "datetime")
#'
#'                      # exclude certain times
#'                      if(do_subset) {
#'                        # print("here")
#'                        data_sub <- data_sub[inrange(datetime, exclude_times$start, exclude_times$end), value := NA_real_]
#'                      }
#'                      data_sub <- data_sub[variable != "datetime"]
#'                      print(plot_ly(data_sub,
#'                                    x = ~datetime,
#'                                    y = ~value,
#'                                    color = ~variable,
#'                                    type = "scatter",
#'                                    mode = "lines"))
#'                    } else {
#'
#'                      data_sub <- self$data[as.integer(seq.int(1, .N, length.out = n_max)), subs, with = FALSE]
#'                      data_sub <- melt(data_sub, id.vars = "datetime")
#'
#'                      # exclude certain times
#'                      if(do_subset) {
#'                        # print("here")
#'                        data_sub <- data_sub[inrange(datetime, exclude_times$start, exclude_times$end), value := NA_real_]
#'                      }
#'                      data_sub <- data_sub[variable != "datetime"]
#'                      print(plot_ly(data_sub,
#'                                    x = ~datetime,
#'                                    y = ~value,
#'                                    color = ~variable,
#'                                    type = "scatter",
#'                                    mode = "lines"))
#'
#'                    }
#'
#'                    invisible(self)
#'                  },
#'
#'                  to_rdata = function(output = NULL, ...){
#'                    if(is.null(output)) {
#'                      output <- gsub(".rsk", ".rda", self$file_name)
#'                    }
#'                    saveRDS(self, output, ...)
#'                  },
#'
#'                  to_csv = function(output = NULL, ...){
#'
#'                    if(is.null(output)) {
#'                      output <- gsub(".rsk", ".gz", self$file_name)
#'                    }
#'                    fwrite(self$data, output, ...)
#'                  },
#'
#'                  # to_fst = function(output = NULL, ...){
#'                  #   if(is.null(output)) {
#'                  #     output <- gsub(".rsk", ".fst", self$file_name)
#'                  #   }
#'                  #   write_fst(self$data, output, ...)
#'                  # },
#'
#'                  to_duckdb = function(output = NULL){
#'                    if(is.null(output)) {
#'                      output <- gsub(".rsk", ".duckdb", x)
#'                    }
#'
#'                    duck_db <- DBI::dbConnect(duckdb::duckdb(), dbdir = output)
#'                    DBI::dbWriteTable(duck_db, "data", self$data)
#'                    DBI::dbWriteTable(duck_db, "coefficients", self$coefficients)
#'                    DBI::dbWriteTable(duck_db, "channels", self$channels)
#'                    DBI::dbWriteTable(duck_db, "instruments", self$instruments)
#'                    DBI::dbWriteTable(duck_db, "errors", self$errors)
#'                    DBI::dbWriteTable(duck_db, "events", self$events)
#'                    DBI::dbWriteTable(duck_db, "continuous", self$continuous)
#'
#'                    DBI::dbDisconnect(duck_db, shutdown=TRUE)
#'                  },
#'
#'                  non_representative = function(raw_head) {
#'
#'                    h_len  <- max(which(raw_head == self$time_1))
#'                    to_rem <- list()
#'
#'                    to_rem[[1]] <- 1:h_len
#'
#'                    for(i in 2:nrow(self$events)) {
#'                      to_rem[[i]] <- (self[["events"]][["sampleIndex"]][i]-3+i) * self$n_chan + h_len + 1:2
#'                    }
#'
#'                    sort(unlist(to_rem, use.names = FALSE))-1
#'
#'
#'                  },
#'
#'
#'                  rbr_time_reference = function() {
#'                    as.numeric(as.POSIXct("2000-01-01 00:00:00", tz = "UTC"))
#'                  })
#'
#'
#' )





#' read_levelogger
#'
#' @param x
#' @param transducer_depth
#' @param well_elevation
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
read_levelogger <- function(x,
                            transducer_depth = NULL,
                            well_elevation   = NULL,
                            ...) {

  r <- XML::xmlRoot(XML::xmlParse(x, encoding = "ISO-8859-1"))

  download_time <- .parse_datetime_levelogger(
    XML::getChildrenStrings(r[['File_info']][['Date']]),
    XML::getChildrenStrings(r[['File_info']][['Time']]))

  version_s <- XML::getChildrenStrings(r[['File_info']][['Created_by']])
  version_f <- XML::getChildrenStrings(r[['Instrument_info']][['Firmware']])

  instrument <- XML::getChildrenStrings(r[['Instrument_info']])
  header     <- XML::getChildrenStrings(r[['Instrument_info_data_header']])

  n_channel  <- as.numeric(instrument[['Channel']])


  channels <- list()
  for(i in 1:n_channel) {
    channels[[i]] <- data.table(t(
      tolower(XML::getChildrenStrings(r[[paste0('Ch', i, '_data_header')]]))))
  }
  channels <- rbindlist(channels)

  len_dat <- 3 + n_channel
  nms  <- as.character(names(r[['Data']][['Log']]))
  data <- vapply(XML::xmlChildren(r[['Data']]),
                 function(x) XML::getChildrenStrings(x, asVector = TRUE, len = len_dat), FUN.VALUE = rep(NA_character_, len_dat))
  data <- as.data.table(matrix(unlist(data),
                               ncol = len_dat,
                               byrow = TRUE,
                               dimnames = list(NULL, nms)))

  cols <- names(data)[!names(data) %chin% c('Date', 'Time')]
  data[, (cols):= lapply(.SD, as.numeric), .SDcols = cols]

  data[, datetime := .parse_datetime_levelogger(Date, Time) + as.numeric(ms)/1000]

  start <- min(data$datetime)
  end   <- max(data$datetime)

  data[, Date := NULL]
  data[, Time := NULL]
  data[, ms := NULL]
  setnames(data, c(channels$Identification, 'datetime'))

  data <- melt(data, id.vars = 'datetime')
  data <- lapply(channels$Identification, function(x) data[variable == x, ][, variable := NULL])


  dat <- data.table(
    file    = x,
    channel = channels$Identification,
    data    = data,
    id      = header['Event_ch'],
    calibration = list(data.table(coef = character(), value = numeric())),
    parameter =  channels$Identification,
    units = channels$Unit,
    version = instrument['Firmware'],
    serial = instrument['Serial_number'],
    model = paste(instrument['Instrument_type'], instrument['Model_number']),
    dt = as.numeric(header['Sample_rate']) * 10
    # n = as.numeric(header['Num_log'])
  )

  if(!is.null(well_elevation)) {
    dat[, well_elevation := well_elevation]
  }

  if(!is.null(transducer_depth)) {
    dat[, transducer_depth := transducer_depth]
  }

  # setcolorder(h, "file", "channel", "data", "id", "calibration",
  #                "parameter", "units", "version", "serial", "model", "dt")

  dat
}


.parse_datetime_levelogger <- function(x, y) {

  as.POSIXct(paste0(x, y), format = '%Y/%m/%d %H:%M:%S', tz = 'UTC')

}
