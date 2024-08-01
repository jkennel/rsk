#' @export
Rsk <- R6Class("Rsk",
               public = list(

  # database name
  file_name     = NA_character_,


  # tables from rsk
  blob         = NULL,
  channels     = NULL,
  errors       = NULL,
  instruments  = NULL,
  events       = NULL,
  coefficients = NULL,
  continuous   = NULL,


  # table of measurements
  data          = NULL,


  time_1        = NULL,
  measurement_interval = NULL,

  names         = NULL,
  units         = NULL,

  pressure_index    = 0,
  temperature_index = 0,
  is_temperature_compensated = FALSE,

  n_obs         = NULL,
  n_chan        = NULL,

  serial        = NULL,
  model         = NULL,

  latitude      = NULL,
  longitude     = NULL,
  elevation     = NULL,

  time_start    = NULL,
  time_end      = NULL,
  has_gaps      = NULL,

  db_tables     = NULL,



#' @description
#' Class to hold data from an .rsk or .duckdb file
#'
#' @param file_name
#' @param subset
#' @param latitude
#' @param longitude
#' @param elevation
#' @param start
#' @param end
#' @param by
#' @param times
#' @param raw
#' @param simplify_names
#'
#' @return
#' @export
#'
#' @examples
initialize = function(file_name,
                      subset = FALSE,
                      latitude = NULL,
                      longitude = NULL,
                      elevation = NULL,
                      start = NULL,
                      end = NULL,
                      by = NULL,
                      times = NULL,
                      raw = FALSE,
                      keep_raw = FALSE,
                      simplify_names = FALSE
) {

  self$file_name <- file_name

  # cannot keep raw data if you aren"t going from the binary table
  if (!raw) {
    keep_raw <- FALSE
  }


  # get data base tables
  if (tools::file_ext(file_name) == "rsk") {
    db <- dbConnect(SQLite(), file_name)
    self$db_tables    <- RSQLite::dbListTables(db)
    if ("downloads" %in% self$db_tables) {
      self$blob <- setDT(dbGetQuery(db, "SELECT * FROM downloads"))
    } else {
      if ("data" %in% self$db_tables) {

        keep_raw <- FALSE
        raw <- FALSE

      } else {
        dbDisconnect(db)

        return(self)
      }
    }
  } else if (tools::file_ext(file_name) == "duckdb") {
    if (raw) {
      stop("raw cannot be TRUE if file is type .duckdb")
    }
    db <- dbConnect(duckdb(), file_name)
  } else {
    dbDisconnect(db)
    stop("File type not currently supported")
  }
  self$db_tables    <- RSQLite::dbListTables(db)

  if ("coefficients" %in% self$db_tables) {
    self$coefficients <- setDT(dbGetQuery(db, "SELECT calibrationID, key, cast(value as REAL) as value FROM coefficients"))
  } else if ("calibrations" %in% self$db_tables) {
    tmp <- setDT(dbGetQuery(db, "SELECT * FROM calibrations"))
    nms <- names(tmp)
    nms <- nms[nms == "calibrationID" | nms %in% paste0("c", 0:20) | nms %in% paste0("x", 0:20) | nms %in% paste0("n", 0:20) ]
    tmp <- tmp[, nms, with = FALSE]
    tmp <- suppressWarnings(melt(tmp, id.vars = "calibrationID"))
    if (ncol(tmp) <= 1) {
      warning("Could not find calibration parameters.  Setting keep_raw and raw to FALSE")
      keep_raw <- FALSE
      raw <- FALSE
    } else {
      setnames(tmp, c("calibrationID", "key", "value"))
    }
    self$coefficients <- na.omit(tmp)
  }
  self$channels     <- setDT(dbGetQuery(db, "SELECT * FROM channels"))
  self$instruments  <- setDT(dbGetQuery(db, "SELECT * FROM instruments"))
  self$errors       <- setDT(dbGetQuery(db, "SELECT * FROM errors"))
  self$events       <- setDT(dbGetQuery(db, "SELECT * FROM events"))

  if ("continuous" %in% self$db_tables) {
    self$continuous   <- setDT(dbGetQuery(db, "SELECT * FROM continuous"))
  } else {
    self$continuous <- setDT(dbGetQuery(db, "SELECT * FROM schedules"))
  }

  # simple constants
  self$n_chan        <- nrow(self$channels)
  self$time_1        <- self$events$tstamp[[1]] / 1000 - self$rbr_time_reference()
  self$measurement_interval <- self$continuous[["samplingPeriod"]] / 1000
  self$model         <- self$instruments[["model"]]
  self$serial        <- self$instruments[["serialID"]]


  co_compensation    <- self$comp_coefficients()
  self$channel_names(keep_raw)
  self$channel_units()


  if (raw) {

    # get raw values
    raw_val <- unlist(lapply(self$blob$data, function(x) {
      readBin(x,
              n = 68000,
              what = "raw")
    }))


    # determine non representative values

    # calculate mv, dbar, temperature
    self$data <- setDT(
      rsk:::rsk_read_bin(raw_val,
                         self$is_temperature(),
                         t(self$base_coefficients()),
                         co_compensation,
                         self$measurement_interval,
                         self$pressure_index,
                         self$temperature_index,
                         keep_raw)
    )
    setnames(self$data, self$names)

    # subset samples
    if (!is.null(by)) {
      self$data <- self$data[as.numeric(datetime) %% by == 0]
    }

    if (simplify_names) {
      self$data <- self$rename_data()
      self$data <- self$simplify_names()
    }

  } else {

    sql_suffix <- rsk_generate_sql_times(start, end, by, times)
    if(sql_suffix != "") {
      sql_text <- paste0(rsk_generate_sql(db, start=NULL, end=NULL, by=NULL),
                         sql_suffix)
    } else {
      sql_text <- rsk_generate_sql(db, start, end, by)
    }

    self$data <- rsk_read_data_table(db, sql_text)
    # setnames(self$data, self$channels$shortName, self$names)
  }


  dbDisconnect(db)


  self$n_obs <- nrow(self$data)
  self$update_time_ranges()

  invisible(self)

},
rename_data = function() {

  data(rbr_channels)

  nms <- names(self$data)
  nms <- intersect(nms, rbr_channels$channel_name)
  nms_dt <- data.table::data.table(channel_name = nms)

  nms_new <- rbr_channels[nms_dt, on = "channel_name"][["field_name"]]

  setnames(self$data, nms, nms_new)

  invisible(self$data)

},
simplify_names = function() {
  # remove non-final columns

  nms <- names(self$data)

  if ("pressure_compensated" %in% nms) {
    self$data[, pressure := NULL]
    setnames(self$data, "pressure_compensated", "pressure", skip_absent = TRUE)
  }

  if ("temperature" %in% nms & "temperature_onboard" %in% nms) {
    self$data[, temperature_onboard := NULL]
  }
  setnames(self$data, "temperature_onboard", "temperature", skip_absent = TRUE)

  wh <- grep("raw", nms)
  data.table::set(self$data, , wh, NULL)

  self$data
  # invisible(self$data[, -wh])

},
update_time_ranges = function() {

  self$time_start <- self$data[1]$datetime
  self$time_end   <- self$data[.N]$datetime
  self$has_gaps   <- ((as.numeric(self$time_end)-
                         as.numeric(self$time_start)) /
                        self$measurement_interval + 1) != self$n_obs

  invisible(self)
},
channel_names = function(keep_raw = FALSE) {

  base_name <- tolower(self$channels[["shortName"]])

  if(keep_raw) {

    if(self$is_temperature_compensated) {
      base_name <- c("datetime",
                     paste0(rep(base_name, each = 2), c("_raw", "")),
                     paste0(base_name[self$pressure_index],"_compensated"))
    } else {
      base_name <- c("datetime", paste0(rep(base_name, each = 2), c("_raw", "")))
    }

  } else {
    if(self$is_temperature_compensated) {
      base_name <- c("datetime",
                     paste0(base_name),
                     paste0(base_name[self$pressure_index], "_compensated"))
    } else {
      base_name <- c("datetime", paste0(base_name))
    }
  }

  self$names <- base_name

  invisible(self)

},
channel_units = function() {

  units <- rep(NA_character_, length(self$names))
  units[grepl("pressure", self$names)] <- "dbar"
  units[grepl("temperature", self$names)] <- "degrees_c"
  units[grepl("_raw", self$names)] <- "mv"

  self$units <- units

  invisible(self)
},

# which channels are temperature
is_temperature = function() {
  nms <- names(self$channels)

  if("longNamePlainText" %in% nms) {
    field_name <- "longNamePlainText"
  } else if ("longName" %in% nms) {
    field_name <- "longName"
  }

  self$channels[[field_name]] == "Temperature"

},

# return a matrix (or vector) of the calibration coefficients (i.e. c0,c1,...)
base_coefficients = function(id = NULL) {

  if(is.null(id)) {

    co <- matrix(
      self$coefficients[grepl("c", key)][["value"]],
      ncol = self$n_chan)

    return(co)

  } else {

    return(self$coefficients[grepl("c", key) & calibrationID == id][["value"]])
  }
},

# return a vector of the compensation coefficients (i.e. x0,x1,...)
comp_coefficients = function() {

  if("key" %in% names(self$coefficients)){
    comp_cols <- self$coefficients[key == "n0"]

    if(nrow(comp_cols) > 0) {
      if(comp_cols$value != -1) {

        self$pressure_index <- comp_cols$calibrationID
        self$temperature_index <- comp_cols$value
        self$is_temperature_compensated <- TRUE

        return(self$coefficients[grepl("x", key)]$value)
      }
    }
  }

  self$pressure_index <- 0L
  self$temperature_index <- 0L

  return(NA_real_)

},

glance = function(n_max = 1e5, type = "calculated", exclude_times = NULL) {

  # do you exclude times
  if(is.null(exclude_times)) {
    do_subset <- FALSE
  } else {
    do_subset <- TRUE
  }

  # subset columns
  if(type == "calculated") {
    subs <- self$names[!grepl("_raw", self$names)]
  } else if(type == "raw") {
    subs <- self$names[grepl("_raw", self$names)]
    subs <- c("datetime", subs)
  } else {
    subs <- self$names
  }
  data_sub <- self$data[, subs, with = FALSE]


  if(self$n_obs < n_max){

    data_sub <- melt(data_sub, id.vars = "datetime")

    # exclude certain times
    if(do_subset) {
      # print("here")
      data_sub <- data_sub[inrange(datetime, exclude_times$start, exclude_times$end), value := NA_real_]
    }
    data_sub <- data_sub[variable != "datetime"]
    print(plot_ly(data_sub,
                  x = ~datetime,
                  y = ~value,
                  color = ~variable,
                  type = "scatter",
                  mode = "lines"))
  } else {

    data_sub <- self$data[as.integer(seq.int(1, .N, length.out = n_max)), subs, with = FALSE]
    data_sub <- melt(data_sub, id.vars = "datetime")

    # exclude certain times
    if(do_subset) {
      # print("here")
      data_sub <- data_sub[inrange(datetime, exclude_times$start, exclude_times$end), value := NA_real_]
    }
    data_sub <- data_sub[variable != "datetime"]
    print(plot_ly(data_sub,
                  x = ~datetime,
                  y = ~value,
                  color = ~variable,
                  type = "scatter",
                  mode = "lines"))

  }

  invisible(self)
},

to_rdata = function(output = NULL, ...){
  if(is.null(output)) {
    output <- gsub(".rsk", ".rda", self$file_name)
  }
  saveRDS(self, output, ...)
},

to_csv = function(output = NULL, ...){

  if(is.null(output)) {
    output <- gsub(".rsk", ".gz", self$file_name)
  }
  fwrite(self$data, output, ...)
},

# to_fst = function(output = NULL, ...){
#   if(is.null(output)) {
#     output <- gsub(".rsk", ".fst", self$file_name)
#   }
#   write_fst(self$data, output, ...)
# },

to_duckdb = function(output = NULL){
  if(is.null(output)) {
    output <- gsub(".rsk", ".duckdb", x)
  }

  duck_db <- DBI::dbConnect(duckdb::duckdb(), dbdir = output)
  DBI::dbWriteTable(duck_db, "data", self$data)
  DBI::dbWriteTable(duck_db, "coefficients", self$coefficients)
  DBI::dbWriteTable(duck_db, "channels", self$channels)
  DBI::dbWriteTable(duck_db, "instruments", self$instruments)
  DBI::dbWriteTable(duck_db, "errors", self$errors)
  DBI::dbWriteTable(duck_db, "events", self$events)
  DBI::dbWriteTable(duck_db, "continuous", self$continuous)

  DBI::dbDisconnect(duck_db, shutdown=TRUE)
},

non_representative = function(raw_head) {

  h_len  <- max(which(raw_head == self$time_1))
  to_rem <- list()

  to_rem[[1]] <- 1:h_len

  for(i in 2:nrow(self$events)) {
    to_rem[[i]] <- (self[["events"]][["sampleIndex"]][i]-3+i) * self$n_chan + h_len + 1:2
  }

  sort(unlist(to_rem))-1


},


rbr_time_reference = function() {
  as.numeric(as.POSIXct("2000-01-01 00:00:00", tz = "UTC"))
})


)


# #need to check 63
#
# library(rsk)
# library(R6)
# library(RSQLite)
# library(data.table)
# library(tools)
#
# fn <- list.files("/media/jonathankennel/Seagate Expansion Drive/guelph_south",
#                  full.names = TRUE)
# fn <- fn[tools::file_ext(fn) == "rsk"]
# i  <- 1
#
#
# file_name <- "/media/jonathankennel/Seagate Expansion Drive/rbr_g360/rd45a 081871_20201026_1223.rsk"
# file_name <- "/media/jonathankennel/Seagate Expansion Drive/rbr_g360/RD130_077623_20190307_1312.rsk"
#
# zz <- Rsk$new(fn[1], raw = TRUE, keep_raw = TRUE)
#
#
# for(i in seq_along(fn)) {
#   zz <- Rsk$new(fn[i], raw = TRUE, keep_raw = TRUE)
#   print(zz$names)
# }
#
#
#
#
# system.time(z <- Rsk$new(fn[i], raw = TRUE, keep_raw = TRUE))
#
# system.time(
#
#   (bb <- rsk_parse_binary(fn[i], keep_raw = TRUE))
# )


# tmp <- fread("/media/jonathankennel/Seagate Expansion Drive/guelph_south/meta/transducer_subsets.csv")
# tmp <- tmp[file_name == basename(fn[i])]
# tmp[, start := as.POSIXct(start, tz = "UTC")]
# tmp[, end := as.POSIXct(end, tz = "UTC")]
#
#
# system.time(z <- Rsk$new(fn[i]))
# print(fn[i])
# z$time_start
# z$time_end
# z$glance(n_max = 1e5, type = "calculated", exclude_times = tmp)



# library(transducer)
# library(rsk)
# library(RSQLite)
# library(data.table)


# file_name <- "/media/jonathankennel/Seagate Expansion Drive/rbr_g360/077753_20211115_1043.rsk"
# file_name <- paste0("/media/jonathankennel/Seagate Expansion Drive/rbr_g360/", well_fn[8])
# Rsk$new(well_fn[9], raw = TRUE, keep_raw = TRUE)


#
#
# db <- dbConnect(SQLite(), file_name)
# # dbListTables(db)
#
# data        <- setDT(dbGetQuery(db, "SELECT * FROM data"))
#
# data[, tstamp := as.double(tstamp)]
# fwrite(data, "data.csv")
#
# tmp <- fread("data.csv")
#
# dbDisconnect(db)
#
#
# z <- Rsk$new(file_name)
# down <- z$data[, list(tstamp = as.integer(datetime)*1000, temperature, pressure_dbar_comp, temperature_onboard)]
# down[, tstamp := as.double(tstamp)]
# fwrite(down, "downloads.csv")
# tmp2 <- fread("downloads.csv")
# tmp2[, id := 1:.N]
# differences <- tmp2[!J(tmp), on = "tstamp"]
# differences[, datetime := as.POSIXct(tstamp/1000, origin = "1970-01-01")]
# fwrite(differences, "differences.csv")
# plot(temperature~tstamp, differences, type = "l")
# plot(pressure_dbar_comp~tstamp, differences, type = "l")
#
#
# plot(pressure_dbar_comp~tstamp, down[seq(1, .N, 100)], type="l")

 # fwrite(down, "downloads.csv")
# library(rsk)
# library(R6)
# library(data.table)
# library(RSQLite)
#
# x <- "/home/jonathankennel/Storage/data/rbr/rd45a 081871_20191118_1213.rsk"
# x <- "/home/jonathankennel/Storage/data/rbr/rd130 077623_20191119_1500.rsk"
# x <- "/home/jonathankennel/Storage/data/rbr/rd45a 081871_20191118_1213.rsk"
# system.time(z1 <- Rsk$new(x))
# system.time(z2 <- Rsk$new(x, raw = TRUE))
# system.time(z <- Rsk$new(x)$data)
# system.time(zz <- transducer::read_rbr(x)$data[[1]])
# profvis::profvis(z <- Rsk$new(x))
# zz <- z$data
#
# plot(V2~datetime, zz[40:2000], type="l")
# points(V6~datetime, zz[40:2000], type="l", col = "red")
# # str(z$base_coefficients(id = 1))
# # str(z$comp_coefficients(id = 2))

# library(data.table)
# library(duckdb)
# library(DBI)
# x <- "/home/jonathankennel/Storage/data/rbr/rd45a 081871_20191118_1213.rsk"
# system.time(z1 <- Rsk$new(x))
# system.time(z1$to_duck())
# x <- "/home/jonathankennel/Storage/data/rbr/rd45a 081871_20191118_1213.duckdb"
# system.time({
# db <- duckdb::dbConnect(duckdb::duckdb(), dbdir = x, read_only = TRUE)
# a <- setDT(dbGetQuery(db, "SELECT * FROM data"))
# dbDisconnect(db, shutdown=TRUE)
# })


