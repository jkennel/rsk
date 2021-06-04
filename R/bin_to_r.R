# library(data.table)
# library(here)
# # system.time(
# # hex_dat <- readBin(x, n = 2e8, what = 'raw')
# # )
#
# dat_dir <-  '/media/jonathankennel/Seagate Expansion Drive/rbr_g360/'
# fn <- list.files(dat_dir, full.names = TRUE, patter = '.rsk')
# fn <- fn[grepl('GSMW2', fn)]
# x  <- fn[11]
#
# # system.time(
# #
# #   q <- rbr_raw_times(raw_tstamp = times,
# #                      raw_index = tmp,
# #                      n_raw = length(dat),
# #                      n_to_rem = length(idx),
# #                      ti = 1.0)
# # )
#
# index_to_mask <- function(x, n_chan) {
#
#   # incomplete rows have non-zero remainder
#   inc_to_mask <- (((diff(x) - 8) %% (4 * n_chan))) / 4
#
#   idx <- list()
#   for (i in seq_along(x)) {
#     if(i > 1)
#       if(inc_to_mask[i-1] > 0) {
#         n <- inc_to_mask[i-1]
#         idx[[i]] <- rep(x[i] %/% 4, each = 2 + n) + (1-n):2
#       } else {
#         idx[[i]] <- rep(x[i] %/% 4, each = 2) + 1:2
#       } else {
#         idx[[i]] <- rep(x[i] %/% 4, each = 2) + 1:2
#       }
#   }
#
#   return(as.integer(unlist(idx)-1))
#
# }
#
#
# get_bin_int_4_byte_unsigned <- function(con) {
#
#   rsk:::raw_to_4byte_unsigned(
#     readBin(con, n = 4, what = 'raw'))
#
# }
#
# get_bin_int_1_byte <- function(con, n) {
#   readBin(con,
#           n = n,
#           what = 'integer',
#           size = 1L,
#           signed = TRUE,
#           endian = 'little')
# }
# get_bin_int_1_byte_unsigned <- function(con, n) {
#   readBin(con,
#           n = n,
#           what = 'integer',
#           signed = FALSE,
#           size = 1L,
#           endian = 'little')
# }
#
# get_bin_int_2_byte <- function(con, n) {
#   readBin(con,
#           n = n,
#           what = 'integer',
#           size = 2L,
#           signed = TRUE,
#           endian = 'little')
# }
#
# get_bin_int_4_byte <- function(con, n) {
#   readBin(con,
#           n = n,
#           what = 'integer',
#           size = 4L,
#           signed = TRUE,
#           endian = 'little')
# }
#
#
# get_bin_double_4_byte <- function(con, n) {
#   readBin(con,
#           n = n,
#           what = 'numeric',
#           size = 4L,
#           signed = TRUE,
#           endian = 'little')
#
# }
#
#
# get_bin_character <- function(con, n, null_ended = TRUE) {
#   raw_hex <- readBin(con,
#                      n = n,
#                      what = 'raw')
#
#   if(null_ended) {
#     wh <- which(raw_hex == '00')
#     return(rawToChar(raw_hex[1:(wh[1]-1)]))
#   }
#
#   else {
#     return(rawToChar(raw_hex))
#   }
#
# }
#
#
# parse_channel <- function(con, i) {
#
#   name              <- get_bin_character(con, 6L, null_ended = FALSE)
#   label             <- get_bin_character(con, 32L)
#   fwtype            <- get_bin_character(con, 32L)
#   firmware          <- get_bin_int_4_byte(con, 1L)
#   channel_ext       <- get_bin_int_2_byte(con, 1L)
#   calibration_date  <- get_bin_int_4_byte_unsigned(con)
#   n_coefficients    <- get_bin_int_1_byte(con, 1L)
#
#   if(grepl('DUET-temp', fwtype)) {
#     key <- paste0('c', 0:3)
#     coefficients      <- get_bin_double_4_byte(con, n_coefficients)
#   } else if(grepl('DUET-prs', fwtype)) {
#     key <- c(paste0('c', 0:3),
#              paste0('x', 0:5),
#              paste0('n', '0'))
#     coefficients  <- c(get_bin_double_4_byte(con, n_coefficients-2),
#                        get_bin_int_4_byte(con, 2)[1])
#   } else if(grepl('SOLO-comp', fwtype)) {
#     key <- c(paste0('c', 0:3))
#     coefficients      <- get_bin_double_4_byte(con, n_coefficients)
#   }
#
#   channel_size      <- get_bin_int_2_byte(con, 1L)
#
#
#   if(channel_size > 0) {
#
#     channel_type        <- get_bin_int_1_byte(con, 1L)
#     channel_type_length <- get_bin_int_2_byte(con, 2L)
#
#     if(channel_type == 2) {
#       # key value
#       tmp <- readBin(con, n = channel_type_length[1]-5, what = 'raw')
#       tmp <- tmp[tmp != '7f']
#       wh  <- c(0, which(tmp == '00'))
#       key_val <- c()
#       for(i in 1:(length(wh) - 1)) {
#         key_val[i] <- readBin(tmp[ (wh[i]+1):wh[i+1] ], what = 'character')
#       }
#
#     } else if (channel_type == 3){
#     }
#   }
#
#   channel <- data.table(
#     channelID = i,
#     shortName = name,
#     label,
#     fwtype,
#     firmware,
#     calibration_date
#   )
#
#   channel[grepl('temp', shortName), longName := 'Temperature']
#   channel[grepl('temp', shortName), longNamePlainText := 'Temperature']
#   channel[grepl('temp', shortName), units := '°C']
#   channel[grepl('temp', shortName), unitsPlainText := 'Degrees_C']
#
#   channel[grepl('pres', shortName), longName := 'Pressure']
#   channel[grepl('pres', shortName), longNamePlainText := 'Pressure']
#   channel[grepl('pres', shortName), units := 'dbar']
#   channel[grepl('pres', shortName), unitsPlainText := 'dbar']
#
#   calibrations <-   data.table(
#     calibrationID = i,
#     k = key,
#     value = coefficients)
#   list(
#     channel,
#     calibrations
#   )
# }
#
#
# rsk_parse_binary <- function(x, keep_raw = FALSE) {
#
#   if(tools::file_ext(x) == 'rsk') {
#
#     db  <- RSQLite::dbConnect(RSQLite::SQLite(), x)
#     blob <- RSQLite::dbGetQuery(db, "SELECT data FROM downloads")[[1]]
#     RSQLite::dbDisconnect(db)
#
#     hex_dat <- unlist(lapply(blob, function(x) {
#       readBin(x,
#               n = 68000,
#               what = 'raw')
#     }))
#
#     con <- rawConnection(hex_dat, "rb") # start with empty raw vector
#
#   } else if(tools::file_ext(x) == 'bin'){
#     con = file(x, "rb")
#   }
#
#
#   # section 01 - metadata section
#
#   id_01            <- get_bin_int_1_byte(con, 1L)
#   length_id        <- get_bin_int_2_byte(con, 1L)
#   version          <- get_bin_int_4_byte(con, 1L)
#   header_length    <- get_bin_int_2_byte(con, 1L)
#
#
#   # section 02 - logger section
#     id_02            <- get_bin_int_1_byte(con, 1L)
#     length_id        <- get_bin_int_2_byte(con, 1L)
#
#     # if length is 503 this is the logger 2 format
#     if(length_id == 503) {
#     } else { # use logger 3 format
#       firmware_type    <- get_bin_int_4_byte(con, 1L)
#       firmware_version <- get_bin_int_4_byte(con, 1L)
#       serial_number    <- get_bin_int_4_byte(con, 1L)
#       logger_model     <- get_bin_character(con, 16L)
#       logger_pn_length <- get_bin_int_2_byte(con, 1L)
#       logger_pn        <- get_bin_character(con, logger_pn_length)
#       power_pn_length  <- get_bin_int_2_byte(con, 1L)
#
#       if(power_pn_length > 0) {
#         power_pn <- readBin(con,
#                             n = power_pn_length,
#                             what = 'raw',
#                             size = 1L)
#       }
#       tmp              <- readBin(con, n = 4, what = 'raw')
#
#     }
#
#
#   }
#
#   #-------------------------------------------------------------------------------
#   instruments <- data.table(
#     instrumentID    = 1,
#     serialID        = serial_number,
#     model           = logger_model,
#     firmwareVersion = firmware_version/1000,
#     firmwareType    = firmware_type,
#     partNumber      = logger_pn
#   )
#   #-------------------------------------------------------------------------------
#
#   # section 03 - deployment section
#   id_03                <- get_bin_int_1_byte(con, 1L)
#   length_id            <- get_bin_int_2_byte(con, 1L)
#   # 0 for rawbin, 1 for calbin
#   dataset_format       <- get_bin_int_4_byte(con, 1L)
#   logger_time          <- get_bin_int_4_byte_unsigned(con)
#   schedule_start_time  <- get_bin_int_4_byte_unsigned(con)
#   schedule_end_time    <- get_bin_int_4_byte_unsigned(con)
#   measurement_interval <- get_bin_int_4_byte(con, 1L)
#   # 1= pending, 2 = logging, or 4 = gated
#   logger_status        <- get_bin_int_4_byte(con, 1L)
#
#   #----------------------------------------------------
#   # currently skip this portion - can come back to this
#   flag <- readBin(con,
#                   n = 4L,
#                   what = 'raw')
#   tmp <- readBin(con,
#                  n = length_id-31,
#                  what = 'raw')
#   #----------------------------------------------------
#
#   epochs <- data.table(
#     deploymentID = 1,
#     startTime = schedule_start_time,
#     endTime = schedule_end_time
#   )
#
#
#
#   # section 04 - Non deployment related settings
#   id_04       <- get_bin_int_1_byte(con, 1L)
#   length_id   <- get_bin_int_2_byte(con, 1L)
#   nd_settings <- get_bin_int_4_byte(con, (length_id - 3L) / 4L)
#
#
#   # section 05 - Channel section
#
#   id_05          <- get_bin_int_1_byte(con, 1L)
#   length_id      <- get_bin_int_2_byte(con, 1L)
#   n_chan         <- get_bin_int_1_byte(con, 1L)
#   channel_offset <- get_bin_int_2_byte(con, n_chan)
#   channel_offset <- channel_offset - channel_offset[1] + 1
#
#   coefficients <- list()
#
#   for(i in 1:n_chan) {
#     coefficients[[i]] <- parse_channel(con, i)
#   }
#
#   channels     <- rbindlist(lapply(coefficients, '[[', 1))
#   coefficients <- rbindlist(lapply(coefficients, '[[', 2))
#   setnames(coefficients, 'k', 'key')
#
#   # section 06 - Channel data
#
#   # why do I need to read this?
#   id <- readBin(con,
#                 n = 2L,
#                 what = 'raw')
#
#
#   if(tools::file_ext(x) == 'bin') {
#     hex_dat <- readBin(con, n = 2e8, what = 'raw')
#     header_length <- 0
#   }
#
#   close(con)
#
#   base_calib <- matrix(coefficients[grepl('c', key)]$value,
#                        nrow = n_chan, byrow = TRUE)
#   base_comp <- coefficients[grepl('x', key)]$value
#
#   pressure_index = coefficients[key == 'n0']$calibrationID
#   temperature_index = coefficients[key == 'n0']$value
#
#   is_temp = c(TRUE, FALSE, TRUE)
#
#   data <- setDT(rsk_read_bin(hex_dat,
#                              is_temp,
#                              base_calib,
#                              base_comp,
#                              measurement_interval = measurement_interval * 1e-3,
#                              pressure_index = pressure_index,
#                              temperature_index = temperature_index,
#                              keep_raw = keep_raw))
#
#   if(keep_raw) {
#
#     if(length(base_comp) > 0) {
#       names <- c('datetime',
#                  paste0(rep(channels$shortName, each = 2), c('_raw', '')),
#                  paste0(channels$shortName[pressure_index],'_compensated'))
#     } else {
#       names <- c('datetime', paste0(rep(channels$shortName, each = 2),
#                                     c('_raw', '')))
#     }
#
#   } else {
#     if(length(base_comp) > 0) {
#       names <- c('datetime',
#                  paste0(channels$shortName),
#                  paste0(channels$shortName[pressure_index], '_compensated'))
#     } else {
#       names <- c('datetime', paste0(channels$shortName))
#     }
#   }
#
#   setnames(data, names)
#
# }
#
#
# system.time(
#   (aa <- rsk_parse_binary(x, keep_raw = TRUE))
# )
#
# system.time(
#   (bb <- rsk_parse_binary(file_name, keep_raw = TRUE))
# )
#
#
#
# library(rsk)
# library(RSQLite)
# library(data.table)
# library(tools)
#
# file_name <- '/media/jonathankennel/Seagate Expansion Drive/rbr_g360/rd45a 081871_20201026_1223.rsk'
#
# system.time(
#
# zz <- Rsk$new(file_name, raw = TRUE, keep_raw = TRUE)
#
# )
#
#
# system.time(
#   (aa <- rsk_parse_binary(file_name, keep_raw = TRUE))
# )
#
# # zz$instruments
# #
# # tmp <- head(zz$blob$data[[1]], 670)
# #
# # readBin(tail(zz$blob$data[[nrow(zz$blob)]], 8),
# #         n = 2L,
# #         what = 'integer',
# #         size = 4L,
# #         signed = TRUE,
# #         endian = 'little')
# #
# # tmp <- tail(zz$blob$data[[1]], 8)
# #
# #
# #
# # which(tmp =='00')
#
#
