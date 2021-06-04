#
# # https://rbr-global.com/wp-content/uploads/2017/04/0001963revB-Logger2-Command-Reference.pdf
#
# # metadata section
#
# # deployment section
# firmware_version     <- get_bin_int_4_byte(con, 1L)
# serial_number        <- get_bin_int_4_byte(con, 1L)
# logger_time          <- get_bin_int_4_byte_unsigned(con)
# schedule_start_time  <- get_bin_int_4_byte_unsigned(con)
# schedule_end_time    <- get_bin_int_4_byte_unsigned(con)
# measurement_interval <- get_bin_int_4_byte(con, 1L)
# output_format        <- get_bin_int_4_byte(con, 1L)
# logger_status        <- get_bin_int_4_byte(con, 1L)
# baud_rate            <- get_bin_int_4_byte(con, 1L)
# feature_flag         <- get_bin_int_4_byte(con, 1L)
# average_interval     <- get_bin_int_4_byte(con, 1L)
# average_length       <- get_bin_int_4_byte(con, 1L)
# burst_interval       <- get_bin_int_4_byte(con, 1L)
# burst_length         <- get_bin_int_4_byte(con, 1L)
# altitude             <- get_bin_double_4_byte(con, 1L)
# thresholding_chan    <- get_bin_int_4_byte(con, 1L)
# thresholding_cond    <- get_bin_int_4_byte(con, 1L)
# thresholding_val     <- get_bin_double_4_byte(con, 1L)
# thresholding_int     <- get_bin_int_4_byte(con, 1L)
# fetch_power_off_delay <- get_bin_int_4_byte(con, 1L)
# temperature_default   <- get_bin_double_4_byte(con, 1L)
# conductivity_default  <- get_bin_double_4_byte(con, 1L)
# pressure_default      <- get_bin_double_4_byte(con, 1L)
# atm_pressure_default  <- get_bin_double_4_byte(con, 1L)
# density_default       <- get_bin_double_4_byte(con, 1L)
# regime       <- get_bin_double_4_byte(con, 1L)
#
