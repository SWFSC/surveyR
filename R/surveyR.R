#' Calculate visual transect width from camera pitch, altitude, and zoom level.
#'
#' @param pitch Camera pitch in degrees (positive or negative; zero equals horizontal).
#' @param altitude Raw altitude from the altimeter or similar device.
#' @param ROV ROV on which camera is installed.
#' @param z Analog voltage output from camera zoom indicator.
#' @param max.sr Maximum slant range or viewing distance in meters.
#' @return A data frame containing the \code{camera_altitude}, \code{slant_range}, and \code{center_width} in meters.
#' @export
calc_width <- function(pitch, altitude, ROV="HDHV", z=0.06, max.sr=6) {
  if (ROV == "HDHV") {
    # Empirically measured viewing angles for the HDHV ROV
    Alpha      <- -122.79 * z^4 + 132.93 * z^3 + 81.984 * z^2 - 142.07 * z + 46.159 # vertical range from 64 deg (full-wide) to 6.6 deg (10X zoom)
    Beta       <- -297.54 * z^4 + 425 * z^3 - 15.571 * z^2 - 194 * z + 74.416 # horiz range from 73.9 deg (full-wide) to 2.26 deg (10X zoom)
    Beta.max   <- 73.95 # max horiz viewing angle at full-wide zoom; measured empirically
    Beta.min   <- 2.26 # min horiz viewing angle at 10x zoom; measured empirically
    Beta.spec  <- 85 # max horiz viewing angle at full-wide zoom; manufacturer spec
    Alpha.max  <- 45.9 # max vertical viewing angle at full-wide zoom; measured empirically
    Alpha.min  <- 1.27 # min vertical viewing angle at 10x zoom; measured empirically
    Alpha.spec <- 64 # max vertical viewing angle at full-wide zoom; manufacturer spec
    Beta.diff  <- 100 * (1 - (tan((Beta.max / 2) * (pi / 180)) / tan((Beta.spec / 2) * (pi / 180)))) # difference in area between spec and empirical FOV
    v.offset   <- 0.4127 # altitude of the camera above the frame; added to altitude from DVL to get camera height
    lens.dist  <- 0.1232 # distance (m) of camera lens from tilt tray axis
    pitch.diff <- 415.5 # difference (deg) between lens pitch and all other pitch values
  } else {
    # Empirically measured viewing angles for the Phantom DS4 ROV
    Alpha      <- 47.3 # max vertical viewing angle at full-wide zoom; measured empirically
    Beta       <- 61.1 # max horiz viewing angle at full-wide zoom; measured empirically
    Beta.spec  <- 85 # max horiz viewing angle at full-wide zoom; manufacturer spec
    Alpha.spec <- 64 # max vertical viewing angle at full-wide zoom; manufacturer spec
    v.offset   <- 0.3429 # altitude of the camera above the frame; added to altitude from DVL to get camera height
    lens.dist  <- 0.1060 # distance (m) of camera lens from tilt tray axis
    pitch.diff <- 410.0 # difference (deg) between lens pitch and all other pitch values
  }
  # Calculate altitude of camera above the seabed (distance between the bottom of the frame and the camera lens) from smoothed altitude
  # Accounts for rotation of the camera module about the tilt tray axis
  alt <- altitude + v.offset + (lens.dist * (sin(pitch + pitch.diff) * (pi / 180))) # altitude of the camera above the frame

  # calculate the slant range (sr) from the camera to the seabed, using smoothed pitch values
  sr <- -alt / sin(pitch * (pi / 180))

  # calculate horizontal width at the center of the camera (or field of view, FOV)
  cw <- 2 * sr * tan((Beta / 2) * (pi / 180))

  # output results as a data frame
  data.frame(camera_altitude = alt, slant_range = sr, center_width = cw)
}

#' Convert latitude and longitude from Winfrog format to decimal degrees.
#'
#' @param x Latitude or longitude values in Winfrog format.
#' @return \code{y} in decimal degrees.
#' @export
winfrog2dd <- function(x) {
  if (length(grep("N", x)) > 0) {
    y <- as.numeric(substr(x, 2, 3)) + as.numeric(substr(x, 5, 11)) / 60
  } else {
    y <- -(as.numeric(substr(x, 2, 5)) + as.numeric(substr(x, 6, 12)) / 60)
  }
  # Return values in decimal degrees
  return(y)
}

#' Calculate dissolved oxygen percent saturation from temperature, salinity, and oxygen concentration.
#'
#' @param sal Salinity (PSU).
#' @param temp Temperature (degrees C).
#' @param oxy.conc Oxygen concentration (\code{micromoles}).
#' @return Percent saturation  \code{oxy_sat} (percent).
#' @export
calc_sat <- function(sal, temp, oxy.conc) {
  # calculate dissolved oxygen saturation from temperature and sal
  # From Weiss, R.F.(1970), The solubility of nitrogen, oxygen and argon in water and seawater.
  # Deep Sea Research 20:291-303

  # Equation details
  # ln[O2] = A1 + A2*(100/T) + A3*ln(T/100) + A4*(T/100) + S*[B1 + B2*(T/100) + B3*((T/100)^2)];
  # equation constants
  a1 <- -173.4292
  a2 <-  249.6339
  a3 <-  143.3483
  a4 <-  -21.8492
  b1 <-   -0.033096
  b2 <-    0.014259
  b3 <-   -0.0017000
  # temperature and salinity data from CTD
  S <- sal
  T <- temp
  K <- T + 273.16
  # calculate theoretical oxygen saturation at a given T and S
  DOSAT_ml <- exp(a1 + a2 * (100 / K) + a3 * (log(K / 100)) + a4 * (K / 100) + S * (b1 + b2 * (K / 100) + b3 * ((K / 100)^2)))
  DOSAT_umoles <- DOSAT_ml / 0.0224
  oxy.sat <- (oxy.conc / DOSAT_umoles) * 100
  return(oxy.sat)
}

#' Calculate depth from pressure in SEAWATER (factors-in gravity and latitude).
#'
#' @param latitude Latitude in decimal degrees.
#' @param pressure Pressure in decibars.
#' @return A data frame containing depth in meters from pressure and adjusted for latitude.
#' @export
calc_depth <- function(latitude, pressure) {
  # Fofonoff NP, Millard RC, Jr. (1983) Algorithms for computation of fundamental properties of seawater. UNESCO Tech Pap Mar Sci 44:53 pp.
  # Assumes water column at 0 deg C, and salinity of 35 PSU
  x <- (sin(latitude / 57.29578))^2 # latitude component of pressure
  p <- pressure * -0.01 # convert pressure to decibars
  g <- 9.780318 * (1.0 + (5.2788e-3 + 2.36e-5 * x) * x) + 1.092e-6 * p # calculation of lat specific gravity
  # calculate ROV depth from pressure
  depth.pressure <- ((((-1.82e-15 * p + 2.279e-10) * p - 2.2512e-5) * p + 9.72659) * p) / g
  return(depth.pressure)
}

#' Convert time from Echoview format to POSIXct
#'
#' @param date Date in Echoview format (Date_M).
#' @param time Time in Echoview format (Time_M).
#' @param tz Time zone.
#' @return A date/time object in POSIXct format.
#' @export
ev2posix <- function(date, time, tz="GMT") {
  # convert date and time string to POSIXct format 20170321 23:14:55
  datetime <- as.POSIXct(paste(date, time), format = "%Y%m%d %H:%M:%OS", tz = tz)
  return(datetime)
}

#' Convert decimal degrees to degrees and decimal minutes
#'
#' @param dd Latitude or longitude in decimal degrees.
#' @return Latitude or longitude in degrees and decimal minutes.
#' @export
dd2decmin <- function(dd, format="latex") {
  deg <- as.integer(dd)
  dd  <- abs(dd)
  decmin <- round((abs(dd) - abs(deg)) * 60, 6)
  if (format == "latex") {
    ddecmin <- paste(deg, "^o^ ", decmin, "'", sep = "")
  } else {
    ddecmin <- paste(deg, "deg ", decmin, "min", sep = "")
  }

  # Add names to output, if present in input
  names(ddecmin) <- names(dd)

  return(ddecmin)
}

#' convert latitude or longitude from SCS format to decimal degrees
#'
#' @param x Latitude or longitude in SCS format.
#' @return Latitude or longitude in decimal degrees.
#' @export
scs2dd <- function (x) {
  if (length(grep("N", x)) > 0) {
    # Remove all non-numeric or decimal characters
    x <- gsub("[^0-9.]", "", x)
    # Parse the remaining characters to extract the latitude
    y <- as.numeric(substr(x, 1, 2)) + signif(as.numeric(substr(x, 3, 9))/60, digits = 6)
  }
  else {
    # Remove all non-numeric or decimal characters
    x <- gsub("[^0-9.]", "", x)
    # Parse the remaining characters to extract the longitude
    y <- -(as.numeric(substr(x, 1, 3)) + signif(as.numeric(substr(x, 4, 10))/60, digits = 6))
  }
  return(y)
}

#' Convert date and time from SCS format to POSIXct
#'
#' @param date Date in SCS format.
#' @param time Time in SCS format.
#' @return A date/time object in POSIXct format.
#' @export
scs2posix <- function(date, time) {
  x <- as.POSIXct(paste(date, time), tz = "GMT", format = "%m/%d/%Y %H:%M:%S")
  return(x)
}

#' Install and load all packages provided from a character vector
#'
#' @param pkgs A character vector containing the names of packages to be loaded and/or installed.
#' @export
load_pkgs <- function(pkgs) {
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new_pkgs) > 0) install.packages(new_pkgs)
  invisible(lapply(pkgs, function(x)
    suppressPackageStartupMessages(library(x, character.only = TRUE))))
}

#' Convert x/y points to latitude/longitude
#'
#' @param x X position, in meters.
#' @param y Y position in meters.
#' @param lat Latitude of known origin.
#' @param lon Longitude of known origin.
#' @param units Distance units in kilometers (default).
#' @return A date/time object in POSIXct format.
#' @export
xy2latlon <- function(x, y, lat, lon, units="km") {
  # Calculate shift in X direction (east positive)
  shift.x <- swfscMisc::destination(lat, lon, 90, x / 1000, units = units)
  df1 <- data.frame(shift.x[grep("lat", names(shift.x))], shift.x[grep("lon", names(shift.x))])
  names(df1) <- c("lat", "lon")
  # Calculate shift in Y direction (north positive)
  shift.y <- swfscMisc::destination(df1$lat, df1$lon, 0, y / 1000, units = units)
  df2 <- data.frame(shift.y[grep("lat", names(shift.y))], shift.y[grep("lon", names(shift.y))])
  names(df2) <- c("lat", "lon")
  # Remove row names
  row.names(df2) <- NULL
  # Return calculated lat/lon data frame
  return(df2)
}

#' Calculate map boundaries from lat/lon input.
#'
#' @param lat Latitude in decimal degrees.
#' @param lon Longitude in decimal degrees.
#' @param pad Percentage that map boundaries are extended beyond the input data.
#' @return A data frame containing the range of \code{lat} and \code{lon}.
#' @export
map_bounds <- function(lat, lon, pad = 0.05) {
  # configure survey plan map
  # determine the lat/lon to add to the data range to achieve the desired frame
  pad.x <- (range(lon)[2] - range(lon)[1]) * pad
  pad.y <- (range(lat)[2] - range(lat)[1]) * pad
  # set limits for desired frame
  range.lat <- c(min(lat) - pad.x, max(lat) + pad.x)
  range.lon <- c(min(lon) - pad.y, max(lon) + pad.y)
  # Return data frame with lat/lon range
  data.frame(range.lat, range.lon)
}

#' Convert numbers less than 9 to words.
#'
#' @param x Number to convert to letters
#' @return If greater than 10, a number, else the number spelled-out in words.
#' @export
num2words <- function(x) {
  ## Function by John Fox found here:
  ## http://tolstoy.newcastle.edu.au/R/help/05/04/2715.html
  ## Tweaks by AJH to add commas and "and"
  helper <- function(x) {
    digits <- rev(strsplit(as.character(x), "")[[1]])
    nDigits <- length(digits)
    if (nDigits == 1) {
      as.vector(ones[digits])
    } else if (nDigits == 2) {
      if (x <= 19) {
        as.vector(teens[digits[1]])
      } else {
        trim(paste(
          tens[digits[2]],
          Recall(as.numeric(digits[1]))
        ))
      }
    } else if (nDigits == 3) {
      trim(paste(
        ones[digits[3]], "hundred and",
        Recall(makeNumber(digits[2:1]))
      ))
    } else {
      nSuffix <- ((nDigits + 2) %/% 3) - 1
      if (nSuffix > length(suffixes)) stop(paste(x, "is too large!"))
      trim(paste(
        Recall(makeNumber(digits[
          nDigits:(3 * nSuffix + 1)
        ])),
        suffixes[nSuffix], ",",
        Recall(makeNumber(digits[(3 * nSuffix):1]))
      ))
    }
  }
  trim <- function(text) {
    # Tidy leading/trailing whitespace, space before comma
    text <- gsub("^\ ", "", gsub("\ *$", "", gsub("\ ,", ",", text)))
    # Clear any trailing " and"
    text <- gsub(" and$", "", text)
    # Clear any trailing comma
    gsub("\ *,$", "", text)
  }
  makeNumber <- function(...) as.numeric(paste(..., collapse = ""))
  # Disable scientific notation
  opts <- options(scipen = 100)
  on.exit(options(opts))
  ones <- c(
    "", "one", "two", "three", "four", "five", "six", "seven",
    "eight", "nine"
  )
  names(ones) <- 0:9
  teens <- c(
    "ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen",
    "sixteen", " seventeen", "eighteen", "nineteen"
  )
  names(teens) <- 0:9
  tens <- c(
    "twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty",
    "ninety"
  )
  names(tens) <- 2:9
  x <- round(x)
  suffixes <- c("thousand", "million", "billion", "trillion")
  if (length(x) > 1) return(trim(sapply(x, helper)))
  helper(x)
}
