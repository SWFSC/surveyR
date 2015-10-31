# functions used by dive_proc.R

calc_cw <- function(pitch,altitude,ROV="HDHV",z=0.06,max.sr=6,smoother=15){
  # smooth altitude
  alt.sm <- as.numeric(ma(altitude,order=smoother))
  # find NAs
  isna <- which(is.na(alt.sm)==TRUE)
  # replace NAs with original values
  alt.sm[isna] <- altitude[isna]

  # smooth pitch
  pitch.sm <- as.numeric(ma(pitch,order=smoother))
  # find NAs
  isna <- which(is.na(pitch.sm)==TRUE)
  # replace NAs with original values
  pitch.sm[isna] <- pitch[isna]

  # Calculate field of view estimates
  if(ROV=="HDHV"){
    # Empirically measured viewing angles for the HDHV ROV
    Alpha	      <- -122.79*z^4 + 132.93*z^3 + 81.984*z^2 - 142.07*z + 46.159# vertical range from 64 deg (full-wide) to 6.6 deg (10X zoom)
    Beta	      <- -297.54*z^4 + 425*z^3 - 15.571*z^2 - 194*z + 74.416 # horiz range from 73.9 deg (full-wide) to 2.26 deg (10X zoom)
    Beta.max    <- 73.95# max horiz viewing angle at full-wide zoom; measured empirically
    Beta.min    <- 2.26 # min horiz viewing angle at 10x zoom; measured empirically
    Beta.spec   <- 85   # max horiz viewing angle at full-wide zoom; manufacturer spec
    Alpha.max   <- 45.9 # max vertical viewing angle at full-wide zoom; measured empirically
    Alpha.min   <- 1.27 # min vertical viewing angle at 10x zoom; measured empirically
    Alpha.spec  <- 64   # max vertical viewing angle at full-wide zoom; manufacturer spec
    Beta.diff   <- 100*(1-(tan((Beta.max/2)*(pi/180))/tan((Beta.spec/2)*(pi/180)))) # difference in area between spec and empirical FOV
    h.offset    <- 0.4127 # altitude of the camera above the frame; added to altitude from DVL to get camera height
  } else {
    # Empirically measured viewing angles for the Phantom DS4 ROV
    Alpha	      <- -122.79*z^4 + 132.93*z^3 + 81.984*z^2 - 142.07*z + 46.159# vertical range from 64 deg (full-wide) to 6.6 deg (10X zoom)
    Beta	      <- -297.54*z^4 + 425*z^3 - 15.571*z^2 - 194*z + 74.416 # horiz range from 73.9 deg (full-wide) to 2.26 deg (10X zoom)
    Beta.max    <- 73.95# max horiz viewing angle at full-wide zoom; measured empirically
    Beta.min    <- 2.26 # min horiz viewing angle at 10x zoom; measured empirically
    Beta.spec   <- 85   # max horiz viewing angle at full-wide zoom; manufacturer spec
    Alpha.max   <- 45.9 # max vertical viewing angle at full-wide zoom; measured empirically
    Alpha.min   <- 1.27 # min vertical viewing angle at 10x zoom; measured empirically
    Alpha.spec  <- 64   # max vertical viewing angle at full-wide zoom; manufacturer spec
    Beta.diff   <- 100*(1-(tan((Beta.max/2)*(pi/180))/tan((Beta.spec/2)*(pi/180)))) # difference in area between spec and empirical FOV
    h.offset    <- 0.4127 # altitude of the camera above the frame; added to altitude from DVL to get camera height
  }
  # Calculate altitude of camera above the seabed (distance between the bottom of the frame and the camera lens) from smoothed altitude
  # Accounts for rotation of the camera module about the tilt tray axis
  alt     <- alt.sm + h.offset + (0.1232 * (sin(pitch.sm + 415.5) * (pi/180)))   # altitude of the camera above the frame

  # calculate the slant range (sr) from the camera to the seabed, using smoothed pitch values
  sr  <- -alt/sin(pitch.sm * (pi/180))

  # replace values when they are negative with zero
  sr[sr < 0] <- NA
  # replace values when they are beyond the theoretical maximum viewing distance
  sr[sr > max.sr] <- NA
  # use linear interpolation to replace missing slant range values
  sr <- as.numeric(na.interp(sr))

  # calculate horizontal width at the center of the camera (or field of view, FOV)
  cw  <- 2 * sr * tan((Beta/2)*(pi/180))

  # output results as a data frame
  data.frame(camera_alt = alt,slant_range = sr,center_width = cw)
}

# convert WinFrog lat/lon (e.g., NDD MM.MMMM,WDDD MM.MMM) to decimal degrees
winfrog2dd <- function(x){
  if(length(grep("N",x))>0){
    as.numeric(substr(x,2,3)) + as.numeric(substr(x,5,11))/60
  } else {
    -(as.numeric(substr(x,2,5)) + as.numeric(substr(x,6,12))/60)
  }
}

# calculate dissolved oxygen saturation from temperature and sal
calc_sat <- function(sal,temp,oxy.conc){
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
  DOSAT_ml      <- exp(a1 + a2*(100/K) + a3*(log(K/100)) + a4*(K/100) + S*(b1 + b2*(K/100)+ b3*((K/100)^2)))
  DOSAT_umoles  <- DOSAT_ml / 0.0224
  oxy.sat       <- (oxy.conc/DOSAT_umoles)*100
  return(oxy.sat)
}

# Calculate depth from pressure in SEAWATER (factors-in gravity and latitude)
calc_depth <- function(latitude,pressure){
  # Fofonoff NP, Millard RC, Jr. (1983) Algorithms for computation of fundamental properties of seawater. UNESCO Tech Pap Mar Sci 44:53 pp.
  # Assumes water column at 0 deg C, and salinity of 35 PSU
  x <- (sin(latitude/57.29578))^2		 # latitude component of pressure
  p <- pressure * -0.01					 # convert pressure to decibars
  g <- 9.780318 * (1.0 + (5.2788e-3 + 2.36e-5 * x) * x) + 1.092e-6 * p  # calculation of lat specific gravity
  # calculate ROV depth from pressure
  depth.pressure 	<- ((((-1.82e-15 * p + 2.279e-10) * p - 2.2512e-5) * p + 9.72659) * p) / g
  return(depth.pressure)
}
