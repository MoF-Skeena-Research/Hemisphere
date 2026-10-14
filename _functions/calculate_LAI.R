calculate_LAI = function(fractions, width = 6){
  if(dim(fractions)[1] != 89) return("incorrect dimensions for fractions") 
  deg2rad = pi/180
  angle    = c(7, 23, 38, 53, 68)                  # angles of LAI2000
  wi       = c(0.034, 0.104, 0.160, 0.218, 0.494)  # weights given by Licor canopy analyzer manual
  
  if(width < 0 | width > 6) width = 6
  T = rep(0,5)
  for(i in -width:width) T = T + gap.fractions[angle + i]
  T = T/(2*width +1)
  
  LAI      = 2 * sum(-log(T) * wi * cos(angle*deg2rad))
  return(LAI)
}