calculate_gap_fraction = function(image = ""){
  deg2rad = pi / 180
  gap.fractions = matrix(0,89)
  steps         = 1:360
  fisheye.cx    = image[[2]]
  fisheye.cy    = image[[3]]
  fisheye.cr    = image[[4]]
  
  for (i in 1:89){
    x = round(fisheye.cx + cos(steps*deg2rad)*i*fisheye.cr/90,0)
    y = round(fisheye.cy + sin(steps*deg2rad)*i*fisheye.cr/90,0)
    for(j in 1:360) gap.fractions[i] = gap.fractions[i] + image[[1]][y[j],x[j]]
  }
  return(gap.fractions/360)
}
