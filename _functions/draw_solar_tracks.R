draw_solar_tracks = function(image ="", lat = 0, lon = 0, time.zone = 0, 
                           d = 150, magn.corr = 0, sun.location = F, h = 0){
  fisheye.cx   = image[[2]]
  fisheye.cy   = image[[3]]
  fisheye.cr   = image[[4]]
  
  deg2rad      = pi / 180
  rad2deg      = 180 / pi
  day2angle    = (pi / 180) * (360 / 365)
  
  ## latitude angle (psi)
  if (lat >  89.99) lat =  89.99
  if (lat < -89.99) lat = -89.99
  psi     = lat * deg2rad
  sin.psi = sin(psi)
  cos.psi = cos(psi)
  
  for(k in 1:length(d)){
    ## solar declination (delta)
    day.angle    = d[k] * day2angle
    sin.delta   = 0.3985 * sin(4.869 + day.angle + 0.03345 * sin(6.224 + day.angle))
    cos.delta   = sqrt(1 - sin.delta*sin.delta)
    
    ## Solar day length
    ## hour angle (eta) for sun.up and sun.down
    cos.eta = -(sin.psi/cos.psi)*(sin.delta/cos.delta)
    if (cos.eta > 1) {
      cos.eta = 1
      eta = pi
    } else if(cos.eta < -1){
      cos.eta = -1
      eta = pi
    } else if(cos.eta == 0){
      eta = 0.5 * pi
    } else {
      eta = atan(sqrt(1 - cos.eta*cos.eta)/cos.eta)
    }
    if(eta < 0) eta = eta + pi
    sun.up   = (round(12 - (eta*12/pi)))*30 - 10 #{add a safety 20 min}
    sun.down = (round(12 + (eta*12/pi)))*30 + 10
    
    for(i in sun.up:sun.down){
      hour = i/30                    # 2 minute steps - size of solar disc in sky
      
      ## hour angle (eta)
      sin.eta  = sin(pi*(hour-12)/12)  #pi/12 = 2pi/24!
      cos.eta  = cos(pi*(hour-12)/12)  #pi/12 = 2pi/24!
      
      ## Solar angle(alfa)
      sin.alfa = sin.psi*sin.delta + cos.psi*cos.delta*cos.eta
      cos.alfa = sqrt(1-(sin.alfa*sin.alfa))
      alfa     = atan(sin.alfa/cos.alfa);
      
      ## solar azimuth (beta)
      sin.beta = -cos.delta*sin.eta/cos.alfa
      cos.beta = -(sin.delta-sin.psi*sin.alfa)/(cos.psi*cos.alfa)
      
      sin.angle = sin(magn.corr*deg2rad)
      cos.angle = cos(magn.corr*deg2rad)
      sin.beta_c = (sin.beta*cos.angle + cos.beta*sin.angle)   ##sin(a+b) = sin(a)*cos(b) + cos(a)*sin(b)
      cos.beta_c = (cos.beta*cos.angle - sin.beta*sin.angle)   ##cos(a+b) = cos(a)*cos(b) - sin(a)*sin(b)
      
      ## draw track points only if alfa > 0 (sun above horizon)
      if(alfa > 0){
        xp   = fisheye.cx - (sin.beta_c * (90-rad2deg*alfa) * fisheye.cr/90)
        yp   = fisheye.cy - (cos.beta_c * (90-rad2deg*alfa) * fisheye.cr/90)
        points(xp, yp, cex = 0.1, col = "red")
      }
    }
    if(sun.location == T){
      ## Equation of time
      EoT = 229.18*(0.000075+0.001868*cos(day.angle)-0.032077*sin(day.angle)-0.014615*cos(2*day.angle)-0.040849*sin(2*day.angle))
      hour = h + (EoT + 15*(lon - (15 * time.zone)))/60  # in hours
      
      ## hour angle (eta)
      sin.eta  = sin(pi*(hour-12)/12)  #pi/12 = 2pi/24!
      cos.eta  = cos(pi*(hour-12)/12)  #pi/12 = 2pi/24!
      
      ## Solar angle(alfa)
      sin.alfa = sin.psi*sin.delta + cos.psi*cos.delta*cos.eta
      cos.alfa = sqrt(1-(sin.alfa*sin.alfa))
      alfa     = atan(sin.alfa/cos.alfa);
      
      ## solar azimuth (beta)
      sin.beta = -cos.delta*sin.eta/cos.alfa
      cos.beta = -(sin.delta-sin.psi*sin.alfa)/(cos.psi*cos.alfa)
      
      sin.angle = sin(magn.corr*deg2rad)
      cos.angle = cos(magn.corr*deg2rad)
      sin.beta_c = (sin.beta*cos.angle + cos.beta*sin.angle)   ##sin(a+b) = sin(a)*cos(b) + cos(a)*sin(b)
      cos.beta_c = (cos.beta*cos.angle - sin.beta*sin.angle)   ##cos(a+b) = cos(a)*cos(b) - sin(a)*sin(b)
      if(alfa > 0){
        xp   = fisheye.cx - (sin.beta_c * (90-rad2deg*alfa) * fisheye.cr/90)
        yp   = fisheye.cy - (cos.beta_c * (90-rad2deg*alfa) * fisheye.cr/90)
        points(xp, yp, cex = 1, pch = 21, col = "yellow", bg = "yellow")
      } 
    }
  } 
}