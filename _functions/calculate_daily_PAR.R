calculate_daily_PAR = function(image = "", lat = 0, lon = 0, time.zone = 0, 
                       d = 150, tau = 0.6, uoc = 0.15, 
                       magn.corr = 0, draw.tracks = F, full.day = T){
  force(image)
  fisheye.cx    = image[[2]]
  fisheye.cy    = image[[3]]
  fisheye.cr    = image[[4]]
  
  deg2rad        = pi / 180
  rad2deg        = 180 / pi
  solar.constant = 1360 * 0.51 * 4.6            ## PPFD (umol m-2 s-1)
  uMol2Mol       = 2 * 60 / 1000000             ## PAR calculated every 2 minutes + conversion from umol to mol
  day2angle      = (pi / 180) * (360 / 365)
  
  ## diffuse factor
  fractions = CalcGapFractions(image = image)
  a         = deg2rad*1:89
  d05       = deg2rad*0.5
  Atot      = sin(a[89] + d05) - sin(a[1] - d05)
  Aa        = (sin(a + d05) - sin(a - d05))*sin(a)
  ISF       = sum(fractions*Aa/Atot)
  
  ## latitude angle (psi)
  if (lat >  89.99) lat =  89.99
  if (lat < -89.99) lat = -89.99
  psi     = lat * deg2rad
  sin.psi = sin(psi)
  cos.psi = cos(psi)
  
  rad.day  = matrix(0, length(d) , 6)     
  
  for(k in 1:length(d)){
    
    ## solar declination (delta)
    day.angle    = d[k] * day2angle
    sin.delta   = 0.3985 * sin(4.869 + day.angle + 0.03345 * sin(6.224 + day.angle))
    cos.delta   = sqrt(1 - sin.delta*sin.delta)
    R_out       = solar.constant + 0.034 * cos(day.angle)   #Radiation at outer atmosphere
    
    ## Equation of time
    if(length(d) == 1){
      EoT = 229.18*(0.000075+0.001868*cos(day.angle)-0.032077*sin(day.angle)-0.014615*cos(2*day.angle)-0.040849*sin(2*day.angle))
      time.offset = (EoT + 15*(lon - (15 * time.zone)))/60  # in hours
    }
    
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
    sun.up   = (round(12 - (eta*12/pi)))*30 #- 5 #{add a safety 10 min}
    sun.down = (round(12 + (eta*12/pi)))*30 #+ 5 #{add a safety 10 min}
    
    nrows = 1+sun.down-sun.up
    radiation = matrix(0,nrows,6) 
    
    for(i in sun.up:sun.down){
      hour = i/30                    # 2 minute steps - size of solar disc in sky
      ## hour angle (eta)
      sin.eta  = sin(pi*(hour-12)/12)
      cos.eta  = cos(pi*(hour-12)/12)
      
      ## Solar angle(alfa) - method see SunLocation
      sin.alfa = sin.psi*sin.delta + cos.psi*cos.delta*cos.eta
      cos.alfa = sqrt(1-(sin.alfa*sin.alfa))
      alfa     = atan(sin.alfa/cos.alfa);
      
      ## solar azimuth (beta)
      sin.beta = -cos.delta*sin.eta/cos.alfa
      cos.beta = -(sin.delta-sin.psi*sin.alfa)/(cos.psi*cos.alfa)
      
      sin.angle = sin(magn.corr * deg2rad)
      cos.angle = cos(magn.corr * deg2rad)
      sin.beta_c = (sin.beta*cos.angle + cos.beta*sin.angle)   ##sin(a+b) = sin(a)*cos(b) + cos(a)*sin(b)
      cos.beta_c = (cos.beta*cos.angle - sin.beta*sin.angle)   ##cos(a+b) = cos(a)*cos(b) - sin(a)*sin(b)
      
      ## use track points only if alfa > 0 (sun above horizon)
      if(alfa > 0){
        ## airmasses
        M        = sqrt(1229 + (614 * sin.alfa)^2) - 614 * sin.alfa 
        
        ## Transmission of direct light through airmasses
        trans    = tau^M
        
        ## direct light
        direct   = solar.constant * trans * sin.alfa
        
        ## diffuse light as Uniform Overcast Sky (UOC)
        # in which 15% of direct light becomes diffuse
        diffuse  = direct * uoc                            #Uniform Overcast Sky (15% of direct light becomes diffuse
        
        ## other options not yet implemented
        #diffuse = R_out * (0.271 - 0.294*Trans)*sin.alfa  #Clear not dust free sky        
        #diffuse = direct * sin.alfa * (1 + 2*sin.alfa)/3) #SOC 
        
        ## pixel location
        xp   = fisheye.cx - round(sin.beta_c * (90-rad2deg*alfa) * fisheye.cr/90)
        yp   = fisheye.cy + round(cos.beta_c * (90-rad2deg*alfa) * fisheye.cr/90)  # on image
        yp2  = fisheye.cy - round(cos.beta_c * (90-rad2deg*alfa) * fisheye.cr/90)  # on plot!
        
        if(draw.tracks == T) points(xp, yp2, cex = 0.1, col = "red")
        
        i2 = i-sun.up+1
        if(length(d) == 1){
          radiation[i2,1] = hour
          radiation[i2,2] = hour + time.offset
        }
        radiation[i2,3] = radiation[i2,3] + direct
        radiation[i2,4] = radiation[i2,4] + diffuse
        radiation[i2,5] = radiation[i2,5] + direct * image[[1]][yp,xp]
        radiation[i2,6] = radiation[i2,6] + diffuse * ISF
      }
    }
    rad.day[k,1] = sum(radiation[,3]) * uMol2Mol / length(d)
    rad.day[k,2] = sum(radiation[,4]) * uMol2Mol / length(d)
    rad.day[k,3] = sum(radiation[,5]) * uMol2Mol / length(d)
    rad.day[k,4] = sum(radiation[,6]) * uMol2Mol / length(d)
  } 
  if(length(d) == 1 & full.day == T){ 
    colnames(radiation) = c("Solar Time","Location Time","DirectAbove", "DiffuseAbove", "DirectUnder", "DiffuseUnder")
    return(radiation)
  } else {
    rad.out = c(sum(rad.day[,1]), sum(rad.day[,2]), sum(rad.day[,3]), sum(rad.day[,4])) 
    names(rad.out) = c("DirectAbove", "DiffuseAbove", "DirectUnder", "DiffuseUnder")
    return(rad.out) 
  }
}