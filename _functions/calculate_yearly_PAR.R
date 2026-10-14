calculate_yearly_PAR = function(image = "", lat = 0, tau = 0.6, 
                        uoc = 0.15, magn.corr = 0){
  force(image)
  fisheye.cx     = image[[2]]
  fisheye.cy     = image[[3]]
  fisheye.cr     = image[[4]]
  
  deg2rad        = pi / 180
  rad2deg        = 180 / pi
  solar.constant = 1360 * 0.51 * 4.6           ## PPFD (umol m-2 s-1)
  uMol2Mol       = 2 * 60 / 1000000
  day2angle      = (pi / 180) * (360 / 365)
  d              = 1:365 
  rad.year        = matrix(0, 365, 5)
  rad.year[,1]    = d 
  
  ## diffuse factor
  fractions = calculate_gap_fraction(image = image)
  a         = deg2rad * 1:89
  d05       = deg2rad * 0.5
  Atot      = sin(a[89] + d05) - sin(a[1] - d05)
  Aa        = (sin(a + d05) - sin(a - d05)) * sin(a)
  ISF       = sum(fractions * Aa / Atot)
  
  ## latitude angle (psi)
  if (lat >  89.99) lat =  89.99
  if (lat < -89.99) lat = -89.99
  psi     = lat * deg2rad
  sin.psi = sin(psi)
  cos.psi = cos(psi)
  
  for(k in 1:365){
    
    ## solar declination (delta) - method see solar.declination()
    day.angle    = d[k] * day2angle
    sin.delta   = 0.3985 * sin(4.869 + day.angle + 0.03345 * sin(6.224 + day.angle))
    cos.delta   = sqrt(1 - sin.delta*sin.delta)
    R_out       = solar.constant + 0.034 * cos(day.angle)   #Radiation at outer atmosphere
    
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
      
      sin.angle = sin(magn.corr*deg2rad)
      cos.angle = cos(magn.corr*deg2rad)
      sin.beta_c = (sin.beta*cos.angle + cos.beta*sin.angle)   ##sin(a+b) = sin(a)*cos(b) + cos(a)*sin(b)
      cos.beta_c = (cos.beta*cos.angle - sin.beta*sin.angle)   ##cos(a+b) = cos(a)*cos(b) - sin(a)*sin(b)
      
      ## use track points only if alfa > 0 (sun above horizon)
      if(alfa > 0){
        ## nr of airmasses direct light has to pass through
        M        = sqrt(1229 + (614 * sin.alfa)^2) - 614 * sin.alfa 
        
        ## Transmission of direct light through airmasses
        trans    = tau^M
        
        ## direct light
        direct  = solar.constant * trans * sin.alfa
        
        ## diffuse light
        ## used is the Uniform Overcast Sky (UOC)
        diffuse = direct * uoc                             #Uniform Overcast Sky (15% of direct light becomes diffuse
        
        ## other methods not yet implemented
        #diffuse = R_out * (0.271 - 0.294*Trans)*sin.alfa  #Clear not dust free sky        
        #diffuse = direct * sin.alfa * (1 + 2*sin.alfa)/3) #SOC 
        
        ## pixel location
        xp   = fisheye.cx - round(sin.beta_c * (90-rad2deg*alfa) * fisheye.cr/90)
        yp   = fisheye.cy + round(cos.beta_c * (90-rad2deg*alfa) * fisheye.cr/90)
        
        i2 = i-sun.up+1
        if(length(d) == 1){
          Radiation[i2,1] = hour
          Radiation[i2,2] = hour + TimeOffset
        }
        radiation[i2,3] = radiation[i2,3] + direct
        radiation[i2,4] = radiation[i2,4] + diffuse
        radiation[i2,5] = radiation[i2,5] + direct * image[[1]][yp,xp]
        radiation[i2,6] = radiation[i2,6] + diffuse * ISF
      }
    }
    rad.year[k,2] = sum(radiation[,3]) * uMol2Mol
    rad.year[k,3] = sum(radiation[,4]) * uMol2Mol
    rad.year[k,4] = sum(radiation[,5]) * uMol2Mol
    rad.year[k,5] = sum(radiation[,6]) * uMol2Mol
  } 
  colnames(rad.year) = c("Day", "DirectAbove", "DiffuseAbove", "DirectUnder", "DiffuseUnder")
  return(rad.year)
}
