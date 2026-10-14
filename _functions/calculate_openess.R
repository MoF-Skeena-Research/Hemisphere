calculate_openess = function(fractions){
  deg2rad        = pi / 180
  a              = deg2rad * 1:89
  d05            = deg2rad * 0.5
  Atot           = sin(a[89] + d05) - sin(a[1] - d05)
  Aa             = sin(a + d05) - sin(a - d05)
  canopy.openess = sum(fractions*Aa/Atot)
  return(canopy.openess)
}