# Soil Water model function for automatic calibration#
SWMsim <- function (         # Model parameter definition
  rE = -0.2,           # extinction factor
  alpha = 3,          # interception threshold parameter, default set to 4.5 mm (Simunek et al, 2008)
  Smax = 85,             # maximum measured soil moisture content upper horizon
  Ic = 52,              # max Infiltration capacity, if exceeding surface runoff Qs generation
  ks1 = 18,						# saturated hydraulic conductivity upper soil horizon
  ks2 = 18,					# saturated hydraulic conductivity lower horizon
  ks3 = 20,             # saturated hydraulic conductivity deepest horizon
  GWmax = 90,					# maximum measured soil moisture content lower horizon		
  Lmax = 280,           # maximum measured soil moisture content deeper horizon
  g1 = 3,            # nonlinear scaling parameter of upper horizon, if g1 = 1, linear case	
  g2 = 4,       # nonlinear scaling parameter of lower horizon, if g2 = 1, linear case
  g3 = 5.5,             # nonlinear scaling parameter of deeper horizon, if g2 = 1, linear case
  PF_Scale = 0.80,     # preferential flowpath parameter
  INTp = 1,          #passive interception storage mixing volume 
  stoSp = 20,       # passive upper storage mixing volume 
  gwSp = 40,         # passive lower storage mixing volume 
  lowSp = 100,          # passive deep storage mixing volume
  k = 0.9,           # seasonality factor CG model
  x = 0.75,         # water vapour mixing ratio CG model
  beta = 0          # extinction factor for root density function
)             				

{  ##### BEGIN OF FUNCTION BODY #####
  
  nrun <- nrow (inp)     # model time step equal to input
  
  tsCount <- rep(NA,nrun)  
  ## Dynamic meteorological and interception variables
  
  PN <- rep(NA,nrun)	            # net precipitation minus interception threshold
  I <- rep(NA,nrun)               # interception = 1 - PN
  D <- rep(NA,nrun)               # interception threshold 
  IMax <- rep(NA,nrun)            # time-varying maximum canopy storage capacity - AJN
  IAvail <- rep(-7,nrun)            # time-varying available canopy storage capacity - AJN
  SCF <- rep(NA,nrun)             # surface cover fraction 
  Ei	<- rep(NA,nrun)							# interception evaporation
  Tr <- rep(NA,nrun)              # transpiration rate
  Es	<- rep(NA,nrun)							# soil evaporation
  Ep	<- rep(NA,nrun)							# potential evaporation
  Tp	<- rep(NA,nrun)							# potential transpiration
  Th <- rep(NA, nrun)              # throughfall
  Tr_Upper <- rep(NA, nrun) 
  Tr_Lower <- rep(NA, nrun)        # Transpiration taken from the lower store
  Tr_Deep <- rep(NA, nrun) 
  
  ## Iteratively filled storage and flux variables
  STO <- rep(NA,nrun)					# upper soil horizon storage in mm
  GW <- rep(NA,nrun)					# lower soil horizon storage in mm
  Sdeep <- rep(NA,nrun)       # deep soil horizon storage in mm
  GWflow <- rep(NA,nrun)			# recharge from lower soil horizon in mm
  Recharge <- rep(NA,nrun)		# recharge loss from Sdeep
  Perc <- rep(NA,nrun)        # vertical flux from STO into GW in mm
  Qs  <- rep(NA,nrun)         # Overland flow if infiltration capacity parameter Ic is exceeded in mm
  #API  <- rep(NA,nrun)         # API vector
  Pref_Flow<- rep(NA, nrun)   # Preferential flow occuring in larger Net Precipitation events
  
  ## Isotope module
  upCP_D <- rep(NA,nrun)					# Precipitation D-concentration
  upCQ_D <- rep(NA,nrun)					# simulated D discharge concentration in permil
  gwCQ_D <- rep(NA,nrun)				# simulated D discharge concentration in permil
  lowCQ_D <- rep(NA,nrun)
  upSTO_D <- rep(NA,nrun)					# storage volume D-concentration
  upCSTO_D <- rep(NA,nrun)				# analog for storage D-concentration
  gwSTO_D <- rep(NA,nrun)	
  gwCSTO_D <- rep(NA,nrun)				# analog for storage D-concentration
  lowSTO_D <- rep(NA,nrun)	
  lowCSTO_D <- rep(NA,nrun)
  Int_D <- rep(NA,nrun)	
  Int_CD <- rep(NA,nrun)	
  fInt_CD <- rep(NA,nrun)
  fupCSTO_D <- rep(NA,nrun)
  WV_D <- rep(NA,nrun)
  
  # Craig-Gordon model
  Idl_H <- rep(NA,nrun) #isotopic content of the residual liquid (sometimes referred to as dS)
  da_H <- rep(NA,nrun) # isotopic content of the ambient athmospheric vapor (atmospheric moisture)
  IdE_H <- rep(NA,nrun) # isotopic content of the evaporating flux
  alphae_H <- rep(NA,nrun)  # fractionation factor at equilibrium (Rliq/Rvap)
  epse_H <- rep(NA,nrun)   # equilibrium fractionation factor (expressed in permil) 
  Tk <- rep(NA,nrun) # Temperature [Kelvin]
  epsk_H <- rep(NA,nrun) # kinetic fractionation factor (expressed in permil) 
  m_H <- rep(NA,nrun)
  dstar_H <- rep(NA,nrun)
  Sdl_H <- rep(NA,nrun) #isotopic content of the residual liquid (sometimes referred to as dS)
  SdE_H <- rep(NA,nrun) # isotopic content of the evaporating flux

  # root fraction and Root Density Distribution Function
#
#  if (beta == 0) {
#    r_Upper <- 0.1
#    r_Lower <- 0.2
#    r_Deep <- 0.7
#  } else {
#    r_Upper <- (1 - exp(-beta * 0.1)) / (1 - exp(-beta * 1.0))
#    r_Lower <- (exp(-beta * 0.1) - exp(-beta * 0.3)) / (1 - exp(-beta * 1.0))
#    r_Deep <- 1 - r_Upper - r_Lower
#  }
 
  # Root Water Withdrawal Efficiency represents the capacity of roots to extract water from a specific soil depth,
  # modeled as an exponential decay function of depth, where a higher decay rate (beta) concentrates root activity near the surface. 

  r_Upper <- exp(-beta * 0.05)
  r_Lower <- exp(-beta * 0.2)
  r_Deep  <- exp(-beta * 0.65)


  #~~~~~~~~~~~~~~ calculate iterative dependent time series ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  for (i in 1:nrun) {					# the loop iterates through the length of the input time series
    
    tsCount[i] <- i  # count time steps 
    #~~~~~~~~~~~~~~ calculate dynamic dependent time series ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Craig-Gordon model
    # get equilibrium fractionation factors from Horita and Wesolowski
    Tk[i] <- AT[i] + 273.15
    
    alphae_H[i] <- exp(1/1000.*(1158.8*Tk[i]^3/10^9-1620.1*Tk[i]^2/10^6+
                                  794.84*Tk[i]/10^3-161.04+2.9992*10^9./Tk[i]^3)) #Horita and Wesolowski (1994) 
    
    epse_H[i] <- (alphae_H[i] - 1)*1000 #permil notation
    
    # get kinetic fractionation factors
    epsk_H[i] <- 0.9755 *(1-0.9755)*1000*(1-RH[i]) # n=0.9755 value from Merlivat (1978) #permil notation
    
    # get atmospheric composition from precipitation-equilibrium assumption (Gibson et al., 2008)
    #  k <- 1 #seasonality factor 
    da_H[i] <- (-58.7 - k*epse_H[i])/(1+k*epse_H[i]*10^-3) 
    
    # compute the useful variables m and dstar ('enrichment slope and limiting isotopic composition)
    m_H[i] <- (RH[i]-10^-3*(epsk_H[i] +epse_H[i]/alphae_H[i]))/(1- RH[i] +10^-3*epsk_H[i]) #'enrichment slope'' (Gibson et al.(2016))
    
    dstar_H[i] <- (RH[i] *da_H[i] +epsk_H[i] +epse_H[i]/alphae_H[i])/
      (RH[i]-10^-3*(epsk_H[i] +epse_H[i]/alphae_H[i])) #this is A/B in Gonfiantini 1986
    
    ## storage initialization
    
    if (i==1) {							# conditional: filling of storage variables with initial volumes...
    #  API[i] <- 0.995*2      # set the first api value to the last rainfall 
      I[i] <- (alpha*LAI[i]*0.5) # interception store initialised half full for current LAI??
      STO[i] <- 20.34						# measured start value				
      GW[i] <- 36.58             # measured start value
      Sdeep[i] <- 86.31            # measured start value
    }
    else {							# ... or i-1 storage fillings
    #  API[i] <- API[i-1]
      I[i] <- I[i-1]
      STO[i] <- STO[i-1]
      GW[i] <- GW[i-1]
      Sdeep[i] <- Sdeep[i-1]
    }
    
  #  API[i] <- 0.995*API[i] + P[i]        # daily decay constant of 0.95 by Hill et al. (2014)
    
    if ( is.na(P[i]) || is.na(PET[i]) || is.na(LAI[i]))   # no input data NAs allowed! - AJN updated; only need P and PET here now, but LAI added too...
    {
      P[i] <- 0
      PET[i] <- 0
      LAI[i] <- 0
    }
    
    # taken outside else statement so NA not assigned to Ep and TP - PN also defined later now... 
    SCF[i] <- 1 - exp(rE * LAI[i])     # surface cover fraction, rExtinct = -0.463 Rutter (1972)
    
    Tp[i] <- (PET[i] * SCF[i])            # fraction potential canopy ET (Ei plus Tr when I = 0) AJN
    Ep[i] <- (PET[i] * (1-SCF[i]))        # fraction potential soil evaporation only AJN
    
    #~~ Interception, simplified Simunek et al (2008) model ##~~#~~#~~
    # Available canopy storage as maximum storage for the current LAI (alpha * LAI[i]) minus the depth already filled
 
    IMax[i] <- (alpha * LAI[i])

    IAvail <- IMax[i] - I[i]

    # Check if change in LAI reduces available storage below current depth in storage (i.e. IAvail < 0)
    if(IAvail<0){

      # Excess water due to change in LAI added to throughfall (MIXING!)
      Th[i] <- I[i] - IMax[i]

      # Then no available interception storage, interception store is set to max for this timestep, no interception
      IAvail <- 0
      I[i] <- IMax[i]
      D[i] <- 0

    } else if(IAvail == 0){

      # If there's exactly 0 available storage - can happen if LAI constant and PET = 0...
      Th[i] <- 0        # No throughfall from max storage capacity change
      D[i] <- 0         # No interception

      I[i] <- IMax[i]   # Interception store is at maximum (already was the case!)

    } else{

      # Calculate the interception based on available storage
      D[i] <- (IAvail) * (1 - (1 / (1 + ((SCF[i] * P[i])/(IAvail)))))

      # Update interception store
      I[i] <- I[i] + D[i]

      # No throughfall
      Th[i] <- 0

    }

    # PN will be P plus excess storage Th if LAI reduces available storage below depth already stored (MIXING IMPLICATIONS)
    # PN will equal P minus interception D if there is available storage, and Th will be 0.
    
    PN[i] <- (P[i] - D[i]) + Th[i]            # net Precip input to storages 
    
    ## water balance calculations
    if(Tp[i]<= I[i]){
      Ei[i]<- Tp[i]           # If sufficient water in store to meet full potential canopy ET (Tr), full volume taken from canopy
      I[i]<- I[i] - Ei[i]     # Interception store reduced by evaporated volume
      Tp[i]<- 0}              # Potential canopy ET now 0
    else{
      Ei[i]<- I[i]            # If potential canopy ET greater than canopy storage all of interception store evaporated
      I[i]<- 0                # Interception store therefore empty
      Tp[i]<- Tp[i] - Ei[i]   # Potential canopy ET updated to reflect used volume - remaining volume available for transpiration
    }
    
    # Hortonion OLF
    if (PN[i] > Ic) {          #surface runoff Qs if Ic is exceeded
      Qs[i] <- PN[i] - Ic
      PN[i] <- Ic}
    else{
      Qs[i] <- 0
      PN[i] <- PN[i]}
    
    if(PN[i] > 4.2){                        # In larger net precipitation events (>P90) preferential flow allowed
      Pref_Flow[i]<- PN[i] * PF_Scale     # Fraction of net precipitation constituting preferential flow determined by calibrated parameter - added below to Recharge[i]
      PN[i] <- PN[i] * (1 - PF_Scale)}    # Net precipitation to upper soil box reduced accordingly
    else{
      Pref_Flow[i]<- 0
      PN[i]<- PN[i] }
    
    STO[i] <- STO[i] + PN[i]				# upper soil horizon storage plus net precipitation
    
    if(is.na(STO[i]))
    {STO[i] <- 0}
    
    if (STO[i] < Smax) {

      abs(Tr_Upper[i] <- (STO[i] / Smax) * Tp[i] * r_Upper) # Transpiration from upper store calculated

      Tr_Upper[i] <- min(Tr_Upper[i], STO[i])

      STO[i] <- STO[i] - Tr_Upper[i]}              # Upper store reduced accordingly
####      Tp[i]<- Tp[i] - Tr_Upper[i]                # Potential transpiration reduced accordingly
    else {                                        # If store equalling or greater than Smax then all potenital transpiration used.
 
      Tr_Upper[i] <- Tp[i] * r_Upper              # Cannot use (STO[i] / Smax) as this will create number >1 and so transpiration would be greater than potential
      
      Tr_Upper[i] <- min(Tr_Upper[i], STO[i])

      STO[i]<- STO[i] - Tr_Upper[i]
####      Tp[i]<- 0
         }   
    
    if (is.na(Tr_Upper[i]))                       # Code to ensure no NA or negative transpiration values
    {Tr_Upper[i] <- 0}
    if (Tr_Upper[i] < 0) 
    { Tr_Upper[i] <- Tr_Upper[i] * -1} 
    else {Tr_Upper[i] <- Tr_Upper[i]}
    
    if(Ep[i] == 0){                     # If no potential evaporation voulme, soil evaporation cannot occur
      Es[i] <- 0}
    else{
      if(STO[i] > Ep[i]){                 # If greater soil storage than potential evporation:
        Es[i]<- (STO[i] / Smax) * Ep[i]   # Soil evaporation is calculated using full potential evaporation volume as multiplier
        Es[i]<- min(Es[i],STO[i])
        STO[i] <- STO[i] - Es[i]
        Ep[i]<- Ep[i] - Es[i]}	
      else{                               # If less soil storage than potential evaporation:
        Es[i]<- (STO[i] / Smax) * STO[i]  # Soil evaporation is calculated using a potential evaporation volume reduced to be identical to the volume of water in storage (therefore == STO[i])
        STO[i] <- STO[i] - Es[i]
        Ep[i]<- Ep[i] - Es[i]}
   #   if (API[i] < 5)   # only Es if median API is exceeded, else Es becomes Tr
   #   {Tr_Upper[i] <- Tr_Upper[i] + Es[i]
   #   Es[i] <- 0}
   #   else {Tr_Upper[i] <- Tr_Upper[i]
    #  Es[i] <- Es[i]}
    }
    
    if (PN[i] > 0)
    { Perc[i] <- (ks1 * (STO[i] / Smax) ^ g1)				# nonlinear calculation of percolation to lower soil storage 
    Perc[i] <- min(Perc[i], STO[i])
 
    STO[i] <- STO[i] - Perc[i]	}
    else {Perc[i] <- 0}
    
    if(is.na(GW[i]))
    {GW[i] <- 0}
    if (is.na(Perc[i]))
    {Perc[i] <- 0}
    
    GW[i] <- GW[i] + Perc[i] + Pref_Flow[i]       # fill lower soil horizon
    
    if(Tp[i] <= 0){                                    # If all potential transpiration taken from upper store no volume can be taken from lower store (<= to be on safe side)
      Tr_Lower[i]<- 0}
    else{                                              
      if (GW[i] < GWmax) {

        abs(Tr_Lower[i] <- (GW[i] / GWmax) * r_Lower * (Tp[i] - Tr_Upper[i]))  # If potenital transpiration volume available the lower store transpiration volume calculated as per upper store
 
        Tr_Lower[i]<- min(Tr_Lower[i], GW[i])
 
        GW[i] <- GW[i] - Tr_Lower[i]}
####        Tp[i]<- Tp[i] - Tr_Lower[i]                     
      else {                                           
        
        Tr_Lower[i] <- r_Lower * (Tp[i] - Tr_Upper[i]) # If GW[i] exceeds GWmax then full potential volume removed (using GW[i]/ GWmax would create value > potenital)
 
        Tr_Lower[i]<- min(Tr_Lower[i], GW[i])
 
        GW[i]<- GW[i] - Tr_Lower[i]
####        Tp[i]<- 0
           }   
      
      if (is.na(Tr_Lower[i]))                          # Code to ensure no NA or negative transpiration values
      {Tr_Lower[i] <- 0}
      if (Tr_Lower[i] < 0) 
      { Tr_Lower[i] <- Tr_Lower[i] * -1} 
      else {Tr_Lower[i] <- Tr_Lower[i]}
    }
    Tr[i]<- Tr_Upper[i] + Tr_Lower[i]                 # Total transpiration volume 
    
    if ( Perc[i] > 0 ) {
      GWflow[i] <- ks2 * (GW[i] / GWmax) ^ g2					# nonlinear recharge calculation based on soil storage volume
      GW[i] <- GW[i] - GWflow[i]}
    
    if ( is.na(GWflow[i]))  						
    {
      GWflow[i] <- 0
    }
    else {GWflow[i] <-  GWflow[i]} 			# simulated percolation
    
    if(is.na(Sdeep[i]))
    {Sdeep[i] <- 0}
    
    Sdeep[i] <- Sdeep[i] + GWflow[i]        # fill lower soil horizon
    
    if(Tp[i] <= 0){                                    # If all potential transpiration taken from upper store no volume can be taken from lower store (<= to be on safe side)
      Tr_Deep[i]<- 0}
    else{                                              
      if (Sdeep[i] < Lmax) {

        abs(Tr_Deep[i] <- (Sdeep[i] / Lmax) * r_Deep * (Tp[i] - Tr_Upper[i] - Tr_Lower[i]))  # If potenital transpiration volume available the lower store transpiration volume calculated as per upper store
        
        Tr_Deep[i]<- min(Tr_Deep[i], Sdeep[i])

        Sdeep[i] <- Sdeep[i] - Tr_Deep[i]}
####        Tp[i]<- Tp[i] - Tr_Deep[i]}                     
      else {     
                                      
        Tr_Deep[i] <- r_Deep * (Tp[i] - Tr_Upper[i] - Tr_Lower[i])  # If GW[i] exceeds GWmax then full potential volume removed (using GW[i]/ GWmax would create value > potenital)
 
        Tr_Deep[i]<- min(Tr_Deep[i], Sdeep[i])

        Sdeep[i]<- Sdeep[i] - Tr_Deep[i]
####        Tp[i]<- 0
           }   
      
      if (is.na(Tr_Deep[i]))                          # Code to ensure no NA or negative transpiration values
      {Tr_Deep[i] <- 0}
      if (Tr_Deep[i] < 0) 
      { Tr_Deep[i] <- Tr_Deep[i] * -1} 
      else {Tr_Deep[i] <- Tr_Deep[i]}
    }
    Tr[i]<- Tr_Upper[i] + Tr_Lower[i] + Tr_Deep[i]                # Total transpiration volume 
    
    if ( GWflow[i] > 0 ) {
      Recharge[i] <- ks3 * (Sdeep[i] / Lmax) ^ g3					# nonlinear recharge calculation based on soil storage volume

      Recharge[i] <- min(Recharge[i], Sdeep[i])

      Sdeep[i] <- Sdeep[i] - Recharge[i]}
    
    if ( is.na(Recharge[i]))  						
    {
      Recharge[i] <- 0
    }
    else {Recharge[i] <- Recharge[i]} 			# simulated recharge
    
    #### ISOTOPE module ####
    
    if (P[i]>0 )							# uniform isotope composition over plot 
    {			upCP_D[i] <- P_D[i]			 }
    else  {		upCP_D[i] <- 0			}
    
    ### Interception storage calculations, no passive storage mixing volume ###
    
    Int_D[i] <- I[i] - D[i] + Ei[i] + Th[i] + INTp
    
    if (i==1){
      Int_CD[i] <- ((Int_D[i] * -50) + (D[i] * upCP_D[i]))/(Int_D[i] + D[i])
      Idl_H[i] <- da_H[i]
      fInt_CD[i] <- -50
    }
    else {
      Int_CD[i] <- ((Int_D[i] * Int_CD[i-1]) + (D[i] * upCP_D[i])) / (Int_D[i] + D[i]) 
      # compute the isotopic composition of the residual liquid
      Idl_H[i] <- (Int_CD[i] - dstar_H[i])*(1-x)^m_H[i] +dstar_H[i] #desiccating water body
      # compute vapor isotopic composition (Craig and Gordon 1965 formula, with notation by Gibson 2016)
      fInt_CD[i] <- ((Int_D[i] * Int_CD[i]) - (Ei[i] * Idl_H[i])) / (Int_D[i] - Ei[i]) 
      IdE_H[i] <- ((Idl_H[i] -epse_H[i])/alphae_H[i] -RH[i] *da_H[i] -epsk_H[i])/(1- RH[i] +10^-3*epsk_H[i]) #permil notation
    }
    
    ### Upper storage calculations with passive storage mixing volume ###
    
    upSTO_D[i] <- STO[i] - (P[i] - D[i]) - Th[i] + Es[i] + Tr_Upper[i] + Perc[i] + stoSp
    
    if (i==1){
      upCSTO_D[i] <- ((upSTO_D[i] * -52) + ((Th[i] * fInt_CD[i]) + (P[i] - D[i]) * upCP_D[i])) / (upSTO_D[i] + P[i] - D[i] + Th[i])
      Sdl_H[i] <- da_H[i]
      fupCSTO_D[i] <- -52
    }
    else {
      upCSTO_D[i] <- ((upSTO_D[i] * upCSTO_D[i-1]) + (Th[i] * fInt_CD[i]) + ((P[i] - D[i])*upCP_D[i]) ) / (upSTO_D[i] + P[i] - D[i] + Th[i])
      
      Sdl_H[i] <- (upCSTO_D[i] - dstar_H[i])*(1-x)^m_H[i] + dstar_H[i] #desiccating water body
      # compute vapor isotopic composition (Craig and Gordon 1965 formula, with notation by Gibson 2016)
      fupCSTO_D[i] <- ((upSTO_D[i] * upCSTO_D[i]) - (Es[i] * Sdl_H[i])) / (upSTO_D[i] - Es[i]) 
      
      SdE_H[i] <- ((Sdl_H[i] -epse_H[i])/alphae_H[i] -RH[i] *da_H[i] -epsk_H[i])/(1- RH[i] + 10^-3*epsk_H[i]) #permil notation
      
    }
    upCQ_D[i] <- fupCSTO_D[i]  # Percolation flux into GW
    
    # calculate concentration plus passive storage volume of GW storage
    
    gwSTO_D[i] <- GW[i] - Perc[i] - Pref_Flow[i] + GWflow[i] + Tr_Lower[i] + gwSp
    
    if (i==1){
      gwCSTO_D[i] <- ((gwSTO_D[i] * -66) + (Perc[i] * upCQ_D[i]))/(gwSTO_D[i] + Perc[i])
    }
    else {
      gwCSTO_D[i] <- ((gwSTO_D[i] * gwCSTO_D[i-1]) + (Perc[i] * upCQ_D[i]) + (Pref_Flow[i] * upCP_D[i])) /(gwSTO_D[i] + Perc[i] + Pref_Flow[i])
    }
    
    gwCQ_D[i] <- gwCSTO_D[i] # Percolation flux
    
    # calculate concentration plus passive storage volume of deep storage
    
    lowSTO_D[i] <- Sdeep[i] - GWflow[i] + Recharge[i] + Tr_Deep[i] + lowSp
    
    if (i==1){
      lowCSTO_D[i] <- ((lowSTO_D[i] * -70) + (GWflow[i] * gwCQ_D[i]))/(lowSTO_D[i] + GWflow[i])
    }
    else {
      lowCSTO_D[i] <- ((lowSTO_D[i] * lowCSTO_D[i-1]) + (GWflow[i] * gwCQ_D[i-1])) /(lowSTO_D[i-1] + GWflow[i])
    }
    
    lowCQ_D[i] <- lowCSTO_D[i] # Recharge flux
  }   
  # output simulated series
  output <- data.frame(SCF[367:nrun], D[367:nrun], Qs[367:nrun], Th[367:nrun], PN[367:nrun], I[367:nrun],
                       Ei[367:nrun], Es[367:nrun], Tr[367:nrun],Tr_Upper[367:nrun], Tr_Lower[367:nrun],
                       Tr_Deep[367:nrun], Perc[367:nrun], Pref_Flow[367:nrun],Recharge[367:nrun],
                       GWflow[367:nrun], STO[367:nrun], GW[367:nrun], Sdeep[367:nrun], 
                       Int_CD[367:nrun],Idl_H[367:nrun], fInt_CD[367:nrun],upCSTO_D[367:nrun],Sdl_H[367:nrun], fupCSTO_D[367:nrun],
                       gwCQ_D[367:nrun], lowCQ_D[367:nrun])
#  write.csv(data.frame(output), file = paste("SWMout_PARtest1.csv",sep=""))
   return (cbind(output))
}     ##### END OF FUNCTION BODY #####
