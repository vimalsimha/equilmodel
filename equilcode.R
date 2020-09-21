source("equilfuncs.R")

#PARAMETERS

#Model Parameters (Can be fitted to data)
etaslope1 = 12
etaslope2 = -2/3
etaparam = 0
trec_param1 = 3
trec_param2 = -1/3
zetaslope = -1.083
zetaparam1 = 0.05
zetaparam2 = 0.25

#Basic (Cosmological and Physical) Parameters
fbaryon = 0.16428
yield = 0.0126 #Solar metal fraction from Asplund et al. 2009
tdepslope = 0
zend = -1
Hub = 70.4 # H_0 - Hubble Parameter
totMass = 0.27 # Omega_Matter
Lambda = 1-totMass #Assume flatness of Universe
t0=0


#EQUILIBRIUM CODE FOR GALAXY FORMATION

Mhalo = 9.225982  #Tuned for final Halo mass of ~1e14.

startflag = 1

H0 = sqrt(8.*pi/3.) 	                                     # work in tipsy units: Omega_crit=G=1

z0 = 9.0*10.0**((Mhalo-8.0)/2.72)-1	     # start roughly at Mphoto
Mhalo = 10**Mhalo
if( zend >= 0 ){
  tend = cosmictime(zend)
} else { 
tend = cosmictime(0)
}

unit_Time=H0*3.086e24/(Hub*1e5)/3.15576e7                 # present time = CosmicTime(0)*unit_Time in yr

etaold = 0.5                                             #this is for computing redshift to time; not massload
Mstar = 1e-10                                    	     # to avoid zeros to begin
t = cosmictime(z0) 
zeqflag = 0                                                # begin galaxy out of equilibrium
dt=1e7
k=1
j=1

Mrecyc = vector()
time = vector()
Deltat = vector()
Mdotrecyc = 0

while (t <= tend)
{
  redshiftz = redshift(t)
  dmh = dt*mhdot(redshiftz,Mhalo)                         # halo mass growth
  if( dmh > 0.01*Mhalo )  dt = dt*0.5                    # if too fast, reduce timestep
  if( dmh < 0.001*Mhalo ) dt = dt*2                     # if too slow, increase timestep
  
  #Compute Equilibrium Relations
  
  sfrate = (mhdot(redshiftz,Mhalo)*fbaryon*zeta(redshiftz,Mhalo)+Mdotrecyc)/(1+eta(redshiftz,Mhalo))
  
  Zg     = yield*sfrate/(mhdot(redshiftz,Mhalo)*fbaryon*zeta(redshiftz,Mhalo))
  
  fg     = fgas(redshiftz,t*unit_Time,Mstar,sfrate/Mstar,Zg)
  
  # Recycling Calculation
  
  trec = (trec_param1*10 **9 *(Mhalo/10 **12 )**(trec_param2))/unit_Time  # Now the trec and t have the same units
  
  Rrecyc = 0.046 *log(t*unit_Time/(2.76 *10 **5 ) + 1 ) #from Leitner & Kravtsov (2011) arXiv:1011.1252
  
  DeltaMstar = sfrate*dt-Rrecyc*dt*sfrate
  Mrecyc[k] = eta(redshiftz,Mhalo)*DeltaMstar
  Mrecyc[1] = 0 
  time[k] = t+trec
  Deltat[k] = dt
  
  Mdotrecyc = 0 
  Mrecyc_avg = 0 
  kount=1
  
  
  
  for (j in 1:k)
  {
    if (k>1 && t>=time[j] && time[j]>(t-dt/unit_Time))
    {
      Mrecyc_avg = Mrecyc_avg+Mrecyc[j]
      kount=kount+1
    }
  }
  
  Mrecyc_avg = Mrecyc_avg/kount
  Mdotrecyc = Mrecyc_avg/dt

  #Grow Halo
  Mhalo = Mhalo+dmh
  
  # check if equilibrium reached
  if( zeqflag == 0 && redshiftz < zeq(fg,Mhalo,redshiftz) ) zeqflag = 1
  # grow stellar mass if in equilibrium, else no growth
  if( zeqflag == 1 ) Mstar = Mstar + dt*sfrate-Rrecyc*dt*sfrate
  
  if(k==1)
  {
    outframe <- data.frame("time"=t*unit_Time/1.e9,"redshift"=redshiftz,"SFR"=log10(sfrate),"Metals"=log10(Zg/yield),"GasFrac"=fg,"M_halo"=Mhalo,"M_star"=Mstar)
  }else{
    outframe <- rbind(outframe,list(t*unit_Time/1.e9,redshiftz,log10(sfrate),log10(Zg/yield),fg,Mhalo,Mstar))  
  }
  
  t = t + dt/unit_Time
  k = k+1
  
}

write.table(outframe,file="outputs.txt")

  
  
  




