#FUNCTIONS FOR EQUILIBRIUM MODEL OF GALAXY FORMATION
# Contains various functions used by the equilibrium model


cosmictime <- function(z) {
#Calculate cosmic time given redshift
  
  tol = 10^{-5}
  it = 0
  aextemp = 1 /(1 +z)
  t = t0
  etalast=t+1
  
  while(abs(t-etalast)/etalast > tol)
  {
    f = sqrt(Lambda/totMass)*aextemp**1.5  + sqrt(aextemp*aextemp*aextemp*Lambda/totMass+1 ) - exp(1.5 *sqrt(Lambda)*H0*t)
    fprime =-1.5 *sqrt(Lambda)*H0*exp(1.5 *sqrt(Lambda)*H0*t)
    etalast = t
    t = t - f/fprime
  }
  etaold = t
  t
}

redshift <- function(t)
#Calculate redshift given time  
{
  tol = 1e-5
  
  it = 0
  eta = sqrt(1-totMass)*1.5 *H0*t
  aex = (sqrt(totMass/(1.-totMass))*sinh(eta))**(2 /3 )
  aex3 = aex*aex*aex
  hubble = H0*sqrt(totMass/aex/aex/aex+Lambda)
  aexhub = aex*hubble
  redshift = 1 /aex - 1 
  redshift
  
}

mhdot <- function(z,Mh){
#Calculate halo mass growth rate given halo mass and redshift  
  mhd = 0.47*10^{-9}*Mh*((Mh/10^{12})^0.15)*(((1+z)/3)^2.25)
  mhd
}

zeta <- function(z,Mh){
  Mphoto = 8-2.72 *(((1 +z)/9 )^0.2 )*log10((1+z)/9 )       # Okamoto+08, Fig 3
  zphoto = (1 + 0.587 *(Mh/(10 ^Mphoto))^(-2 ))^(-1.5 )        # Okamoto+08, eq 1
  
  Mquench = (0.96 +zetaparam1*z)*10^12 
  
  if (Mh>Mquench) {
    zquench = (Mh/Mquench)^(zetaslope)
  } else {
    zquench=1
  }
  
  zgrav = 0.47  * (((1 +z)/4 )^0.38 ) * (Mh/10^{12})^(-0.25 )       
  if(zgrav > 1) zgrav = 1 
  
  Mwind = 10 -zetaparam2*z
  zwind = 1 -exp(-sqrt(Mh/(10^Mwind)))     
  
  zeta = zphoto*zquench*zgrav*zwind
  zeta = zeta*0.5 	                           # fudge factor to match ampl of M*/Mh in Behroozi+12
  
  zeta
}

eta <- function(z,Mh){
  eta = (Mh/10 ^(etaslope1))^(etaslope2+etaparam*z)
  eta
}

sfr <- function(z,Mh){
#Calculate star formation rate as a function of halo mass and redshift.  
  sfr = mhdot(z,Mh)*fbaryon*zeta(z,Mh)/(1+eta(z,Mh))
  sfr
}

zeq <- function(fg,Mh,z){
  zeq = ((5 *fg*(1 +eta(z,Mh)))^(4 /3 )*(Mh/10^{12})^(-0.2))
  if(zeq < 2) zeq = 2
  zeq
}

fgas <- function(z,t,Ms,ssfr,Zg){
# Gas Fraction  
  tdep = t * 0.4 *(Ms/10^{10})^(-0.3 )	                  
  if( tdepslope > 0  ){ 
    Rdisk = 3 *(1 +z)^(-0.4 )*(Ms/5e10^{10})^(0.3333 )     
    Sigma = 183 *(Ms*ssfr/Rdisk/Rdisk)^(0.71 )             
    tdep = tdep * (Sigma+10 *0.02 /Zg)/Sigma
  }
  if( Zg < 0.0189  && tdepslope < 0  ) tdep = t * (Zg/0.0189 )**tdepslope * 0.4 *(Ms/1e10)**(-0.3 )
  
  fgas = 1 /(1 +1 /(tdep*ssfr))
  fgas
}

zgas <- function(z,Mh)  {
  zgas = yield*sfr(z,Mh)/(mhdot(z,Mh)*fbaryon*zeta(z,Mh))
  zgas
}
