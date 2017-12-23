#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "Helium.H"
/****************************/
/****************************/
void Helium::setPlasmaParameters(plasmaParameters& a_ppars, Vector<Real>& a_minData, Vector<Real>& a_maxData,  string& a_bgSpec, Vector<string>& a_charSpec)
{
  
  Real nreac = 2;
  Real nspec = 3;
  Real ncomp = nspec + 1;
  // Reaction Indexes 
  a_ppars.reactions.resize(nreac);
  a_ppars.actReac.resize(nreac,1.0);   
  for (int rnr = 0; rnr < a_ppars.reactions.size(); rnr++)
    {
      a_ppars.reactions[rnr]= rnr;
    }
  Vector<int>   r = a_ppars.reactions;
   
  a_ppars.reacRateCoef.resize(nreac, 12);
  a_ppars.reactant.resize(nreac, vector<Real>(ncomp));  

  a_bgSpec = string("He");
  string arr[] = { "He", "He+", "e-", "ee" };
  a_charSpec = vector<string>(arr, end(arr));

  //floor and ceiling
  a_minData.resize(a_charSpec.size(),1e-30);
  a_maxData.resize(a_charSpec.size(),1e30);
  Real ne = 1e0;
  a_ppars.initNe = ne;
  a_minData[2]= ne/1e5;
  a_maxData[2]= ne*1e4;

  // Species Indexes
  a_ppars.species.resize(nspec); 
  for (int snr = 0; snr < a_ppars.species.size(); snr++)
    {
      a_ppars.species[snr]= snr;
    }
	
  Vector<int>   s = a_ppars.species;
	
  // first vector component is type of fit 
  // next parameters are for fit:
  // log f(x)=sumj pj/T^j, j=0...9 for type -1 (2nd param is threshold Temp)
  // f(x)=a*T^b for type 0 
  // f(x)=a*T^b*exp(c/T) for type 1 
  
  // Electron Temperature Dependent Reaction Rates
  // R0 Ionization He + e -> He+ + 2e

  a_ppars.reacRateCoef[r[0]][0] = 0;
  a_ppars.reacRateCoef[r[0]][1] = 1.0; // rion=alp*mu*E calculated in reactionRates as f(Te,E)
  a_ppars.reacRateCoef[r[0]][2] = 0;
  a_ppars.reacRateCoef[r[0]][3] = 0;
  a_ppars.reacRateCoef[r[0]][4] = 0;
  a_ppars.reacRateCoef[r[0]][5] = 0;
  a_ppars.reacRateCoef[r[0]][6] = 0;
  a_ppars.reacRateCoef[r[0]][7] = 0;
  a_ppars.reacRateCoef[r[0]][8] = 0;
  a_ppars.reacRateCoef[r[0]][9] = 0;
  a_ppars.reacRateCoef[r[0]][10] = 0;
  a_ppars.reacRateCoef[r[0]][11] = 0;

  // Source Terms
  a_ppars.reactant[r[0]][s[0]] = 0; // here rate = rion * ne => ng coeff = 0
  a_ppars.reactant[r[0]][s[1]] = 0;
  a_ppars.reactant[r[0]][s[2]] = 1; // e- 
  a_ppars.reactant[r[0]][nspec] = 0; //M 

  // R1 Recombination He+ + e -> He
 
  a_ppars.reacRateCoef[r[1]][0] = 0;
  a_ppars.reacRateCoef[r[1]][1] = 1.0; // recb calculated in reactionRates as f(Ng)
  a_ppars.reacRateCoef[r[1]][2] = 0;
  a_ppars.reacRateCoef[r[1]][3] = 0;
  a_ppars.reacRateCoef[r[1]][4] = 0;
  a_ppars.reacRateCoef[r[1]][5] = 0;
  a_ppars.reacRateCoef[r[1]][6] = 0;
  a_ppars.reacRateCoef[r[1]][7] = 0;
  a_ppars.reacRateCoef[r[1]][8] = 0;
  a_ppars.reacRateCoef[r[1]][9] = 0;
  a_ppars.reacRateCoef[r[1]][10] = 0;
  a_ppars.reacRateCoef[r[1]][11] = 0;

  // Source Terms
  a_ppars.reactant[r[1]][s[0]] = 0;
  a_ppars.reactant[r[1]][s[1]] = 1; // He+
  a_ppars.reactant[r[1]][s[2]] = 1; // e-
  a_ppars.reactant[r[1]][nspec] = 0; //M 
	
  // Mobility & Diffusion Coefficients
  // Mobility of ions and electrons = first 10 vector components
  // mu=f(x)/Ng, with f(x) = exp(p9*x^9 + p8*x^8 + p7*x^7 + p6*x^6 + p5*x^5 + p4*x^4 + p3*x^3 + p2*x^2 + p1*x + p0) & x=log(E/Ngface)
  // For diffusion of all species: g(x)=fact*f(x)^alp*Tvar^eta/nk^beta
  // where:
  // Tvar is Te for electrons ([0]=1) and Tg otherwise ([0]=0)
  // fact = kb/q/abs(Zk) for electrons/ions, sqrt(Ru/(2*mk))/sigmak ([11]) assuming gk=sqrt(2*Ru*Tk/mk) 
  // eta = 1 for electrons/ions, 0.5 for neutrals  (hardwired per [12])
  // alp = 1 for electrons/ions, 0 for neutrals ([13])
  // beta = 1 for neutrals, 0 otherwise ([14])

  // unified atomic mass unit [g/(NA mol)]
  Real ua = 1.660538782e-27; // [kg]
  Real Navo = 6.02214179e23;
  Real kgpmol = 1e-3;// /Navo;
  a_ppars.NA = Navo;
  a_ppars.mass.resize(nspec);
  // masses in [kg] 
  a_ppars.mass[s[0]]= 4.0026*kgpmol; //He
  a_ppars.mass[s[1]]= 4.0026*kgpmol; //He+
  a_ppars.mass[s[2]]= 9.10938215e-31 * a_ppars.NA;//e-
  
  a_ppars.kb = 1.3806504e-23;
  a_ppars.q = 1.602176487e-19; // [C]
  Real mbg = a_ppars.mass[s[0]];

  // Neutral momentum transfer x-section [m^2]
  Real sigbg = 1E-18; // [m^2] MDCheck [A. V. Phelps average vs. energy range]
  a_ppars.porder = 9;	  
  a_ppars.transpCoef.resize(nspec, Vector<Real>(15,0));	
  a_ppars.netReac.resize(nspec, Vector<Real>(nreac,0));	
		
  // Neutral	
  // He
  a_ppars.transpCoef[s[0]][0] = 0;
  a_ppars.transpCoef[s[0]][1] = 0;	
  a_ppars.transpCoef[s[0]][2] = 0;
  a_ppars.transpCoef[s[0]][3] = 0;
  a_ppars.transpCoef[s[0]][4] = 0;
  a_ppars.transpCoef[s[0]][5] = 0;
  a_ppars.transpCoef[s[0]][6] = 0;
  a_ppars.transpCoef[s[0]][7] = 0;
  a_ppars.transpCoef[s[0]][8] = 0;
  a_ppars.transpCoef[s[0]][9] = 0;
  a_ppars.transpCoef[s[0]][10] = 0;
  a_ppars.transpCoef[s[0]][11] = sqrt(a_ppars.kb*a_ppars.NA/mbg)/sigbg;   
  a_ppars.transpCoef[s[0]][12] = 0.5;
  a_ppars.transpCoef[s[0]][13] = 0;
  a_ppars.transpCoef[s[0]][14] = 1;

  // Ion	
  // He+
  // 4th degree polynom for ln(mob*ng)=f(ln(te)) 
  a_ppars.transpCoef[s[1]][0] = 0;
  a_ppars.transpCoef[s[1]][1] =  1.12731289E+04;	
  a_ppars.transpCoef[s[1]][2] =  1.02784380E+03;
  a_ppars.transpCoef[s[1]][3] =  3.52079699E+01;
  a_ppars.transpCoef[s[1]][4] =  5.34690992E-01;
  a_ppars.transpCoef[s[1]][5] =  3.03832413E-03;
  a_ppars.transpCoef[s[1]][6] =  0;
  a_ppars.transpCoef[s[1]][7] =  0;
  a_ppars.transpCoef[s[1]][8] =  0;
  a_ppars.transpCoef[s[1]][9] =  0;
  a_ppars.transpCoef[s[1]][10] = 0;
  a_ppars.transpCoef[s[1]][11] = a_ppars.kb/a_ppars.q; 
  a_ppars.transpCoef[s[1]][12] = 1;
  a_ppars.transpCoef[s[1]][13] = 1;
  a_ppars.transpCoef[s[1]][14] = 0;

  // Electrons
  // 8th degree polynom for ln(mob*ng)=f(ln(te))
  a_ppars.transpCoef[s[2]][0] = 1;
  a_ppars.transpCoef[s[2]][1] = 1.17247204E+02;
  a_ppars.transpCoef[s[2]][2] = -7.25346504E+01;
  a_ppars.transpCoef[s[2]][3] = 3.91135676E+01;
  a_ppars.transpCoef[s[2]][4] = -1.16742629E+01;
  a_ppars.transpCoef[s[2]][5] = 2.09224759E+00;
  a_ppars.transpCoef[s[2]][6] = -2.30502137E-01;
  a_ppars.transpCoef[s[2]][7] = 1.52502606E-02;
  a_ppars.transpCoef[s[2]][8] = -5.54542542E-04;
  a_ppars.transpCoef[s[2]][9] = 8.50025754E-06;
  a_ppars.transpCoef[s[2]][10] = 0;
  a_ppars.transpCoef[s[2]][11] = a_ppars.kb/a_ppars.q; 
  a_ppars.transpCoef[s[2]][12] = 0; 
  a_ppars.transpCoef[s[2]][13] = 1;
  a_ppars.transpCoef[s[2]][14] = 0;

  // net molecular production
  //(R0) e+He->2e+He+
  a_ppars.netReac[s[0]][r[0]] = -1;//He 
  a_ppars.netReac[s[1]][r[0]] = 1; //He+ 
  a_ppars.netReac[s[2]][r[0]] = 1; //e

  //(R1) e+He+->He
  a_ppars.netReac[s[0]][r[1]] = 1; //He 
  a_ppars.netReac[s[1]][r[1]] = -1; //He+
  a_ppars.netReac[s[2]][r[1]] = -1; //e

  // Wall Electron Secondary Emission Coefficient and Energy
  // Tungsten electrode
	
  a_ppars.secEmElec.resize(nspec,0e0);
  a_ppars.secEmDiel.resize(nspec,0e0);
  a_ppars.secEmEnergy.resize(nspec,0e0);
	
  a_ppars.secEmElec[s[0]]= 0;
  a_ppars.secEmElec[s[1]]= 0.26;// on W electrode [Braithwaite Plasma Sources Sci. Techn. 2000];
  a_ppars.secEmElec[s[2]]= 0;

  a_ppars.secEmDiel[s[0]]= 0;
  a_ppars.secEmDiel[s[1]]= 0; // MDCheck
  a_ppars.secEmDiel[s[2]]= 0;
		
  a_ppars.secEmEnergy[s[0]]= 0;
  a_ppars.secEmEnergy[s[1]]= 1.9; //[eV]
  a_ppars.secEmEnergy[s[2]]= 0;

  // Electron Energy inelastic collision losses

  a_ppars.eInelEn.resize(nreac);
  a_ppars.gInelEn.resize(nreac);

  a_ppars.eInelEn[r[0]] = 15.4;
  a_ppars.eInelEn[r[1]] = 0; //MDCheck


  // Gas Energy inelastic collision losses
  a_ppars.gInelEn[r[0]] = 0; 
  a_ppars.gInelEn[r[1]] = 0; //MDcheck should be non-zero, in this example OK because approx of Tg=300K

  // Charge numbers
  a_ppars.chargeNr.resize(nspec,0);
  a_ppars.chargeNr[s[0]] = 0;
  a_ppars.chargeNr[s[1]] = 1;
  a_ppars.chargeNr[s[2]] = -1;

  // fraction of Joule heating transfer to the gas
  a_ppars.etaT = 0.25; //MDCorr check literature [47] 

  // 7th degree polynom log10(sig) = f(log10(te))
  a_ppars.sigebg.resize(8,0);
  a_ppars.sigebg[0] = 6.350614E+00; 
  a_ppars.sigebg[1] = -5.955160E+01; 
  a_ppars.sigebg[2] = 5.467754E+01; 
  a_ppars.sigebg[3] = -2.615789E+01;
  a_ppars.sigebg[4] = 7.106106E+00;
  a_ppars.sigebg[5] = -1.101330E+00; 
  a_ppars.sigebg[6] = 9.048729E-02;
  a_ppars.sigebg[7] = -3.055798E-03; 

  // Electron elastic collision x-section 
  // triple Gaussian fit
  a_ppars.sigeo2.resize(9,0);
  a_ppars.sigeo2[0] = -2.605857E+02; 
  a_ppars.sigeo2[1] = -9.824257E+00; 
  a_ppars.sigeo2[2] = 3.437898E+01; 
  a_ppars.sigeo2[3] = 2.377912E-01;
  a_ppars.sigeo2[4] = -4.073995E+01;
  a_ppars.sigeo2[5] = 8.932246E-01; 
  a_ppars.sigeo2[6] = 1.200053E+02;
  a_ppars.sigeo2[7] = -2.861041E+01; 
  a_ppars.sigeo2[8] = 1.696869E+01;

  a_ppars.stickFac2bg.resize(nspec,0);
  a_ppars.stickFac2bg[s[0]] = 0;//He
  a_ppars.stickFac2bg[s[1]] = 1;//He+
  a_ppars.stickFac2bg[s[2]] = 1;//e-

  a_ppars.stickFac2o2.resize(nspec,0);
  a_ppars.stickFac2o2[s[0]] = 0;
  a_ppars.stickFac2o2[s[1]] = 0;
  a_ppars.stickFac2o2[s[2]] = 0;
}
