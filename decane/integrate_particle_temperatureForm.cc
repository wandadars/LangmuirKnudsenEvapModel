#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
using namespace std;

struct NewParticle
{
  double vel[3] ;
  double temp[2] ;
  double pos[2] ;
  float mass[2] ;  
  double density ;
  double c ;
  double diameter ;
  double mdot ;
  double time ;


  //Debugging quantities
  double Term1;
  double Term2;
  void init(double dens, double cp, double T, double D)
  {
    vel[0] = vel[1] = vel[2] = 0 ;
    temp[0] = T ;
    temp[1] = T ;
    pos[0] = 0 ;
    pos[1] = 0 ;  
    density = dens ; 
    c = cp ;

    time=0.0 ;
    mdot=0.0 ;
    //using initial diameter, store initial mass of droplet
    mass[0] = mass[1] = density*(M_PI*pow(D,3)/6.0) ;


    //Debugging quantities
    double Term1=0;
    double Term2=0;

  }

  double get_diameter()
  {
    return max(pow(6.0*mass[0]/(M_PI*density),1./3.),1e-20) ;
  }
  double get_temperature()
  {
    return temp[0];
  }
  double get_density()
  {
    return density;
  }
  double get_mass_fraction()
  {
    return 1.0;
  }
  double get_c()
  {
    return c;
  }
  void printInfo()
  {
    cout << "Temperature: "<< get_temperature()<< "\t Diameter: "<<get_diameter()<<"\tVelocity "<<vel[0]<<"\t Position(x) "<<pos[0]<<endl;
  }

  void writeInfo()
  {
    ofstream ofs;
    ofs.open("data.txt", ios::out | ios::app);
    ofs<< get_temperature()<< "\t"<<get_diameter()<<"\t"<<vel[0]<<"\t"<<pos[0]<<"\t"<<mdot<<"\t"<<Term1<<"\t"<<Term2<<"\t"<<Term1+Term2<<"\t"<<mass[0]<<endl;
    ofs.close();

  }

};



//Function prototypes
double Langmuir_Knudsen_mdot(
		   double , // Particle diameter
		   double , // particle temperature
		   double , // saturation pressure
		   double , // Reynolds Number
           double , //gas viscosity
           double , //gas constant pressure heat capacity
           double , //gas thermal conductivity
           double , //carrier gas phase pressure
           double , //carrier gas phase gas constant
           double , // schmidt number
		   double , // vapor gas constant
           double  );

double integrateMassEnergy(
              NewParticle &,
		      float ,
              double ,
		      float ,
              float ,  //carrier fluid velocity
              double , //carrier density
              double , //carrier viscosity
              double , //heat capacity of carrier gas
              double , //thermal conductivity of carrier gas
              double , //pressure of gas phase
              double , //specific gas constant of carrier phase
              double , //Mass fraction of droplet vapor in carrier gas
              double , //Specific gas constant for vapor of droplet phase
              double )  ; //Schmidt number of carrier gas and droplet phase



void integrateMomentum(
            NewParticle & ,
		    float ,
            double , //carrier gas velocity
            double , // Volume allocated to the droplet
            double , //carrier gas viscosity
            double  //carrier gas density
		                )  ;



int main(void)
{

  //Numerics Section
  float dt = 1e-3; //timestep size(seconds)
  int IMAX = 10000; //Number of timesteps to take

  //Define background gas quantities(Properties evaludated at the wet-bulb temperature of 420K
  //Twb = 137*(T_B/373.15K)^(0.68)*log10(T_G) - 45K
  float T_g = 1000; //Temperature of gas (Kelvin)
  double P_g = 101325; //Pressure of gas (Pascals)
  double Y_g = 0.00; //Mass Fraction of vapor in carrier gas
  double mu_g = 2.410E-5; //Viscosity of gas (Pascals*seconds)
  double Sc_g = 2.7; //Schmidt number for gas diffusing into air http://webserver.dmt.upm.es/~isidoro/dat1/Mass%20diffusivity%20data.pdf
  double R_g = 286.9; //Specific gas constant for carrir(air), J/kgK
  float u_g = 0.58; //background fluid velocity(m/s) 
  double lambda_g = 0.03505; //carrier gas thermal conductivity
  double cp_g = 1015; //heat capacity of gas J/kgK

  double fvolpp = 8e-3; //Volume of container that droplet is in m^3

  //Define droplet quantities
  double D_p = 2.0e-3; //initial droplet diameter(meters)
  double T_p = 315; //Kelvin
  double r_p = 730; //kg/m^3 for water
  double cp_p = 2218; //heat capacity of the liquid particle J/kgK

  //Defined droplet vapor properties
  double R_v = 58.55; //gas constant for vapor phase

  //Compute the gas density
  double r_g = P_g/(R_g*T_g);
  cout<<"r_g: "<<r_g<<endl;

  //Initialize a particl
  NewParticle droplet;
  droplet.init(r_p,cp_p,T_p,D_p);


  double mdot=0.0;
  double m_air = r_g*fvolpp;
  double m_water = Y_g*m_air/(1-Y_g); //initial mass of water

  double time = 0.0;
  for(int i = 0;i<IMAX;i++)
  {

    droplet.printInfo();

    mdot = integrateMassEnergy(droplet,
                        dt,
                        fvolpp,
                        T_g,
                        u_g,
                        r_g,
                        mu_g,
                        cp_g,
                        lambda_g,
                        P_g,
                        R_g,
                        Y_g,
                        R_v,
                        Sc_g                   
                       );

    droplet.time = time;
    time = time + dt;

    //Update Y_g to account for droplet evaporation into the background space
    double dm = -mdot*dt;
    cout <<"Mdot: "<<mdot<<"\tAdded Water mass: "<<dm<<endl;
    //m_water= m_water + dm;
    Y_g = m_water/(m_water+m_air);
    cout <<"New Y_g = "<<Y_g<<endl;

    droplet.printInfo();

    /*
    integrateMomentum(droplet,
                      dt,
                      u_g,
                      fvolpp,
                      mu_g,
                      r_g
                    );

    */
    
    if(droplet.get_diameter() < 1e-8)
    {
        break;
    }



  droplet.writeInfo();
  
  }
  
  cout << "Program Finished" << endl;

  return 0;
}



  // Langmuir-Knudsen I non-equilibrium vaporization model


  inline double wilkeRuleProperty(double vapor_prop,double vapor_R, double gas_prop, double gas_R, double mol_fraction)
{
    //#Evaluate the Wilke rule for computing thermophysical properties
    const double theta= vapor_R/gas_R;

    omega_vg = (pow(pow(1.0+vapor_prop/gas_prop,1.0/2.0)*pow(theta,1.0/4.0),2.0))/(pow(8.0*(1.0+(1.0/theta)),1.0/2.0))
    omega_gv = (pow(pow(1.0+gas_prop/vapor_prop,1.0/2.0)*pow(1.0/theta,1.0/4.0),2.0)/((8.0*pow(1.0+(theta),1.0/2.0))

    mixture_prop= (ref_mol_frac*vapor_prop)/(ref_mol_frac + (1.0-ref_mol_frac)*omega_vg) + ((1.0-ref_mol_frac)*gas_prop)/(ref_mol_frac*omega_gv + (1.0-ref_mol_frac))

    return mixture_prop

}


  // Knudsen layer thickness
  inline double computeLK(double Tp, double R_v,double mu_g, double Sc_g, double P_g) 
  {
    double alpha_e = 1.0 ;
    return mu_g*sqrt(2.*M_PI*Tp*R_v)/(alpha_e*Sc_g*P_g) ;
  }

  inline double computeYsneq(
           double Xseq,    // equilibrium mole fraction
           double Lk,      // Knudsen Layer thickness
		   double D,       // particle diameter 
		   double theta2,  // vapor gas constant ratio
		   double Pr_g,    // Prandtl number
           double Re_b )   // Blowing Reynolds number
  {
    double beta = 0.5*Pr_g*Re_b ;
    double Xsneq = max(Xseq-2.0*beta*Lk/D,0.0) ;
    double Ysneq = Xsneq/(Xsneq+(1.0-Xsneq)*theta2) ;
    return Ysneq ;
  }

  double Langmuir_Knudsen_mdot(
		   double D, // Particle diameter
		   double T_p,   // particle temperature
		   double Psat, // saturation pressure
		   double Re, // Reynolds Number
           double mu_g, //gas viscosity
           double cp_g, //gas constant pressure heat capacity
           double lambda_g, //gas thermal conductivity
           double P_g, //carrier gas phase pressure
           double R_g, //carrier gas phase gas constant
           double Sc_g, // schmidt number
		   double R_v, // vapor gas constant
           double Yinf ) 
    {

      const double Pr_g = mu_g*cp_g/lambda_g ; // Gas Prandtl Number
      const double Sh = 2. + 0.552*sqrt(Re)*pow(Sc_g,1./3.) ;
      double Re_b = 0.0 ; //Blowing Reynolds number 
      double Re_b0 = Re_b ;
      double Xseq = min(Psat/P_g, 1.0) ; //Molar mass fraction
      double theta2 = R_v/R_g ;
      double Yseq = min(Xseq/max(Xseq+(1.-Xseq)*theta2,1e-30),1.0) ;
      double yMin = min(Yseq,Yinf) ;
      double yMax = max(Yseq,Yinf) ;

      // Iterate to converge on Re_b
      // This part could be optimized
      double Lk = computeLK(T_p,R_v,mu_g,Sc_g,P_g) ;
    
      { // Estimate starting Re_b
        double Ysneq = 0.5*(yMin+yMax) ;
        const double BMneq = (Ysneq-Yinf)/max(1.-Ysneq,1e-30) ;
        const double Hm = log(max(1.+BMneq,1e-40)) ;
        Re_b = Hm*Sh/(Sc_g) ;
        Re_b0 = Re_b ;
      }
      
      for(int i=0;i<40;++i) 
      {

        double Ysneq = computeYsneq(Xseq,Lk,D,theta2,Pr_g,Re_b) ;
        // bound Ysneq so that it lies between Yseq and Yinf
        Ysneq = max(yMin,min(yMax,Ysneq)) ;
        const double BMneq = (Ysneq-Yinf)/max(1.-Ysneq,1e-30) ;
        const double Hm = log(max(1.+BMneq,1e-40)) ;
        Re_b0 = Re_b ;
        Re_b = Hm*Sh/(Sc_g) ;
        const double factor = min(0.8, .5*D/Lk) ; // Damping factor

        if(i>=39)
        {
            cout<<"Mdot Calculation failed to converge"<<endl;
        }
        if(fabs(Re_b-Re_b0)<1e-8)
        {
          break ;
        }

        // Relax update to help convergence
        Re_b = factor*Re_b+(1.-factor)*Re_b0 ;

     
        //Chris debug
        const double beta = 0.5*Pr_g*Re_b ;
        if(i > -1) 
        {
          cout << i << "  Re_b=" << Re_b << "  Ysneq="
              << computeYsneq(Xseq,Lk,D,theta2,Pr_g,Re_b) 
              << "  Hm=" << Hm << "  BMneq=" << BMneq
              << "  Lk=" << Lk << "  Lk/D=" << Lk/D
              << "  factor=" << factor<<" beta="<< beta
              << endl;
         }
          
              //mdot = -Hm*D*Sh*mu_g/Sc_g ;
              //Re_b = -mdot/(D*mu_g*M_PI) ;
      }


      // back out mdot from blowing reunolds number
      const double  mdot = -Re_b*D*mu_g*M_PI ;
      return mdot ;
    }





  double integrateMassEnergy(
              NewParticle &p,
		      float dt,
              double fvolpp,
		      float fluid_T,
              float fluid_v, //carrier fluid velocity
              double rho_g, //carrier density
              double mu_g, //carrier viscosity
              double cp_g, //heat capacity of carrier gas
              double lambda_g, //thermal conductivity of carrier gas
              double P_g, //pressure of gas phase
              double R_g, //specific gas constant of carrier phase
              double Yinf, //Mass fraction of droplet vapor in carrier gas
              double R_v, //Specific gas constant for vapor of droplet phase
              double Sc_g)  { //Schmidt number of carrier gas and droplet phase


    // time advance coefficients for 2nd order implicit BDF
    const double beta = 2./3. ;
    const double alpha1 = -4./3. ;
    const double alpha2 = 1./3. ;
    const double dtbeta = dt*beta ;
    

    // time advance coefficients for 1st order implicit BDF
    //const double  beta = 1.0 ;
    //const double  alpha1 = -1.0 ;
    //const double  alpha2 = 0 ;
    //const double  dtbeta = dt*beta ;

    double D = p.get_diameter() ;
    const double T_p = p.get_temperature() ;

    ///const double Psat = PsatF.get_Psat(T_p) 
    double hvap = 288.7e3; //Chris heat of vaporization of decane at 420K J/kgK, moved from lower section because Psat needs it up here
    const double T_B = 447.7; //Boiling point of decane at 1atm in Kelvin
    const double P_atm = 101325; //Atmosphereic pressure in pascals
    const double Psat = P_atm*exp((hvap/R_v)*(1/T_B - 1/T_p));
    //const double Psat = 49.007e3; //hardcoded saturatoin pressure @ 420K 
    
    const float dv = fabs(fluid_v - p.vel[0]) ;
    const double Re = rho_g*D*dv/mu_g ;

    double mdot = Langmuir_Knudsen_mdot(
                    D, 
                    T_p, 
                    Psat, 
                    Re,
                    mu_g,
                    cp_g,
                    lambda_g,
                    P_g,
                    R_g, 
                    Sc_g, 
                    R_v,
                    Yinf) ;

    //Hardcoded mdots for debugging
    //mdot = -1.1879e-6*sqrt(-1.5125e-9*p.time + 1.21e-6); //Analytical mdot from miller's paper
    //mdot = -8.71137e-10;
    
    //Chris-store mdot
    p.mdot = mdot;

    //Hardcoded value for mcot for isolating energy equation from mass loss equation
    //mdot = -1.3e-9; //kg/s

    //cout << "Computed mdot: "<<mdot<<endl;
    //if(T_p < PsatF.get_Ttrip()) // Check if colder than triple point
    if(T_p < 273.16) // Chris Hardcoded triple point for water Check if colder than triple point
      mdot = max(mdot,0.0) ; // only allow condensation

    // Compute blowing Reynolds Number
    const double Re_b = max(-mdot/(D*mu_g*M_PI),0.0) ;
    const double Pr_g = mu_g*cp_g/lambda_g ;
    const double Nu = 2. + 0.552*sqrt(Re)*pow(Pr_g,1./3.) ;
    const double beta_b = 0.5*Pr_g*Re_b ;
    // Compute modification to heat transfer due to blowing
    const double f2 = ((beta_b > 1e-6)?beta_b/(exp(beta_b)-1.0):1.0) ;

    const double mf = p.get_mass_fraction() ;
    const double rho_p = p.get_density() ;

    double T1 = p.temp[0] ;
    float massp1 = dtbeta*(mdot) - alpha1*p.mass[0] - alpha2*p.mass[1] ;

    if(massp1 < 0.01*p.mass[0])
      massp1 = 0 ;
    
    // update mass of droplet
    p.mass[1] = p.mass[0] ;
    p.mass[0] = massp1 ;
    // Update droplet diameter
    D = max(p.get_diameter(),1e-10) ;

    //double hvap = 288.7e3; //Chris heat of vaporization of decane at 420K J/kgK

    bool find_hi = false ;
    double hlo = 0 ;
    double hhi = 1e100 ;
    //const double heatratio = p.mass[0]/max(mfluid*cp_g,1e-30) ;
    double oldfunc=1e300 ;
    const double fudge2 = 1.0;
    const double fixedsrc = alpha1*p.temp[0]+alpha2*p.temp[1] - dtbeta*hvap*mdot/(max(p.mass[0],1e-30f)*p.get_c()*fudge2);
    const double tau_p = rho_p*D*D/(18.*mu_g) ;
    cout << "rho_p: "<<rho_p<<"\tD: "<<D<<"\tmu_g: "<<mu_g<<endl;
    double fudge=1.0 ;
    const double Ecoef = dtbeta*fudge*f2*Nu*cp_g/(3.*Pr_g*tau_p*p.get_c()) ;
    
    const int ITMAX=100 ;
    int iter = 0 ;
      cout<<"f2: "<<f2<<"\tNu: "<<Nu<<"\tCp_g: "<<cp_g<<"\tPr: "<<Pr_g<<"\ttau_p: "<<tau_p<<"\tC_liquid: "<<p.get_c()<<"\tTerm1Coeff: "<<fudge*f2*Nu*cp_g/(3.*Pr_g*tau_p*p.get_c())<<"\tTerm2Coeff: "<< hvap*mdot/(max(p.mass[0],1e-30f)*p.get_c())<<"\tBeta: "<<beta_b<<"\tMass: "<<p.mass[0]<<"\tMdot: "<<p.mdot<<"\tXseq: "<<Psat/P_g<<endl;
    for(iter=0;iter<ITMAX;++iter) {

      cout<<"T: "<<T1<<endl;
      
      // Update fluid temperature to consider convection heat transfer
      // Temperature of fluid must range between fluid and particle initial
      // temperatures to be physical

      //double ptemp = energyp.get_temp(h1,mf) ;
      double ptemp = T1;
      //cout <<"ptemp: "<<ptemp<<"\th1: "<<h1<<endl;
      //      const double maxtemp = max(ptemp,double(fluid_T)) ;
      //      const double mintemp = min(ptemp,double(fluid_T)) ;
      //      double fluid_T_corr = min(maxtemp,
      //				max(mintemp,fluid_T - (h1-p.h[0])*heatratio)) ;

      double func = T1 + fixedsrc - Ecoef*(fluid_T-ptemp)  ;
      //cout<< "func: "<<func<<endl;

      if((iter > 2) && (fabs(func) < 1e-6*max(double(p.temp[0]),T1) || hhi-hlo < 1e-5))
        break ;
      if(func < 0) {
        hlo = T1 ;
      }
      if(func > 0) {
        hhi = T1 ;
        find_hi = true ;
      }
    
      T1 += -func;


      double factor=0.5;
      // Relax update to help convergence
      T1 = factor*T1+(1.0-factor)*ptemp ;


      // Try to bracket if possible by overshooting increasing predictions if
      // upper bound not found
      if(!find_hi && iter > 4 && func < 0 )
        T1 = 1.1*T1 ;
      // If out of bounds or not converging fast enough, switch to
      // bisection
//    if(find_hi && (h1 > max(hlo,hhi) || h1 < min(hlo,hhi)|| fabs(func/oldfunc) > .5)) {
      if(find_hi && (T1 > max(hlo,hhi) || T1 < min(hlo,hhi))) { // Removed last test to prevent FPE when oldfunc=0. JW 01/03/2017
        // Use bisection if outside of known root bounds
        T1 = 0.5*(hlo+hhi) ;
      }

      oldfunc = func ;

    }

    if(iter == ITMAX) {
      cerr << "particle temperature update failed to converge!" <<endl ;
    }
 
    cout<<"FluidT: "<<fluid_T<<"\tNew Droplet T: "<<T1<<endl;
    cout<<"f2: "<<f2<<"\tNu: "<<Nu<<"\tCp_g: "<<cp_g<<"\tPr: "<<Pr_g<<"\ttau_p: "<<tau_p<<"\tC_liquid: "<<p.get_c()<<"\tD: "<<p.get_diameter()<<endl;
    
    //These are the two terms in the energy equation
    double debug_term_1 = (fudge*f2*Nu*cp_g/(3.*Pr_g*tau_p*p.get_c()))*(fluid_T-T1);
    double debug_term_2 = hvap*mdot/(p.mass[0]*p.get_c()*fudge2);
    cout<<"Term 1: "<<debug_term_1<<"\tTerm 2: "<<debug_term_2<<endl;

    //Store debugging quantities
    p.Term1 = debug_term_1;
    p.Term2 = debug_term_2;


  
    p.temp[1] = p.temp[0] ; // Update T^n-1
    p.temp[0] = T1 ;   // Now T = T^n+1

    const double T_pn = p.get_temperature() ;
    //const double nPsat = PsatF.get_Psat(T_pn) ;
    const double nPsat = P_atm*exp((hvap/R_v)*(1/T_B - 1/T_pn));

    if(nPsat > P_g) { // boiling limit temperature to boiling point
      cout <<"Inside Boiling Temperature Limiting section"<<endl;
      
      // Find boiling temperature at this pressure
      double Tmax = T_pn ;
      //double Tmin = PsatF.get_Ttrip() ;
      double Tmin = 273.16; //Chris hardcoded triple point temperature of water ;
      double Tboil = 0.9*Tmax+0.1*Tmin ;
      int MAXITER=50 ;
      for(int i=0;i<MAXITER;++i) {
        double dPdT ;
        //double f = PsatF.get_Psat(Tboil,dPdT)-P_g ;
        double f = P_atm*exp((hvap/R_v)*(1/T_B - 1/Tboil)) - P_g ; //Chris hardcoded saturation pressure at boiling temperature
        if(fabs(f)<1e-4*P_g)
          break ;
        if(f > 0.0) {
          Tmax = Tboil ;
        }
        if(f < 0.0) {
          Tmin = Tboil ;
        }

        Tboil -= f/dPdT ;
        if(Tboil > Tmax || Tboil < Tmin)
          Tboil = 0.5*(Tmax+Tmin) ;
      }
      //      cerr << "boiling temperature switch activated, Tboil =" << Tboil
      //           << endl ;
      double mf = p.get_mass_fraction() ;
      //double hnew = energyp.get_h(Tboil,mf) ;
      //reset temperature to boiling temperature
      p.temp[0] = Tboil ;
    }


    return mdot;
  } 





  void integrateMomentum(
            NewParticle &p,
		    float dt,
            double fluid_v, //carrier gas velocity
            double fvolpp, // Volume allocated to the droplet
            double mu_g, //carrier gas viscosity
            double rfluid //carrier gas density
		                ) 
  {
    
    //integration constants
    const double beta = 2./3. ;
    const double alpha1 = -4./3. ;
    const double alpha2 = 1./3. ;
    const double dtbeta = dt*beta ;

    double vel1 = p.vel[0] ;
    double pos1 = dtbeta*vel1 - alpha1*p.pos[0]-alpha2*p.pos[1] ;
    const double rp = p.get_density() ;
    const double D = p.get_diameter() ;
    double mdot = (p.mass[0]-p.mass[1])/dt ;
    
    const double mfluid = rfluid*fvolpp+1e-30 ; // mass of fluid around particle
    double fixedsrc = -alpha1*p.vel[0]-alpha2*p.vel[1]  ;
    const double volp = M_PI*D*D*D/6.0 ;
    const double volpp = fvolpp ;
    // enhance drag function for large volume fraction
    const double alphav = min(2.0,volp/max(volpp,1e-30)) ;
    
    const double fp_vf = max(pow(8.*alphav,6.)-0.001,0.0) ;

    //Integration loop
    int i ;
    for(i=0;i<20;++i) {
      //Update fluid velocity based on delta particle momentum
      if(i>0)
        fluid_v = fluid_v - ((vel1-p.vel[0])*p.mass[0]/mfluid ) ;

      const float dv = fabs(fluid_v - vel1) ;
      const double Re = rfluid*D*dv/mu_g ;
      // blowing Reynolds number
      const double Re_b = fabs(mdot/(D*mu_g*M_PI)) ;
      const double a = 0.09+0.077*exp(-0.4*Re) ;
      const double b = 0.4+0.77*exp(-0.04*Re) ;
      const double denom = 1.0+a*pow(Re_b,b) ;

      const float fpblow = (1. + 0.0545*Re + 
                            0.1*sqrt(Re)*(1. - 0.03*Re))/denom + fp_vf ;
      // Clift-Gauvin drag function (Crowe, 1998)
      const float fpcg = 1.+0.15*pow(Re,0.687)+
        0.0175*Re/(1.+4.25e4*pow(Re+1e-20,-1.16)) + fp_vf ;
      // Choose drag function based on reynolds number.  For high reynolds
      // number use Clift Gauvin, otherwise use blowing reynolds number 
      const float fp = (Re< 100.0)?fpblow:fpcg ;
      const float taup = rp*pow(D,2)/(18.*mu_g*fp) ;
      const float vcoef = dtbeta/taup ;

      //      vel1 =  (vcoef*fluid_v + fixedsrc)/(1.+vcoef) ;
      float f = (vcoef*fluid_v + fixedsrc)/(1.+vcoef) - vel1 ;
      const float df = -vcoef*p.mass[0]/(mfluid*(1.+vcoef))-1.0 ;
      vel1 = vel1 - f/df ;
      pos1 = dtbeta*vel1 - alpha1*p.pos[0]-alpha2*p.pos[1] ;

      // If iterated at least 2 times, check for convergence
      if(i>1 && (fabs(f))/
         (fabs(df)*(.1+fabs(vel1)))< 1e-5) {
        break ;
      }

    }
    
    // Now advance the particle momentum in time
    p.vel[2] = p.vel[1] ;
    p.vel[1] = p.vel[0] ;
    p.vel[0] = vel1 ;
    p.pos[1] = p.pos[0] ;
    p.pos[0] = pos1 ;
  }



