import math

class Particle(object):
    def __init__(self, dens, cp, temp, diameter):
        self.vel = 3*[0] 
        self.temp = 2*[temp] 
        self.pos = 2*[0] 
        self.density = dens 
        self.c = cp
        self.diameter = diameter         
        self.mass = 2*[self.density * math.pi * self.diameter ** 3 / 6]  
        self.mdot = [] 
        self.time = [] 

        #Debugging quantities
        self.Term1 = 0.0
        self.Term2 = 0.0

    def get_c(self):
        return self.c

    def get_density(self):
        return self.density

    def get_temperature(self):
        return self.temp[0]

    def get_diameter(self):
        self.diameter = max( (6.0*self.mass[0] / (math.pi*self.density) ) ** (1./3.), 1e-20)
        return self.diameter

    def get_mass_fraction(self):
        return 1.0

    def print_info(self):
        print('Temperature: {0:<8.4f}\tDiameter: {1:<8.4e}\tVelocity: {2:<8.4e}\tPosition(x): {3:<8.4e}\n'.format(self.temp[0], self.diameter, self.vel[0], self.pos[0]))

    def write_info(self):
        f = open('data.txt', 'w+') 
        format_string = '{0:<8.4f}\t{1:<8.4f}\t{2:<8.4f}\t{3:<8.4f}\t{4:<8.4f}\t{5:<8.4e}\t{6:<8.4f}\t{7:<8.4f}\t{8:<8.4f}\t{9:<8.4e}\n'
        f.write(format_string.format(self.time, self.temp[0], self.diameter, self.vel[0], self.pos[0], self.mdot, self.Term1, self.Term2, self.Term1 + self.Term2, self.mass[0]))
        f.close()


class saturation_data(object):
    def __init__(self, file_name):
        self.data_array = []
        self.data_file_name = file_name

    def count_data_file_lines(self):
        number_of_lines = 0
        with open(self.data_file_namefname) as f:
            for i, l in enumerate(f):
                number_of_lines += 1
        return number_of_lines;

    def count_data_file_columns(self):
        print('Filename to open: ', self.data_file_name)
        f = open(self.data_file_name, 'r')
        return len(f.readline().rstrip().split()) 

    def read_table(self, file_name):
        num_lines = self.count_data_file_lines()
        num_columns = self.count_data_file_columns()
        print('num Lines', num_lines, ' Num Columns: ', num_columns) 
        try:
            with open(self.data_file_name) as f:
                print(f.readline())
        except IOError: 
            print('Unable to open file')

        #read lines into data_array


class UserData(object):
    def __init__(self):
        #Numerics Section
        self.dt = 1.0e-3 #timestep size(seconds)
        self.num_timesteps = 10000 #Number of timesteps to take

        #Background gas quantities(Properties evaludated at the wet-bulb temperature of 420K
        #Twb = 137*(T_B/373.15K)^(0.68)*math.log10(T_G) - 45K
        self.T_g = 1000  #Temperature of gas (Kelvin)
        self.P_g = 101325  #Pressure of gas (Pascals)
        self.Y_g = 0.005  #Mass Fraction of vapor in carrier gas
        self.mu_g = 2.410e-5  #Viscosity of gas (Pascals*seconds)
        self.Sc_g = 2.7  #Schmidt number for gas diffusing into air http://webserver.dmt.upm.es/~isidoro/dat1/Mass%20diffusivity%20data.pdf
        self.R_g = 286.9  #Specific gas constant for carrir(air), J/kgK
        self.u_g = 0.6  #Background fluid velocity(m/s) 
        self.lambda_g = 0.043958  #Carrier gas thermal conductivity W/mK
        self.cp_g = 1045  #Heat capacity of gas J/kgK

        self.fvolpp = 8e-3 #Volume of container that droplet is in m^3

        #Define droplet quantities
        self.D_p = 2.0e-3 #Initial droplet diameter(meters)
        self.T_p = 315  #Kelvin
        self.r_p = 681  #kg/m^3 
        self.cp_p = 2413  #Heat capacity of the liquid particle J/kgK

        #Defined droplet vapor properties
        self.R_v = 58.55  #gas constant for vapor phase

    def read_user_data(self):
        pass


def main():
  input_data = UserData()

  #Numerics Section
  dt = input_data.dt #timestep size(seconds)
  num_timesteps = input_data.num_timesteps #Number of timesteps to take

  #Define background gas quantities(Properties evaludated at the wet-bulb temperature of 420K
  #Twb = 137*(T_B/373.15K)^(0.68)*math.log10(T_G) - 45K
  T_g = input_data.T_g #Temperature of gas (Kelvin)
  P_g = input_data.P_g #Pressure of gas (Pascals)
  Y_g = input_data.Y_g #Mass Fraction of vapor in carrier gas
  mu_g = input_data.mu_g #Viscosity of gas (Pascals*seconds)
  Sc_g = input_data.Sc_g #Schmidt number for gas diffusing into air http://webserver.dmt.upm.es/~isidoro/dat1/Mass%20diffusivity%20data.pdf
  R_g = input_data.R_g #Specific gas constant for carrir(air), J/kgK
  u_g = input_data.u_g #background fluid velocity(m/s) 
  lambda_g = input_data.lambda_g #carrier gas thermal conductivity
  cp_g = input_data.cp_g #heat capacity of gas J/kgK

  fvolpp = input_data.fvolpp #Volume of container that droplet is in m^3

  #Define droplet quantities
  D_p = input_data.D_p #initial droplet diameter(meters)
  T_p = input_data.T_p #Kelvin
  r_p = input_data.r_p #kg/m^3 
  cp_p = input_data.cp_p #heat capacity of the liquid particle J/kgK

  #Defined droplet vapor properties
  R_v = input_data.R_v #gas constant for vapor phase

  #Compute the gas density
  r_g = P_g/(R_g*T_g)
  print('Gas density: {0:<10.6f}'.format(r_g))

  #Initialize a particle
  droplet = Particle(r_p,cp_p,T_p,D_p) 

  mdot=0.0
  m_air = r_g*fvolpp
  m_water = Y_g*m_air / (1 - Y_g)  #initial mass of water

  time = 0.0
  for i in range(num_timesteps):
    print('Timestep: {0:<5d}'.format(i+1))
    droplet.print_info()
    mdot = integrateMassEnergy(droplet, dt, fvolpp, T_g, u_g, r_g, mu_g, cp_g, lambda_g, P_g, R_g, Y_g, R_v, Sc_g)

    droplet.time = time
    time = time + dt

    #Update Y_g to account for droplet evaporation into the background space
    dm = -mdot*dt
    print('Mdot: {0:<8.4e}\tAdded Water Mass: {1:<8.4e}'.format(mdot, dm))
    #m_water= m_water + dm
    Y_g = m_water / (m_water + m_air)
    print('New Y_g = {0:<8.4f}'.format(Y_g))

    droplet.print_info()
    """
    integrateMomentum(droplet,
                      dt,
                      u_g,
                      fvolpp,
                      mu_g,
                      r_g
                    )

    """
    
    if droplet.get_diameter() < 1e-8:
        break

    droplet.write_info()
    
  print('Program Finished')



def  wilke_rule_property(vapor_prop, vapor_R, gas_prop, gas_R, mol_frac):
    #Evaluate the Wilke rule for computing thermophysical properties
    theta = vapor_R / gas_R

    numerator = ( (1.0 + vapor_prop / gas_prop) ** (1.0/2.0) * theta ** (1.0/4.0) ) ** 2.0
    denominator = (8.0 * (1.0 + 1.0 / theta)) ** (1.0/2.0)
    omega_vg = numerator / denominator 
    
    numerator = ( (1.0 + gas_prop/vapor_prop) ** (1.0/2.0) * (1.0/theta) ** (1.0/4.0) ) ** 2.0 
    denominator = (8.0 * (1.0 + theta) ) ** (1.0/2.0)
    omega_gv = numerator / denominator 

    mixture_prop = (mol_frac * vapor_prop) / (mol_frac + (1.0 - mol_frac) * omega_vg) + ((1.0 - mol_frac) * gas_prop) / (mol_frac * omega_gv + (1.0 - mol_frac))

    return mixture_prop

def estimate_re_b(y_min, y_max, Y_inf, Sh, Sc_g):
    # Estimate starting Re_b
    Y_sneq = 0.5 * (y_min + y_max) 
    BMneq = (Y_sneq - Y_inf) / max(1.0 - Y_sneq, 1e-30) 
    Hm = math.log(max(1.0 + BMneq, 1e-40)) 
    Re_b = Hm * Sh / Sc_g
    return Re_b

  #Knudsen layer thickness
def computeLK( T_p, R_v, mu_g, sc_g, p_g):
    alpha_e = 1.0 
    return mu_g * math.sqrt(2.0 *math.pi * T_p * R_v) / (alpha_e * sc_g * p_g) 

def computeYsneq( Xseq, Lk, D, theta2, Pr_g, Re_b):   
    """
    Xseq    // equilibrium mole fraction
    Lk      // Knudsen Layer thickness
    D       // particle diameter
    theta2  // vapor gas constant ratio
    Pr_g    // Prandtl number
    Re_b    // Blowing Reynolds number
    """
    beta = 0.5 * Pr_g * Re_b 
    Xsneq = max(Xseq - 2.0 * beta * Lk / D , 0.0) 
    Ysneq = Xsneq /(Xsneq + (1.0 - Xsneq) * theta2) 
    return Ysneq 

  #Langmuir-Knudsen I non-equilibrium vaporization model
def Langmuir_Knudsen_mdot(D, T_p, Psat, Re, mu_g, cp_g, lambda_g, P_g, R_g, Sc_g,  R_v, Yinf): 
      """
      D       # Particle diameter
      T_p     # particle temperature
      Psat    # saturation pressure
      Re      # Reynolds Number
      mu_g    #gas viscosity
      cp_g    #gas constant pressure heat capacity
      lambda_g  #gas thermal conductivity
      P_g     #carrier gas phase pressure
      R_g     #carrier gas phase gas constant
      Sc_g    # schmidt number
      R_v     # vapor gas constant
      Yinf    #Vapor mass fraction in background ambient
      """
      Pr_g = mu_g * cp_g / lambda_g # Gas Prandtl Number
      Sh = 2.0 + 0.552 * math.sqrt(Re) * Sc_g ** (1.0/3.0) 
      Re_b = 0.0  #Blowing Reynolds number 
      Re_b0 = Re_b 
      Xseq = min(Psat / P_g, 1.0)  #Molar mass fraction
      theta2 = R_v / R_g ;
      Yseq = min(Xseq /max(Xseq + (1.0 - Xseq) * theta2, 1e-30), 1.0) 
      yMin = min(Yseq, Yinf) 
      yMax = max(Yseq, Yinf) 

      # Iterate to converge on Re_b
      # This part could be optimized
      Lk = computeLK(T_p,R_v,mu_g,Sc_g,P_g) 
      Re_b0 = estimate_re_b(yMin, yMax, Yinf, Sh, Sc_g)
      
      max_solver_iterations = 40
      for i in range(max_solver_iterations): 
        Ysneq = computeYsneq(Xseq,Lk,D,theta2,Pr_g,Re_b) 
        #Bound Ysneq so that it lies between Yseq and Yinf
        Ysneq = max(yMin, min(yMax, Ysneq)) 
        BMneq = (Ysneq - Yinf) / max(1.0 - Ysneq, 1e-30) 
        Hm = math.log(max(1.0 + BMneq, 1e-40)) 
        Re_b0 = Re_b 
        Re_b = Hm * Sh / Sc_g 
        factor = min(0.8, 0.5 * D / Lk)  #Damping factor

        if i >= 39:
            print('Mdot Calculation failed to converge')

        if abs(Re_b - Re_b0) < 1.0e-6:
          break 

        #Relax update to help convergence
        Re_b = factor * Re_b + (1.0 - factor) * Re_b0 
     
        #Chris debug
        beta = 0.5 * Pr_g * Re_b ;
        if i > -1: 
            format_string = 'i= {0:<4d}   Re_b= {1:<8.4f}   Ysneq= {2:<8.4f}   Hm= {3:<8.4f}   ' \
                            'BMneq= {4:<8.4f}   Lk= {5:<8.4e}   Lk/D= {6:<8.4e}   factor= {7:<8.4e}   beta= {8:<8.4f}'
            print(format_string.format(i, Re_b, computeYsneq(Xseq,Lk,D,theta2,Pr_g,Re_b), Hm, BMneq, Lk, Lk/D, factor, beta))
          
      # Back out mdot from blowing reunolds number
      mdot = -Re_b * D * mu_g * math.pi 
      return mdot 


def integrateMassEnergy( p, dt, fvolpp, fluid_T, fluid_v, rho_g, mu_g, cp_g, lambda_g, P_g, R_g, Yinf, R_v,  Sc_g):  
    """
    p  #Particle
    dt  #Timestep
    fvolpp  #Volume around particle
    fluid_T  #Carrier fluid temperature
    fluid_v  #carrier fluid velocity
    rho_g  #carrier density
    mu_g  #carrier viscosity
    cp_g  #heat capacity of carrier gas
    lambda_g  #thermal conductivity of carrier gas
    P_g  #pressure of gas phase
    R_g  #specific gas constant of carrier phase
    Yinf  #Mass fraction of droplet vapor in carrier gas
    R_v  #Specific gas constant for vapor of droplet phase
    Sc_g)  #Schmidt number of carrier gas and droplet phase
    """

    #Time advance coefficients for 2nd order implicit BDF
    beta = 2.0/3.0 
    alpha1 = -4.0/3.0 
    alpha2 = 1.0/3.0 
    dtbeta = dt * beta 

    #time advance coefficients for 1st order implicit BDF
    # beta = 1.0 
    # alpha1 = -1.0 
    # alpha2 = 0  

    dtbeta = dt * beta 

    D = p.get_diameter() 
    T_p = p.get_temperature() 

    # Psat = PsatF.get_Psat(T_p) 
    hvap = 299.7e3  #Chris heat of vaporization of decane at 420K J/kgK, moved from lower section because Psat needs it up here
    T_B = 447.27  #Boiling point of decane at 1atm in Kelvin
    P_atm = 101325  #Atmosphereic pressure in pascals
    Psat = P_atm * math.exp((hvap / R_v) * (1.0 / T_B - 1.0 / T_p))  #Clausius Clapyron
    # Psat = 49.007e3 #hardcoded saturation pressure @ 420K 
    
    dv = abs(fluid_v - p.vel[0]) 
    Re = rho_g * D * dv / mu_g 

    mdot = Langmuir_Knudsen_mdot(D, T_p, Psat, Re, mu_g, cp_g, lambda_g, P_g, R_g, Sc_g, R_v, Yinf) 

    #Hardcoded mdots for debugging
    #mdot = -1.1879e-6*math.sqrt(-1.5125e-9*p.time + 1.21e-6)  #Analytical mdot from miller's paper
    #mdot = -8.71137e-10
    
    #store mdot for output
    p.mdot = mdot

    #Hardcoded value for mdot for isolating energy equation from mass loss equation
    #mdot = -1.3e-9  #kg/s

    print('Computed mdot: {0:<8.4e}'.format(mdot))
    #//if(T_p < PsatF.get_Ttrip()) // Check if colder than triple point
    if T_p < 273.16: #Chris Hardcoded triple point for water Check if colder than triple point
      mdot = max(mdot, 0.0)  # only allow condensation

    # Compute blowing Reynolds Number
    Re_b = max(-mdot / (D * mu_g * math.pi), 0.0) 
    Pr_g = mu_g * cp_g / lambda_g 
    Nu = 2.0 + 0.552 * math.sqrt(Re) * Pr_g ** (1.0/3.0) 
    beta_b = 0.5 * Pr_g * Re_b 
    #Compute modification to heat transfer due to blowing
    if beta_b > 1.0e-6:
        f2 = beta_b / (math.exp(beta_b) - 1.0)
    else:
        f2 = 1.0

    mf = p.get_mass_fraction() 
    rho_p = p.get_density() 

    T1 = p.temp[0] 
    massp1 = dtbeta *mdot - alpha1 * p.mass[0] - alpha2 * p.mass[1] 

    if massp1 < 0.0:
        massp1 = 0.0
    elif massp1 < 0.01 * p.mass[0]:
      massp1 = 0 
    
    #Update mass of droplet
    p.mass[1] = p.mass[0] 
    p.mass[0] = massp1 
    #Update droplet diameter
    D = max(p.get_diameter(), 1e-10) 

    #double hvap = 288.7e3 //Chris heat of vaporization of decane at 420K J/kgK

    find_hi = False 
    hlo = 0 
    hhi = 1e100 
    #heatratio = p.mass[0]/max(mfluid*cp_g,1e-30) 
    oldfunc = 1e300 
    fudge2 = 1.0;
    fixedsrc = alpha1 * p.temp[0] + alpha2 * p.temp[1] - dtbeta * hvap * mdot / (max(p.mass[0], 1.0e-30) * p.get_c() * fudge2)
    tau_p = rho_p * D * D / (18.0 * mu_g) 
    print('rho_p: {0:<8.4f}\tD: {1:<8.4e}\tmu_g: {2:<8.4e}'.format(rho_p, D, mu_g) )
    fudge = 1.0 
    Ecoef = dtbeta * fudge * f2 * Nu * cp_g /(3.0 * Pr_g * tau_p * p.get_c()) 
    
    term_1 = fudge * f2 * Nu * cp_g / (3.0 * Pr_g * tau_p * p.get_c())
    term_2 = hvap * mdot / (max(p.mass[0], 1.0e-30) * p.get_c())
    format_string = 'f2: {0:<8.4f}   Nu: {1:<8.4f}   Cp_g: {2:<8.4f}   Pr: {3:<8.4f}   tau_p: {4:<8.4f}   ' \
                    'C_liquid: {5:<8.4f}   Term1Coeff: {6:<8.4e}   Term2Coeff: {7:<8.4e}   Beta: {8:<8.4f}\n' \
                    'Mass: {9:<8.4e}   Mdot: {10:<8.4e}   Xseq: {11:<8.4e}   Psat: {12:<8.4e}'

    print(format_string.format(f2, Nu, cp_g, Pr_g, tau_p, p.get_c(), term_1, term_2, beta_b, p.mass[0], p.mdot, Psat/P_g, Psat))


    max_iterations = 100 
    for i in range(max_iterations):
      print('T: {0:<10.4f}'.format(T1)) 
      
      # Update fluid temperature to consider convection heat transfer
      # Temperature of fluid must range between fluid and particle initial
      # temperatures to be physical

      #double ptemp = energyp.get_temp(h1,mf) ;
      ptemp = T1
      #print('ptemp: ', ptemp, '\th1: ', h1)
      #      maxtemp = max(ptemp,double(fluid_T)) 
      #      mintemp = min(ptemp,double(fluid_T)) 
      #      fluid_T_corr = min(maxtemp, max(mintemp, fluid_T - (h1 - p.h[0]) * heatratio)) 

      func = T1 + fixedsrc - Ecoef * (fluid_T - ptemp)  
      #print('func: ', func) 

      if i > 2 and ( abs(func) < 1e-6 * max(p.temp[0], T1) or hhi - hlo < 1e-5):
        break 
      if func < 0: 
        hlo = T1 
      
      if func > 0:
        hhi = T1 
        find_hi = True 
    
      T1 += -func


      factor = 0.5
      #Relax update to help convergence
      T1 = factor * T1 + (1.0 - factor) * ptemp 

      # Try to bracket if possible by overshooting increasing predictions if
      # upper bound not found
      if ~find_hi and  i > 4 and func < 0:
        T1 = 1.1*T1 
      # If out of bounds or not converging fast enough, switch to
      # bisection
#     if(find_hi && (h1 > max(hlo,hhi) || h1 < min(hlo,hhi)|| fabs(func/oldfunc) > .5)) {
      if find_hi and (T1 > max(hlo, hhi) or T1 < min(hlo, hhi)):
        # Use bisection if outside of known root bounds
        T1 = 0.5*(hlo + hhi) 

      oldfunc = func 


    if i > max_iterations -1: 
      print('particle temperature update failed to converge!')
    
 
    print('FluidT: {0:<8.4f}\tNew Droplet T: {1:<8.4f}'.format(fluid_T, T1))
    print('f2: {0:<8.4f}\tNu: {1:<8.4f}\tCp_g: {2:<8.4f}\tPr: {3:<8.4f}\ttau_p: {4:<8.4e}\tC_liquid: {5:<8.4f}\tD: {6:<8.4e}'.format(f2, Nu, cp_g, Pr_g, tau_p, p.get_c(), p.get_diameter())) 
    
    #These are the two terms in the energy equation
    debug_term_1 = (fudge * f2 * Nu * cp_g / (3.0 * Pr_g * tau_p * p.get_c())) * (fluid_T - T1)
    debug_term_2 = hvap * mdot / (p.mass[0] * p.get_c() * fudge2)
    print('Term 1: {0:<10.6e}\tTerm 2: {1:<10.6e}'.format(debug_term_1, debug_term_2))

    #Store debugging quantities
    p.Term1 = debug_term_1
    p.Term2 = debug_term_2

  
    p.temp[1] = p.temp[0]  # Update T^n-1
    p.temp[0] = T1    # Now T = T^n+1

    T_pn = p.get_temperature() 
    #nPsat = PsatF.get_Psat(T_pn) 
    nPsat = P_atm * math.exp((hvap / R_v) *(1.0 / T_B - 1.0 / T_pn))

    if nPsat > P_g:   # boiling limit temperature to boiling point
      print('Boiling temperature limiter active') 
      
      # Find boiling temperature at this pressure
      Tmax = T_pn 
      # Tmin = PsatF.get_Ttrip() 
      Tmin = 273.16 #Chris hardcoded triple point temperature of water 
      Tboil = 0.9 * Tmax + 0.1 * Tmin 
      max_iterations = 50 
      for i in range(max_iterations):
        # f = PsatF.get_Psat(Tboil, dPdT) - P_g 
        f = P_atm * math.exp((hvap / R_v) * (1.0 / T_B - 1.0 / Tboil)) - P_g  #Chris hardcoded saturation pressure at boiling temperature
        if abs(f) < 1e-4 * P_g:
          break 
        if f > 0.0:
          Tmax = Tboil 
        if f < 0.0: 
          Tmin = Tboil 
        
        Tboil -= f / dPdT 
        if Tboil > Tmax or Tboil < Tmin:
          Tboil = 0.5 * (Tmax + Tmin) 
      
      # print('boiling temperature switch activated, Tboil =', Tboil)
      mf = p.get_mass_fraction() 
      # hnew = energyp.get_h(Tboil,mf) 
      #reset temperature to boiling temperature
      p.temp[0] = Tboil 

    return mdot





def integrateMomentum(p, dt, fluid_v, fvolpp, mu_g, rfluid):
    """
    NewParticle &p,
    float dt,
    double fluid_v, //carrier gas velocity
    double fvolpp, // Volume allocated to the droplet
    double mu_g, //carrier gas viscosity
    double rfluid //carrier gas density
    """

    #integration constants
    beta = 2./3. 
    alpha1 = -4./3. 
    alpha2 = 1./3. 
    dtbeta = dt * beta 

    vel1 = p.vel[0] 
    pos1 = dtbeta * vel1 - alpha1 * p.pos[0] - alpha2 * p.pos[1] 
    rp = p.get_density() 
    D = p.get_diameter() 
    mdot = (p.mass[0] - p.mass[1]) / dt 
    
    mfluid = rfluid * fvolpp + 1e-30 # mass of fluid around particle
    fixedsrc = -alpha1 * p.vel[0] - alpha2 * p.vel[1]  
    volp = math.pi * D * D * D / 6.0 
    volpp = fvolpp 
    # enhance drag function for large volume fraction
    alphav = min(2.0, volp / max(volpp, 1e-30)) 
    
    fp_vf = max((8.0 * alphav) ** 6.0 - 0.001, 0.0) 

    #Integration loop
    max_iterations = 20
    for i in range(max_iterations): 
      #Update fluid velocity based on delta particle momentum
      if i > 0: #Past first iteration
        fluid_v = fluid_v - ((vel1 - p.vel[0]) * p.mass[0] / mfluid ) 

      dv = abs(fluid_v - vel1) 
      Re = rfluid * D * dv / mu_g 
      # blowing Reynolds number
      Re_b = abs(mdot / (D * mu_g * math.pi)) 
      a = 0.09 + 0.077 * math.exp(-0.4 * Re) 
      b = 0.4 + 0.77 * math.exp(-0.04 * Re) 
      denom = 1.0 + a * Re_b **b 

      fpblow = (1. + 0.0545 * Re + 0.1 * math.sqrt(Re) * (1.0 - 0.03 * Re)) / denom + fp_vf 
      # Clift-Gauvin drag function (Crowe, 1998)
      fpcg = 1.0 + 0.15 * Re ** 0.687 + 0.0175 * Re / (1.0 + 4.25e4 * (Re+1e-20) **-1.16) + fp_vf 
      # Choose drag function based on reynolds number.  For high reynolds
      # number use Clift Gauvin, otherwise use blowing reynolds number 
      if Re < 100:
        fp = fpblow
      else:
        fp = fpcg
      taup = rp * D ** 2 / (18.0 * mu_g * fp) 
      vcoef = dtbeta / taup 

      #      vel1 =  (vcoef*fluid_v + fixedsrc)/(1.+vcoef) 
      f = (vcoef * fluid_v + fixedsrc) / (1.0 + vcoef) - vel1 
      df = -vcoef * p.mass[0] / (mfluid * (1.0 + vcoef)) - 1.0 
      vel1 -= -f/df 
      pos1 = dtbeta * vel1 - alpha1 * p.pos[0] - alpha2 * p.pos[1] 

      # If iterated at least 2 times, check for convergence
      if i > 1 and abs(f) / (abs(df) * (0.1 + abs(vel1))) < 1.0e-5 : 
        break 
      
    # Now advance the particle momentum in time
    p.vel[2] = p.vel[1] 
    p.vel[1] = p.vel[0] 
    p.vel[0] = vel1 
    p.pos[1] = p.pos[0] 
    p.pos[0] = pos1 




if __name__ == "__main__":
    # execute only if run as a script
    main()
