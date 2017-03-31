

#Simple script to evaluate the Wilke rule for computing thermophysical properties

vapor_prop= 0.0188613
vapor_R=58.55 #J/kgK
gas_prop=0.0261682
gas_R=286.98
ref_mol_frac=0.0414

theta=vapor_R/gas_R

omega_vg = (((1+vapor_prop/gas_prop)**(1/2)*(theta)**(1/4))**2)/((8*(1+(1/theta)))**(1/2))
omega_gv = (((1+gas_prop/vapor_prop)**(1/2)*(1/theta)**(1/4))**2)/((8*(1+(theta)))**(1/2))

mixture_prop= (ref_mol_frac*vapor_prop)/(ref_mol_frac + (1-ref_mol_frac)*omega_vg) + ((1-ref_mol_frac)*gas_prop)/(ref_mol_frac*omega_gv + (1-ref_mol_frac))

print mixture_prop
