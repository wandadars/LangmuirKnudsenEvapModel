import matplotlib.pyplot as plt
import numpy as np
import pylab

pylab.rcParams['agg.path.chunksize'] =  500000
dt = 1.0e-3

datalist = pylab.loadtxt('data.txt')
N=len(datalist)
print(N)


time = datalist[:,0] 

pylab.plot( time, datalist[:,1] ,'*')
pylab.xlabel("Time")
pylab.ylabel("Temperature")
#pylab.ylim([273,298])

pylab.savefig('temperature.png',bbox_inches='tight')


plt.close()
pylab.plot( time, datalist[:,2] )
pylab.xlabel("Time")
pylab.ylabel("Diameter")
pylab.savefig('diameter.png',bbox_inches='tight')


plt.close()
pylab.plot( time, datalist[:,2]**2 )
pylab.xlabel("Time")
pylab.ylabel("D^2 (Squared Diameter)")
pylab.savefig('diameter_squared.png',bbox_inches='tight')


plt.close()
pylab.plot( time, -datalist[:,5] )
pylab.xlabel("Time")
pylab.ylabel("mdot")
pylab.savefig('mdot.png',bbox_inches='tight')


plt.close()
pylab.plot( time, datalist[:,6], label="Term1" )
pylab.plot( time, datalist[:,7], label="Term2")
pylab.xlabel("Time")
pylab.ylabel("Term Size")
pylab.legend()
#pylab.ylim([-10,10])
pylab.savefig('energy_equation_terms.png',bbox_inches='tight')


plt.close()
pylab.plot( time, datalist[:,8], label="RHS of Energy Equation" )
pylab.xlabel("Time")
pylab.ylabel("Term Size")
pylab.legend()
#pylab.ylim([-0.2,0.2])
pylab.savefig('net_equation.png',bbox_inches='tight')


plt.close()
pylab.plot( time, datalist[:,9], label="Droplet Mass" )
pylab.xlabel("Time")
pylab.ylabel("Mass(kg)")
pylab.legend()
#pylab.ylim([-0.2,0.2])
pylab.savefig('mass.png',bbox_inches='tight')
