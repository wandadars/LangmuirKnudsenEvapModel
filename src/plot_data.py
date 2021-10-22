import matplotlib.pyplot as plt
import os
import numpy as np
import pylab

pylab.rcParams['agg.path.chunksize'] =  500000

datalist = pylab.loadtxt('data.txt', comments='#')
N=len(datalist)
print('Number of timestep data: ' + str(N))

current_dir = os.getcwd()

output_path = os.path.join(current_dir, 'output')
print('Outputting plots to: ' + output_path)
os.makedirs(output_path, exist_ok=True)

time = datalist[:,0] 
print('Max Time: ' + str(max(time)))


temperature = datalist[:,1]
plot_file_name = 'temperature.png'
output_filename = os.path.join(output_path, plot_file_name)
pylab.plot( time, temperature, '*')
pylab.xlabel("Time")
pylab.ylabel("Temperature")
pylab.ylim([min(temperature), max(temperature)])
pylab.savefig(output_filename, bbox_inches='tight')
plt.close()

diameter = datalist[:,2]
plot_file_name = 'diameter.png'
output_filename = os.path.join(output_path, plot_file_name)
pylab.plot( time, diameter )
pylab.xlabel("Time")
pylab.ylabel("Diameter")
pylab.savefig(output_filename, bbox_inches='tight')
plt.close()

plot_file_name = 'diameter_squared.png'
output_filename = os.path.join(output_path, plot_file_name)
pylab.plot( time, 1e6*datalist[:,2]**2 )
pylab.xlabel("Time")
pylab.ylabel("D^2 (Squared Diameter)")
pylab.savefig(output_filename, bbox_inches='tight')
plt.close()

plot_file_name = 'mdot.png'
output_filename = os.path.join(output_path, plot_file_name)
pylab.plot( time, -datalist[:,5] )
pylab.xlabel("Time")
pylab.ylabel("mdot")
pylab.savefig(output_filename, bbox_inches='tight')
plt.close()

plot_file_name = 'energy_equation_terms.png'
output_filename = os.path.join(output_path, plot_file_name)
pylab.plot( time, datalist[:,6], label="Term1" )
pylab.plot( time, datalist[:,7], label="Term2")
pylab.xlabel("Time")
pylab.ylabel("Term Size")
pylab.legend()
pylab.ylim([-10, 10])
pylab.savefig(output_filename, bbox_inches='tight')
plt.close()

plot_file_name = 'net_equation.png'
output_filename = os.path.join(output_path, plot_file_name)
pylab.plot( time, datalist[:,8], label="RHS of Energy Equation" )
pylab.xlabel("Time")
pylab.ylabel("Term Size")
pylab.legend()
pylab.ylim([-0.2, 0.2])
pylab.savefig(output_filename, bbox_inches='tight')
plt.close()

plot_file_name = 'mass.png'
output_filename = os.path.join(output_path, plot_file_name)
pylab.plot( time, datalist[:,9], label="Droplet Mass" )
pylab.xlabel("Time")
pylab.ylabel("Mass(kg)")
pylab.legend()
pylab.savefig(output_filename, bbox_inches='tight')
