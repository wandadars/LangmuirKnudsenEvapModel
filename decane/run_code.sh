#!/bin/bash

gcc -o dropet_evap_test integrate_particle_temperatureForm.cc -lstdc++ -lm 
rm -rf out.txt  data.txt
./dropet_evap_test >& out.txt 
python plot_data.py 
xdg-open temperature.png 
xdg-open diameter_squared.png

