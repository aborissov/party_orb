COMMON wkdirs, wkdir_global, tt, FRon, len
@ IDL/StartCFD.pro
;.r IDL/getdata.pro  	    ; reads in data
;.r IDL/getndata.pro  	    ; reads in nonrelativistic data
.r IDL/getrdata.pro  	    ; reads in relativistic data
.r IDL/iniorbplot.pro	    ; plot initial array of positions
.r IDL/legend.pro   	    ; creates quick legend
.r IDL/mk_vector.pro	    ; makes pointy 3d arrows
.r IDL/ODEINT.pro   	    ; ODE integrator
.r IDL/orb__define.pro	    ; used to define 3d orbs
.r IDL/particletrack.pro    ; for 1 particle highlights track in 3d
;.r IDL/pvels.pro    	    ; creates 3d vector plot
;.r IDL/remove.pro   	    ; removes elements from array
;.r IDL/regionzoom.pro	    ; like particletrack but zoomable
.r IDL/symbol_obj.pro	    ; ?
.r IDL/quickbifur.pro	    ; a program to (coarsely) look for zeros in 1d
.r IDL/StartLARE.pro
;.r ODEINT
.r fieldwrapping
.r IDL/plotloop.pro
.r IDL/get_finish_data.pro
.r IDL/plot_finish_data1.pro
.r IDL/plot_finish_data2.pro
@ StartSDF.pro
.r IDL/ptrack_field.pro
.r IDL/AB_hist.pro


Q0 = 1.60217646d-19 ; proton charge [C]
M0 = 9.10938188d-31 ; electron mass [kg]
kb = 1.3806503d-23  ; Boltzmann's constant [J/K]
wkdir_global="Data"

.r trilinear_3D tracefield_3D 
.r destaggerB
.r xplot3dJT



