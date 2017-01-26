; 1 VARIABLE
pro plot_finish_data1,df1

c = 3.0E8
m = 9.11E-31
q = 1.6E-19
gamma = df1.ke_ev/(m*c*c)*q+1.0
vtot = sqrt(c*c - c*c/(gamma*gamma))
zbound = 1.1*max(abs(df1.z))
vbound = 1.1*max(abs(df1.vpar))
thetaboundmin = 1.1*min(df1.theta)
thetaboundmax = 1.1*max(df1.theta)
betaboundmin = 1.1*min(cos(df1.theta))
betaboundmax = 1.1*max(cos(df1.theta))

print,zbound,vbound
n1 = 80

name1 = '$B = 1$'

;pdfz1 = histogram(df1.z,locations = z1bins,max = zbound,min = -zbound,binsize = zbound/n1)
;pdfv1 = histogram(df1.vpar,locations = v1bins,max = vbound,min = -vbound,binsize = vbound/n1)
;pdfth = histogram(df1.theta,locations = thbins,nbins = n1,min = 0.0, max = !pi)
pdfbt = histogram(cos(df1.theta),locations = btbins,nbins = n1,min = -1.0, max = 1.0)
;pdfvt = histogram(vtot,locations = vtbins,nbins = n1)
;pdfg = histogram(df1.ke_ev/1.0E3,locations = gbins,nbins = n1)
pdfgl = histogram(alog10(df1.ke_ev/1.0E3),locations = glbins,nbins = n1)

;z1 = plot(z1bins,pdfz1,stairstep = 1,thick = 2,xtitle = 'z',ytitle = 'count',font_size = 16,name = name1)
;v1 = plot(v1bins,pdfv1,stairstep = 1,thick = 2,xtitle = 'v_par',ytitle = 'count',font_size = 16,name = name1)
;th = plot(thbins,pdfth,stairstep = 1,thick = 2,xtitle = 'theta',ytitle = 'count',font_size = 16,name = name1)
bt = plot(btbins,pdfbt,stairstep = 1,thick = 2,xtitle = 'beta',ytitle = 'count',font_size = 16,name = name1)
;vt = plot(vtbins,pdfvt,stairstep = 1,thick = 2,xtitle = 'vtot',ytitle = 'count',font_size = 16,name = name1)
;g = plot(gbins,pdfg,stairstep = 1,thick = 2,xtitle = 'Energy (keV)',ytitle = 'count',font_size = 16,name = name1)
gl = plot(glbins,pdfgl,stairstep = 1,thick = 2,xtitle = 'log(Energy (keV))',ytitle = 'count',font_size = 16,name = name1)

kb = 8.5E-5
Temp = 1.0E6
E = gbins
;y = E*exp(-E/(kb*Temp)*1.0E3)
y = 2.0/(sqrt(!pi))*sqrt(E)/((kb*Temp)*sqrt(kb*Temp))*exp(-E/(kb*Temp)*1.0E3)
y = y/max(y)*max(pdfg)

;gc = plot(E,y,thick = 2, color = 'red', /overplot)

end
