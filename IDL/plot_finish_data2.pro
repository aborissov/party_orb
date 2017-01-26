; 2 VARIABLES
pro plot_finish_data2,df1,df2

c = 3.0E8
m = 9.11E-31
q = 1.6E-19
gamma1 = df1.ke_ev/(m*c*c)*q+1.0
gamma2 = df2.ke_ev/(m*c*c)*q+1.0
vtot1 = sqrt(c*c - c*c/(gamma1*gamma1))
vtot2 = sqrt(c*c - c*c/(gamma2*gamma2))
zbound1 = 1.1*max(abs(df1.z))
zbound2 = 1.1*max(abs(df2.z))
zbound = max(zbound1,zbound2)
vbound1 = 1.1*max(abs(df1.vpar))
vbound2 = 1.1*max(abs(df2.vpar))
vbound = max(vbound1,vbound2)
thetaboundmin1 = 1.1*min(df1.theta)
thetaboundmin2 = 1.1*min(df2.theta)
thetaboundmin = min(thetaboundmin1,thetaboundmin2)
thetaboundmax1 = 1.1*max(df1.theta)
thetaboundmax2 = 1.1*max(df2.theta)
thetaboundmax = max(thetaboundmax1,thetaboundmax2)
betaboundmin1 = 1.1*min(cos(df1.theta))
betaboundmin2 = 1.1*min(cos(df2.theta))
betaboundmin = min(betaboundmin1,betaboundmin2)
betaboundmax1 = 1.1*max(cos(df1.theta))
betaboundmax2 = 1.1*max(cos(df2.theta))
betaboundmax = max(betaboundmax1,betaboundmax2)

n1 = 80

name1 = '$\lambda_0 = 10^4 m$'		
name2 = '$\lambda_0 = 10^2 m$'		

pdfz1 = histogram(df1.z,locations = z1bins,max = zbound,min = -zbound,binsize = zbound/n1)
pdfv1 = histogram(df1.vpar,locations = v1bins,max = vbound,min = -vbound,binsize = vbound/n1)
pdfth1 = histogram(df1.theta,locations = th1bins,nbins = n1,min = 0.0, max = !pi)
pdfbt1 = histogram(cos(df1.theta),locations = bt1bins,nbins = n1,min = -1.0, max = 1.0)
pdfvt1 = histogram(vtot1,locations = vt1bins,nbins = n1)
pdfg1 = histogram(df1.ke_ev/1.0E3,locations = g1bins,nbins = n1)

pdfz2 = histogram(df2.z,locations = z2bins,max = zbound,min = -zbound,binsize = zbound/n1)
pdfv2 = histogram(df2.vpar,locations = v2bins,max = vbound,min = -vbound,binsize = vbound/n1)
pdfth2 = histogram(df2.theta,locations = th2bins,nbins = n1,min = 0.0, max = !pi)
pdfbt2 = histogram(cos(df2.theta),locations = bt2bins,nbins = n1,min = -1.0, max = 1.0)
pdfvt2 = histogram(vtot2,locations = vt2bins,nbins = n1)
pdfg2 = histogram(df2.ke_ev/1.0E3,locations = g2bins,nbins = n1)

;z1 = plot(z1bins,pdfz1,stairstep = 1,thick = 2,xtitle = 'z',ytitle = 'count',font_size = 16,name = name1)
;z2 = plot(z2bins,pdfz2,stairstep = 1,thick = 2,name = name2,/overplot,color = 'red')
;lz = legend(target = [z1,z2])

;v1 = plot(v1bins,pdfv1,stairstep = 1,thick = 2,xtitle = 'v_par',ytitle = 'count',font_size = 16,name = name1)
;v2 = plot(v2bins,pdfv2,stairstep = 1,thick = 2,name = name2,/overplot,color = 'red')
;lv = legend(target = [v1,v2])

;th1 = plot(th1bins,pdfth1,stairstep = 1,thick = 2,xtitle = 'theta',ytitle = 'count',font_size = 16,name = name1)
;th2 = plot(th2bins,pdfth2,stairstep = 1,thick = 2,name = name2,/overplot,color = 'red')
;lth = legend(target = [th1,th2])

bt1 = plot(bt1bins,pdfbt1,stairstep = 1,thick = 2,xtitle = 'beta',ytitle = 'count',font_size = 16,name = name1)
bt2 = plot(bt2bins,pdfbt2,stairstep = 1,thick = 2,name = name2,/overplot,color = 'red')
lbt = legend(target = [bt1,bt2])

;vt1 = plot(vt1bins,pdfvt1,stairstep = 1,thick = 2,xtitle = 'vtot',ytitle = 'count',font_size = 16,name = name1)
;vt2 = plot(vt2bins,pdfvt2,stairstep = 1,thick = 2,name = name2,/overplot,color = 'red')
;lvt = legend(target = [vt1,vt2])

g1 = plot(g1bins,pdfg1,stairstep = 1,thick = 2,xtitle = 'Energy (keV)',ytitle = 'count',font_size = 16,name = name1)
g2 = plot(g2bins,pdfg2,stairstep = 1,thick = 2,name = name2,/overplot,color = 'red')
lg = legend(target = [g1,g2])

kb = 8.5E-5
Temp = 1.0E6
E = g1bins
;y = g1bins*exp(-g1bins/(kb*Temp)*1.0E3)
y = 2.0/(sqrt(!pi))*sqrt(E)/((kb*Temp)*sqrt(kb*Temp))*exp(-E/(kb*Temp)*1.0E3)
y = y/max(y)*max(pdfg1)

;gc = plot(g1bins,y,thick = 2, color = 'blue', /overplot)

end
