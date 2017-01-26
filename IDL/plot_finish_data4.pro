; 2 VARIABLES
pro plot_finish_data4,df1,df2,df3,df4

c = 3.0E8
m = 9.11E-31
q = 1.6E-19
gamma1 = df1.ke_ev/(m*c*c)*q+1.0
gamma2 = df2.ke_ev/(m*c*c)*q+1.0
gamma3 = df3.ke_ev/(m*c*c)*q+1.0
gamma4 = df4.ke_ev/(m*c*c)*q+1.0
vtot1 = sqrt(c*c - c*c/(gamma1*gamma1))
vtot2 = sqrt(c*c - c*c/(gamma2*gamma2))
vtot3 = sqrt(c*c - c*c/(gamma3*gamma3))
vtot4 = sqrt(c*c - c*c/(gamma4*gamma4))
zbound1 = 1.1*max(abs(df1.z))
zbound2 = 1.1*max(abs(df2.z))
zbound3 = 1.1*max(abs(df3.z))
zbound4 = 1.1*max(abs(df4.z))
zbound = max([zbound1,zbound2,zbound3,zbound4])
vbound1 = 1.1*max(abs(df1.vpar))
vbound2 = 1.1*max(abs(df2.vpar))
vbound3 = 1.1*max(abs(df3.vpar))
vbound4 = 1.1*max(abs(df4.vpar))
vbound = max([vbound1,vbound2,vbound3,vbound4])
thetaboundmin1 = 1.1*min(df1.theta)
thetaboundmin2 = 1.1*min(df2.theta)
thetaboundmin3 = 1.1*min(df3.theta)
thetaboundmin4 = 1.1*min(df4.theta)
thetaboundmin = min([thetaboundmin1,thetaboundmin2,thetaboundmin3,thetaboundmin4])
thetaboundmax1 = 1.1*max(df1.theta)
thetaboundmax2 = 1.1*max(df2.theta)
thetaboundmax3 = 1.1*max(df3.theta)
thetaboundmax4 = 1.1*max(df4.theta)
thetaboundmax = max([thetaboundmax1,thetaboundmax2,thetaboundmax3,thetaboundmax4])
betaboundmin1 = 1.1*min(cos(df1.theta))
betaboundmin2 = 1.1*min(cos(df2.theta))
betaboundmin3 = 1.1*min(cos(df3.theta))
betaboundmin4 = 1.1*min(cos(df4.theta))
betaboundmin = min([betaboundmin1,betaboundmin2,betaboundmin3,betaboundmin4])
betaboundmax1 = 1.1*max(cos(df1.theta))
betaboundmax2 = 1.1*max(cos(df2.theta))
betaboundmax3 = 1.1*max(cos(df3.theta))
betaboundmax4 = 1.1*max(cos(df4.theta))
betaboundmax = max([betaboundmax1,betaboundmax2,betaboundmax3,betaboundmax4])

n1 = 80

name1 = '$\lambda_0 = 10^6 m$'		
name2 = '$\lambda_0 = 10^8 m$'		
name3 = '$\lambda_0 = 10^4 m$'		
name4 = '$\lambda_0 = 10^2 m$'		

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

pdfz3 = histogram(df3.z,locations = z3bins,max = zbound,min = -zbound,binsize = zbound/n1)
pdfv3 = histogram(df3.vpar,locations = v3bins,max = vbound,min = -vbound,binsize = vbound/n1)
pdfth3 = histogram(df3.theta,locations = th3bins,nbins = n1,min = 0.0, max = !pi)
pdfbt3 = histogram(cos(df3.theta),locations = bt3bins,nbins = n1,min = -1.0, max = 1.0)
pdfvt3 = histogram(vtot3,locations = vt3bins,nbins = n1)
pdfg3 = histogram(df3.ke_ev/1.0E3,locations = g3bins,nbins = n1)

pdfz4 = histogram(df4.z,locations = z4bins,max = zbound,min = -zbound,binsize = zbound/n1)
pdfv4 = histogram(df4.vpar,locations = v4bins,max = vbound,min = -vbound,binsize = vbound/n1)
pdfth4 = histogram(df4.theta,locations = th4bins,nbins = n1,min = 0.0, max = !pi)
pdfbt4 = histogram(cos(df4.theta),locations = bt4bins,nbins = n1,min = -1.0, max = 1.0)
pdfvt4 = histogram(vtot4,locations = vt4bins,nbins = n1)
pdfg4 = histogram(df4.ke_ev/1.0E3,locations = g4bins,nbins = n1)

;z1 = plot(z1bins,pdfz1,stairstep = 1,thick = 2,xtitle = 'z',ytitle = 'count',font_size = 16,name = name1)
;z2 = plot(z2bins,pdfz2,stairstep = 1,thick = 2,name = name2,/overplot,color = 'red')
;z3 = plot(z3bins,pdfz3,stairstep = 1,thick = 2,name = name3,/overplot,color = 'green')
;z4 = plot(z4bins,pdfz4,stairstep = 1,thick = 2,name = name4,/overplot,color = 'purple')
;lz = legend(target = [z2,z1,z3,z4])

;v1 = plot(v1bins,pdfv1,stairstep = 1,thick = 2,xtitle = 'v_par',ytitle = 'count',font_size = 16,name = name1)
;v2 = plot(v2bins,pdfv2,stairstep = 1,thick = 2,name = name2,/overplot,color = 'red')
;v3 = plot(v3bins,pdfv3,stairstep = 1,thick = 2,name = name3,/overplot,color = 'green')
;v4 = plot(v4bins,pdfv4,stairstep = 1,thick = 2,name = name4,/overplot,color = 'purple')
;lv = legend(target = [v2,v1,v3,v4])

;th1 = plot(th1bins,pdfth1,stairstep = 1,thick = 2,xtitle = 'theta',ytitle = 'count',font_size = 16,name = name1)
;th2 = plot(th2bins,pdfth2,stairstep = 1,thick = 2,name = name2,/overplot,color = 'red')
;th3 = plot(th3bins,pdfth3,stairstep = 1,thick = 2,name = name3,/overplot,color = 'green')
;th4 = plot(th4bins,pdfth4,stairstep = 1,thick = 2,name = name4,/overplot,color = 'purple')
;lth = legend(target = [th2,th1,th3,th4])

bt1 = plot(bt1bins,pdfbt1,stairstep = 1,thick = 2,xtitle = 'beta',ytitle = 'count',font_size = 16,name = name1)
bt2 = plot(bt2bins,pdfbt2,stairstep = 1,thick = 2,name = name2,/overplot,color = 'red')
bt3 = plot(bt3bins,pdfbt3,stairstep = 1,thick = 2,name = name3,/overplot,color = 'green')
bt4 = plot(bt4bins,pdfbt4,stairstep = 1,thick = 2,name = name4,/overplot,color = 'purple')
lbt = legend(target = [bt2,bt1,bt3,bt4])

;vt1 = plot(vt1bins,pdfvt1,stairstep = 1,thick = 2,xtitle = 'vtot',ytitle = 'count',font_size = 16,name = name1)
;vt2 = plot(vt2bins,pdfvt2,stairstep = 1,thick = 2,name = name2,/overplot,color = 'red')
;vt3 = plot(vt3bins,pdfvt3,stairstep = 1,thick = 2,name = name3,/overplot,color = 'green')
;vt4 = plot(vt4bins,pdfvt4,stairstep = 1,thick = 2,name = name4,/overplot,color = 'purple')
;lvt = legend(target = [vt2,vt1,vt3,vt4])

g1 = plot(g1bins,pdfg1,stairstep = 1,thick = 2,xtitle = 'Energy (keV)',ytitle = 'count',font_size = 16,name = name1)
g2 = plot(g2bins,pdfg2,stairstep = 1,thick = 2,name = name2,/overplot,color = 'red')
g3 = plot(g3bins,pdfg3,stairstep = 1,thick = 2,name = name3,/overplot,color = 'green')
g4 = plot(g4bins,pdfg4,stairstep = 1,thick = 2,name = name4,/overplot,color = 'purple')
lg = legend(target = [g2,g1,g3,g4])

kb = 8.5E-5
Temp = 1.0E6
E = g1bins
;y = g1bins*exp(-g1bins/(kb*Temp)*1.0E3)
y = 2.0/(sqrt(!pi))*sqrt(E)/((kb*Temp)*sqrt(kb*Temp))*exp(-E/(kb*Temp)*1.0E3)
y = y/max(y)*max(pdfg1)

;gc = plot(g1bins,y,thick = 2, color = 'blue', /overplot)

end
