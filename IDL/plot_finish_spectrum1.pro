function edist,energy,temperature
  m = 5.11E5
  kb = 8.6E-5
  np = 100.0
  ;return, np*(m/(2.0*!pi*kb*temperature))^1.5*8.0*!pi*energy/m*exp(-energy/(kb*temperature))
  return, np*2.0/(sqrt(!pi))*sqrt(energy)/((kb*temperature)*sqrt(kb*temperature))*exp(-energy/(kb*temperature))
end

function betadist,theta
  np = 100.0
  return, np/sin(theta)
end

; 1 VARIABLE
pro plot_finish_spectrum1,df1,ds1

c = 3.0E8
m = 9.11E-31
q = 1.6E-19
Temp = 1.0E6
kb = 8.6E-5
gamma = df1.ke_ev/(m*c*c)*q+1.0
thetaboundmin = 1.1*min(df1.theta)
thetaboundmax = 1.1*max(df1.theta)
betaboundmin = 1.1*min(cos(df1.theta))
betaboundmax = 1.1*max(cos(df1.theta))

np = (size(ds1.x))(1)
weights = fltarr(np)

;print, size(weights), size(edist(ds1.ke_ev(*),Temp)), size(betadist(ds1.theta(*)))
weights(*) = edist(ds1.ke_ev(*),Temp)*betadist(ds1.theta(*))
;h = AB_hist(alog10(df1.ke_ev/1.0E3),100,weights(0:np-2))
h = AB_hist(df1.ke_ev/1.0E3,100,weights(0:np-2))

;p = plot(h.whbins,h.whist,stairstep = 1, thick = 2, xtitle = 'Energy(keV)', ytitle = count, font_size = 16)
p = plot(alog10(h.whbins),alog10(h.whist), thick = 2, xtitle = 'Energy(keV)', ytitle = count, font_size = 16)
;print, alog10(h.whbins), h.whbins

;E = h.whbins
;;y = E*exp(-E/(kb*Temp)*1.0E3)
;y = 2.0/(sqrt(!pi))*sqrt(E)/((kb*Temp)*sqrt(kb*Temp))*exp(-E/(kb*Temp)*1.0E3)
;y = y/max(y)*max(h.whist)
;
;gc = plot(E,y,thick = 2, color = 'red', /overplot)

;print,zbound,vbound
;n1 = 80
;
;name1 = '$B = 1$'
;
;pdfbt = histogram(cos(df1.theta),locations = btbins,nbins = n1,min = -1.0, max = 1.0)
;pdfg = histogram(df1.ke_ev/1.0E3,locations = gbins,nbins = n1)
;pdfgl = histogram(alog10(df1.ke_ev/1.0E3),locations = glbins,nbins = n1)
;
;bt = plot(btbins,pdfbt,stairstep = 1,thick = 2,xtitle = 'beta',ytitle = 'count',font_size = 16,name = name1)
;;g = plot(gbins,pdfg,stairstep = 1,thick = 2,xtitle = 'Energy (keV)',ytitle = 'count',font_size = 16,name = name1)
;gl = plot(glbins,pdfgl,stairstep = 1,thick = 2,xtitle = 'log(Energy (keV))',ytitle = 'count',font_size = 16,name = name1)
;
;kb = 8.5E-5
;E = gbins
;;y = E*exp(-E/(kb*Temp)*1.0E3)
;y = 2.0/(sqrt(!pi))*sqrt(E)/((kb*Temp)*sqrt(kb*Temp))*exp(-E/(kb*Temp)*1.0E3)
;y = y/max(y)*max(pdfg)
;
;;gc = plot(E,y,thick = 2, color = 'red', /overplot)

end
