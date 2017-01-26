pro ptrack_field,pn,frame, global_v = global_v, local_ue = local_ue, epar = epar, rho = rho, eta = eta

p = getrdata(pn,/all)
f = getdata(frame)

maxjz=0
minjz=0
mycol=255
l=64
leta=11
al=16
bl=12
flev=2.8*findgen(l)/(l-1.)-2.8
tlev = findgen(l)/(l-1)*(1.205 - 0.45)+0.45
etalev = findgen(leta)/(leta-1.0)*0.01
;etalev = findgen(leta)/(leta-1.0)*2.0E-3 - 1.0E-3
eparlev = findgen(l)/(l-1.0)*2.0E-3-1.0E-3
rholev = findgen(l)/(l-1.0)*0.7 + 0.6
;eparlev = findgen(l)/(l-1.0)*2.0 - 1.0
;alev=9.3*findgen(al)/(al-1.)-9.3
alev=30.0*findgen(al)/(al-1.)-30.0
vlev=0.3*findgen(l)/(l-1.)
mblev=10^(1.2*findgen(bl)/(bl-1.)-0.7)
acol=alev
nvp=20
norm = 1.0E1

p1=[0.1,0.15,0.83,0.85]
p2=[0.87,0.15,0.9,0.85]

CURRENT=1
ETA_CONTOUR = 0
EPARALLEL = 0
RHO = 0
TEMP=0
VELOCITY=0
MYBETA=0

nx=f.grid.npts[0]
ny=f.grid.npts[1]
addaz, f

E_parallel = fltarr(nx-1,ny-1)

w = window(dimensions = [1400,800])

IF CURRENT THEN BEGIN
 f.jz(*,0)=0.0d0
 f.jz(*,ny-1)=0.0d0
 loadct, 10, /silent
 name_c1 = string(f.time,  format='("jz/az at ", f9.3)')
 colours = indgen(l)*(256./DOUBLE(l))
 c1 = contour(f.jz, f.grid.x, f.grid.y, /fill, ytitle = '$y$', xtitle = '$x$', title = name_c1,c_color = colours,c_value = flev,axis_style = 2,/overplot,position = p1,rgb_table = 1)
 c2 = contour(f.az, f.grid.x, f.grid.y, c_value = alev, /overplot, color = 'black',c_label_show = 0)
 pl1 = plot(p.x/norm,p.y/norm,/overplot,vert_colors = p.t/max(p.t)*255, rgb_table = 11,thick = 3)
 ;cb = colorbar(range = [min(f.jz),max(f.jz)], position = p2, title = 'jz')
 cb = colorbar(position = p2, title = 'jz', range = [min(flev),max(flev)], rgb_table = 1,orientation = 1,textpos = 1)
 loadct, 3, /silent
ENDIF
IF ETA_CONTOUR THEN BEGIN
 loadct, 10, /silent
 name_c1 = string(f.time,  format='("eta/az at ", f9.3)')
 colours = indgen(leta)*(256./DOUBLE(leta))
 c1 = contour(f.eta, f.x, f.y, /fill, ytitle = '$y$', xtitle = '$x$', title = name_c1,c_color = colours,c_value = etalev,axis_style = 2,/overplot,position = p1,rgb_table = 57)
 c2 = contour(f.az, f.grid.x, f.grid.y, c_value = alev, /overplot, color = 'black',c_label_show = 0)
 pl1 = plot(p.x/norm,p.y/norm,/overplot,vert_colors = p.t/max(p.t)*255, rgb_table = 11,thick = 3)
 cb = colorbar(range = [min(etalev),max(etalev)], position = p2, title = 'eta',taper = 0,rgb_table = 57,orientation = 1,textpos = 1)
 loadct, 3, /silent
ENDIF
IF EPARALLEL THEN BEGIN
 for j = 1,nx-2 do begin
  for k = 1,ny-2 do begin
   bx = (f.bx(j,k) + f.bx(j-1,k))/2.0
   by = (f.bx(j,k) + f.bx(j,k-1))/2.0
   bz = f.bz(j,k)
   vx = (f.vx(j,k) + f.vx(j-1,k) + f.vx(j,k-1) + f.vx(j-1,k-1))/4.0
   vy = (f.vy(j,k) + f.vy(j-1,k) + f.vy(j,k-1) + f.vy(j-1,k-1))/4.0
   vz = (f.vz(j,k) + f.vz(j-1,k) + f.vz(j,k-1) + f.vz(j-1,k-1))/4.0
   ;E = f.eta(j,k)*[f.jx(j,k),f.jy(j,k),f.jz(j,k)] - crossp([vx,vy,vz],[bx,by,bz])
   ;E_parallel(j-1,k-1) = (E(0)*bx + E(1)*by + E(2)*bz)/sqrt(bx*bx + by*by + bz*bz)
   E_parallel(j-1,k-1) = f.eta(j,k)*(f.jx(j,k)*bx + f.jy(j,k)*by + f.jz(j,k)*bz)/(bx*bx + by*by + bz*bz)
  endfor
 endfor
 loadct, 10, /silent
 name_c1 = string(f.time,  format='("E_parallel/az at ", f9.3)')
 colours = indgen(leta)*(256./DOUBLE(leta))

 rind = indgen(256)
 bind = indgen(256)
 bind = 2*(bind - 128)
 rind(128:255) = 0
 bind(0:127) = 0
 ctable = colortable([[rind,0,0],[0,0,bind],[rind*0,0,0]])

 print, 'max, min Epar = ', max(E_parallel), min(E_parallel)
 ;c1 = contour(E_parallel(1:nx-2,1:ny-2), f.x(1:nx-2), f.y(1:ny-2), /fill, ytitle = 'y', xtitle = 'x', title = name_c1,c_color = colours,c_value = eparlev,axis_style = 2,/overplot,position = p1,rgb_table = 57)
 c1 = contour(E_parallel(1:nx-2,1:ny-2), f.x(1:nx-2), f.y(1:ny-2), /fill, ytitle = '$y$', xtitle = '$x$', title = name_c1,axis_style = 2,/overplot,position = p1,rgb_table = 57)
 ;c1 = contour(E_parallel(1:nx-2,1:ny-2), f.x(1:nx-2), f.y(1:ny-2), /fill, ytitle = 'y', xtitle = 'x', title = name_c1,c_color = colours,c_value = eparlev,axis_style = 2,/overplot,position = p1,rgb_table = ctable)
 c2 = contour(f.az(1:nx-2,1:ny-2), f.grid.x(1:nx-2), f.grid.y(1:ny-2), c_value = alev, /overplot, color = 'black',c_label_show = 0)
 pl1 = plot(p.x/norm,p.y/norm,/overplot,vert_colors = p.t/max(p.t)*255, rgb_table = 11,thick = 3)
 cb = colorbar(range = [min(E_parallel),max(E_parallel)], position = p2, title = 'E parallel (V/m)',taper = 0,rgb_table = 57,orientation = 1,textpos = 1)
 ;cb = colorbar(range = [min(eparlev),max(eparlev)], position = p2, title = 'E parallel',taper = 0,rgb_table = 57,orientation = 1,textpos = 1)
 ;cb = colorbar(range = [min(eparlev),max(eparlev)], position = p2, title = 'E parallel',taper = 0,rgb_table = ctable,orientation = 1,textpos = 1)
 loadct, 3, /silent
ENDIF

IF RHO THEN BEGIN
 loadct, 10, /silent
 name_c1 = string(f.time,  format='("rho/az at ", f9.3)')
 colours = indgen(l)*(256./DOUBLE(l))
 c1 = contour(f.rho, f.x, f.y, /fill, ytitle = 'y', xtitle = 'x', title = name_c1,c_color = colours,c_value = rholev,axis_style = 2,/overplot,position = p1,rgb_table = 57)
 c2 = contour(f.az, f.grid.x, f.grid.y, c_value = alev, /overplot, color = 'black',c_label_show = 0)
 pl1 = plot(p.x/norm,p.y/norm,/overplot,vert_colors = p.t/max(p.t)*255, rgb_table = 11,thick = 3)
 cb = colorbar(range = [min(rholev),max(rholev)], position = p2, title = 'rho',taper = 0,rgb_table = 57,orientation = 1,textpos = 1)
 loadct, 3, /silent
ENDIF

newx=congrid(f.grid.x, nvp)
newy=congrid(f.grid.y, nvp)
newvx=congrid(f.vx, nvp,nvp)
newvy=congrid(f.vy, nvp,nvp)

;print, max(newvx)/mvx, max(newvy)/mvy
;print, where(newjz eq max(newjz))
;maxjz=[maxjz,max(f.jz)]
 minjz=[minjz,min(f.jz)]


if keyword_set(global_v) then v1 = vector(newvx,newvy, newx, newy, /overplot)
p_len = (size(p.x))(1)
vec_skip = 100
if keyword_set(local_ue) then begin
 v2 = vector(reform(p.UE(0,1:p_len - 1:vec_skip)),reform(p.UE(1,1:p_len - 1:vec_skip)),p.x(1:p_len - 1:vec_skip)/norm,p.y(1:p_len - 1:vec_skip)/norm,data_location = 0,/overplot)
 v3 = vector(reform(p.UE(0,1:p_len - 1:vec_skip)),reform(p.UE(1,1:p_len - 1:vec_skip)),p.x(1:p_len - 1:vec_skip)/norm,p.y(1:p_len - 1:vec_skip)/norm,data_location = 0, position = [0.13,0.1,0.9,0.87],xtitle = '$x$', ytitle = '$y$')
 pl1 = plot(p.x/norm,p.y/norm,/overplot,vert_colors = p.t/max(p.t)*255, rgb_table = 11,thick = 3)
 c1 = colorbar(range = [min(p.t),max(p.t)],position = [0.13,0.95,0.9,0.98],title = 'time',rgb_table = 11)
 ue = sqrt(p.ue(0,*)^2 + p.ue(1,*)^2 + p.ue(2,*)^2)
 pl2 = plot(p.x/norm,p.y/norm,vert_colors = ue/max(ue)*255, rgb_table = 11,thick = 3)
 c2 = colorbar(range = [min(ue),max(ue)],position = [0.13,0.95,0.9,0.98],title = '$E \times B$ drift velocity (m/s)',rgb_table = 11)
endif
if keyword_set(epar) then begin
 pl2 = plot(p.x/norm,p.y/norm,vert_colors = bytscl(p.epar), rgb_table = 11,thick = 3,position = [0.13,0.1,0.9,0.87],xtitle = '$x$', ytitle = '$y$')
 c2 = colorbar(range = [min(p.epar),max(p.epar)],position = [0.13,0.95,0.9,0.98],title = 'E parallel (V/m)')
endif
if keyword_set(rho) then begin
 pl2 = plot(p.x/norm,p.y/norm,vert_colors = bytscl(p.rho), rgb_table = 11,thick = 3,position = [0.13,0.1,0.9,0.87],xtitle = '$x$', ytitle = '$y$')
 c2 = colorbar(range = [min(p.rho),max(p.rho)],position = [0.13,0.95,0.9,0.98],title = 'rho')
endif
if keyword_set(eta) then begin
 pl2 = plot(p.x/norm,p.y/norm,vert_colors = bytscl(p.eta), rgb_table = 11,thick = 3,position = [0.13,0.1,0.9,0.87],xtitle = '$x$', ytitle = '$y$')
 if (min(p.eta) eq max(p.eta)) then begin
   c2 = colorbar(range = [min(p.eta),max(p.eta)+1.0],position = [0.13,0.95,0.9,0.98],title = 'eta')
 endif else begin
   c2 = colorbar(range = [min(p.eta),max(p.eta)],position = [0.13,0.95,0.9,0.98],title = 'eta')
 endelse
endif

end
