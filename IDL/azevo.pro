;new program to watch the simple evolution of az in l2d runs
@JTblue_red

nt=40
;loadct, 3
maxjz=0
minjz=0
mycol=255
l=64
al=16
bl=12
flev=2.8*findgen(l)/(l-1.)-2.8
tlev = findgen(l)/(l-1)*(1.205 - 0.45)+0.45
;alev=0.5*findgen(al)/(al-1.)-0.5
;alev=reverse((1.-alog10((10.^2.0-10.^1.0)*findgen(al)/(al-1.)+10.0^0)-1)/4.7)
alev=9.3*findgen(al)/(al-1.)-9.3
vlev=0.3*findgen(l)/(l-1.)
mblev=10^(1.2*findgen(bl)/(bl-1.)-0.7)
acol=alev
nvp=20

p1=[0.15,0.15,0.7,0.9]
p2=[0.15,0.85,0.9,0.9]
window, 0
!p.background=255
;JTblue_red, pos=0.75
;loadct, 10, /silent


CURRENT=1
TEMP=0
VELOCITY=0
MYBETA=0
 ;mr=0
 ;FOR i=0, nt DO BEGIN
 ; ds=getdata(i,/vx,/vy)
 ; mlen=sqrt(ds.vx*ds.vx+ds.vy*ds.vy)
 ; IF (max(mlen) gt mr) THEN mr=max(mlen)
 ;ENDFOR

 
FOR i=0, nt DO BEGIN
 ds=getdata(i)
 nx=ds.grid.npts[0]
 ny=ds.grid.npts[1]
 addaz, ds

 ;shade_surf, ds.jz
 ;wait, 0.3
 ;continue
 IF CURRENT THEN BEGIN
  ds.jz(*,0)=0.0d0
  ds.jz(*,ny-1)=0.0d0
  loadct, 10, /silent
  CONTOUR, ds.jz, ds.grid.x, ds.grid.y, /fill, ytitle='y', xtitle='x', title=string(ds.time,  format='("jz/az at ", f9.3)'), $
  position=p1, /iso, c_colors=indgen(l)*(256./DOUBLE(l)), levels=flev, col=0
  CONTOUR, ds.az, ds.grid.x, ds.grid.y, levels=alev, /noerase, position=p1, col=0, /iso, c_colors=0
  loadct, 3, /silent
 ENDIF
 IF TEMP THEN BEGIN
  ;print,max(ds.temperature),min(ds.temperature)
  ds.temperature(*,0)=0.0d0
  ds.temperature(*,ny-2)=0.0d0
  loadct, 10, /silent
  CONTOUR, ds.temperature, ds.grid.x(0:nx-2), ds.grid.y(0:ny-2), /fill, ytitle='y', xtitle='x', title=string(ds.time,  format='("temperature/az at ", f9.3)'), $
  position=p1, /iso, c_colors=indgen(l)*(256./DOUBLE(l)), levels=tlev, col=0
  CONTOUR, ds.az, ds.grid.x, ds.grid.y, levels=alev, /noerase, position=p1, col=0, /iso, c_colors=0
  loadct, 3, /silent
 ENDIF
 ;modv=sqrt(ds.vx*ds.vx+ds.vy*ds.vy+ds.vz*ds.vz)
 ;print, minmax(modv)
 ;continue
 IF VELOCITY THEN BEGIN
  modv=sqrt(ds.vx*ds.vx+ds.vy*ds.vy+ds.vz*ds.vz)
  loadct, 3, /silent
  CONTOUR, modv, ds.grid.x, ds.grid.y, /fill, ytitle='y', xtitle='x', title=string(ds.time,  format='("|v| at ", f9.3)'), $
  position=p1, /iso, c_colors=indgen(l)*(256./DOUBLE(l)), levels=vlev, col=0
  CONTOUR, ds.az, ds.grid.x, ds.grid.y, levels=alev, /noerase, position=p1, col=0, /iso, c_colors=255
  loadct, 8, /silent
 ENDIF
 IF MYBETA THEN BEGIN
  mbeta=2.0d0*ds.pressure(0:nx-2,0:nx-2)/(ds.bx(0:nx-2,0:nx-2)*ds.bx(0:nx-2,0:nx-2)+ds.by(0:nx-2,0:nx-2)*ds.by(0:nx-2,0:nx-2)+ds.bz(0:nx-2,0:nx-2)*ds.bz(0:nx-2,0:nx-2))
  loadct, 15, /silent
  CONTOUR, mbeta, ds.grid.x(0:nx-2), ds.grid.y(0:nx-2), /fill, ytitle='y', xtitle='x', title=string(ds.time,  format='("beta at ", f9.3)'), $
  position=p1, /iso, c_colors=indgen(bl)*(256./DOUBLE(bl)), levels=mblev, col=0
  CONTOUR, ds.az, ds.grid.x, ds.grid.y, levels=alev, /noerase, position=p1, col=0, /iso, c_colors=255
  loadct, 8, /silent
 ENDIF
 newx=congrid(ds.grid.x, nvp)
 newy=congrid(ds.grid.y, nvp)
 newvx=congrid(ds.vx, nvp,nvp)
 newvy=congrid(ds.vy, nvp,nvp)
 
 ;print, max(newvx)/mvx, max(newvy)/mvy
 ;print, where(newjz eq max(newjz))
 ;maxjz=[maxjz,max(ds.jz)]
; minjz=[minjz,min(ds.jz)]
  
 
 velovect, newvx,newvy, newx, newy, /overplot,  position=p1, /iso, col=150

 ;IF CURRENT THEN BEGIN
 ; loadct, 10, /silent
 ; cgCOLORBAR, ncolors=mycol+1, /vertical, position=p2, $
 ; format='(f9.2)', RANGE=[min(flev),max(flev)];, $
 ; write_png, string(i,format='("img/j",I3.3,".png")'), TVRD(/TRUE)
 ;ENDIF
 ;IF VELOCITY THEN BEGIN
 ; loadct, 3, /silent
 ; cgCOLORBAR, ncolors=mycol+1, /vertical, position=p2, $
 ; format='(f9.2)', RANGE=[min(vlev),max(vlev)];, $
 ; write_png, string(i,format='("img/v",I3.3,".png")'), TVRD(/TRUE)
 ;ENDIF
 ;IF MYBETA THEN BEGIN
 ; loadct, 15, /silent
 ; cgCOLORBAR, ncolors=mycol+1, /vertical, position=p2, $
 ; format='(e8.1)', RANGE=[min(mblev),max(mblev)], /ylog, divisions=bl/2;, $
 ; ;write_png, string(i,format='("img/j",I3.3,".png")'), TVRD(/TRUE)
 ;ENDIF 

 ;IF CURRENT THEN BEGIN
 ; loadct, 10, /silent
 ; COLORBAR, ncolors=mycol+1, /vertical, position=p2, $
 ; format='(f9.2)', RANGE=[min(flev),max(flev)];, $
 ; write_png, string(i,format='("img/j",I3.3,".png")'), TVRD(/TRUE)
 ;ENDIF
 ;IF VELOCITY THEN BEGIN
 ; loadct, 3, /silent
 ; COLORBAR, ncolors=mycol+1, /vertical, position=p2, $
 ; format='(f9.2)', RANGE=[min(vlev),max(vlev)];, $
 ; write_png, string(i,format='("img/v",I3.3,".png")'), TVRD(/TRUE)
 ;ENDIF
 ;IF MYBETA THEN BEGIN
 ; loadct, 15, /silent
 ; COLORBAR, ncolors=mycol+1, /vertical, position=p2, $
 ; format='(e8.1)', RANGE=[min(mblev),max(mblev)], /ylog, divisions=bl/2;, $
 ; ;write_png, string(i,format='("img/j",I3.3,".png")'), TVRD(/TRUE)
 ;ENDIF 
wait, 1.0

ENDFOR


END
