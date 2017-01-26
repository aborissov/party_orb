pro epar,f

nx = (size(f.x))(1)
ny = (size(f.y))(1)
mjx = fltarr(nx,ny)
mjy = fltarr(nx,ny)
mjz = fltarr(nx,ny)
mvx = fltarr(nx,ny)
mvy = fltarr(nx,ny)
mvz = fltarr(nx,ny)
current = fltarr(nx,ny)
mbx = fltarr(nx,ny)
mby = fltarr(nx,ny)
eparallel = fltarr(nx,ny)

for j = 0,nx-1 do begin
   for k = 0,ny-1 do begin
      mjx(j,k) = (f.jx(j,k) + f.jx(j+1,k) + f.jx(j,k+1) + f.jx(j+1,k+1))/4.0
      mjy(j,k) = (f.jy(j,k) + f.jy(j+1,k) + f.jy(j,k+1) + f.jy(j+1,k+1))/4.0
      mjz(j,k) = (f.jz(j,k) + f.jz(j+1,k) + f.jz(j,k+1) + f.jz(j+1,k+1))/4.0
      current(j,k) = sqrt(mjx(j,k)*mjx(j,k) + mjy(j,k)*mjy(j,k) + mjz(j,k)*mjz(j,k))
      mvx(j,k) = (f.vx(j,k) + f.vx(j+1,k) + f.vx(j,k+1) + f.vx(j+1,k+1))/4.0
      mvy(j,k) = (f.vy(j,k) + f.vy(j+1,k) + f.vy(j,k+1) + f.vy(j+1,k+1))/4.0
      mvz(j,k) = (f.vz(j,k) + f.vz(j+1,k) + f.vz(j,k+1) + f.vz(j+1,k+1))/4.0
   endfor
endfor

for j = 0,nx-1 do begin
   mbx(j,*) = (f.bx(j+1,*) + f.bx(j,*))/2.0
endfor
for j = 0,ny-1 do begin
   mby(*,j) = (f.by(*,j+1) + f.by(*,j))/2.0
endfor

for j = 0,nx-1 do begin
   eparallel(j,*) = f.eta(j,*)*(mjx(j,*)*mbx(j,*) + mjy(j,*)*mby(j,*) + mjz(j,*)*f.bz(j,*))
endfor

c1 = contour(eparallel,f.x,f.y,/fill,rgb_table = 3, position = [0.1,0.1,0.9,0.88],xtitle = 'x $\times 10^5 m$',ytitle = 'y $\times 10^5 m$')
cb = colorbar(range = [min(eparallel),max(eparallel)],rgb_table = 3,title = 'epar', position = [0.1,0.89,0.9,0.93],textpos = 1)


end
