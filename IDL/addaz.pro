PRO addaz,data

on_error, 2

close,1

IF N_ELEMENTS(data) NE 0 THEN BEGIN

    nx=data.grid.npts[0]
    ny=data.grid.npts[1]

    dyb = shift(data.grid.y,-1) - data.grid.y
    dxb = shift(data.grid.x,-1) - data.grid.x

    unit = 0.0D

    Az=dblarr(nx,ny)
    ;Az=data.bx

    Az[*] = unit

    FOR iy = 1, ny-2 DO BEGIN
        Az[0,iy] = data.bx[0,iy] * dyb[iy] + Az[0,iy-1]
    ENDFOR

    FOR ix = 1, nx-2 DO BEGIN
        Az[ix,0] = Az[ix-1,0] - data.by[ix,0] * dxb[ix] 
    ENDFOR

    FOR iy = 1, ny-2 DO BEGIN
        FOR ix = 1, nx-2 DO BEGIN
            Az[ix,iy] = data.bx[ix,iy] * dyb[iy] + Az[ix,iy-1]
        ENDFOR
    ENDFOR

data=CREATE_STRUCT(data,{az: az})

ENDIF ELSE BEGIN
    print, "Invalid variable passed"
    print, "Use: addcurrents, <data sturcture>"
ENDELSE

END


