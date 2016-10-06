FUNCTION derivs, x, y

common somename, tt, FRon, len, an
;------------------------------------------------------------------
;******************************************************************
; right hand side of field line ODE
;******************************************************************
;
;JT B-field
len=1e7
bscl=0.01
Lscl=1.0d0*len
dydx = dblarr(3)

; constants:
b0=1.0d0*bscl    ; 	B0 - eq field B strength
b1=20.0d0*b0 ; 	B1 - flux ring B strength

z0=5.0d0*lscl	; 	null position (+/-Z0)
xc=0.0d0*lscl   ; 	flux ring centre in x
yc=0.0d0*lscl  	;	    	    	    y    	    	    
zc=0.0d0*lscl 	;	    	    	    z
a=0.1d*z0   	; a, relates to radius of flux ring
oa=1.0d0/a   
ll=0.2d*z0
oll=1.0d0/ll	; l, relates to height
;oan=L/a   
;oln=L/z0	    ; l, relates to height


Bx=b0*y[0]*(y[2]-3*z0)/Lscl/Lscl
By=b0*y[1]*(y[2]+3*z0)/Lscl/Lscl
Bz=b0*0.5d*(z0*z0-y[2]*y[2]+y[1]*y[1]+y[0]*y[0])/Lscl/Lscl

 ;now the general form of the flux ring
fring=exp(-(y[0]-xc)*(y[0]-xc)*oa*oa-(y[1]-yc)*(y[1]-yc)*oa*oa-(y[2]-xc)*(y[2]-xc)*oll*oll)
Bfx=-2.0d0*b1*oa*(y[1]-yc)*fring
Bfy=2.0d0*b1*oa*(y[0]-xc)*fring
Bfz=0.0d0
 
; the field line tracker can be called at various positions and times.
; make it so that there is a distinction between including/not including flux ring for ANY time. 
 IF FRon THEN BEGIN
  ;print, 'on'	 
  dydx[0]=Bx+tt*Bfx
  dydx[1]=By+tt*Bfy
  dydx[2]=Bz+tt*Bfz
 ENDIF ELSE BEGIN
  ;print, 'off'
  dydx[0]=Bx
  dydx[1]=By
  dydx[2]=Bz
 ENDELSE 
 
return, dydx
end
;##################################################################
; Runge-Kutta routines with adaptive step size control taken
; from Numerical Recipes
; rk4 from idl lib
;##################################################################


PRO RKQC, y, dydx,x,htry,eps,yscal,hdid,hnext,string1

common somename, tt, FRon,len, an
;----------------------------------------------------
; *****************************************************************
;  RKQC is the stepper routine which makes one adaptive step 
;  forward
; *****************************************************************
;-----------------------------------------------------------------
; assign some fixed values to variables
;-----------------------------------------------------------------

pgrow = -0.2
pshrnk= -0.25
fcor  = 1.0/15.0
safety = 0.9
errcon = 6.0e-4
ytemp = dblarr(3)
ysav  = ytemp
dysav = ysav

;-----------------------------------------------------------------
; save initial values for step since it is needed twice
;-----------------------------------------------------------------

xsav=x
ysav=y
dysav=dydx

;-----------------------------------------------------------------
; try stepsize htry first and define hh as stepsize for the 
; two half-steps
;-----------------------------------------------------------------

h = htry
jump1: hh= 0.5 *h


;-----------------------------------------------------------------
; take two half-steps first
;-----------------------------------------------------------------

ytemp=rk4(ysav,dysav,xsav,hh,string1)

x=xsav+hh

dydx=derivs(x,ytemp)
y=rk4(ytemp,dydx,x,hh,string1)

;-----------------------------------------------------------------
; one step next
;-----------------------------------------------------------------

x=xsav+h
ytemp=rk4(ysav,dysav,xsav,h,string1)

;-----------------------------------------------------------------
; now evaluate the accuracy
;-----------------------------------------------------------------

errmax = 0.
ytemp = y-ytemp

errmax= max(abs(ytemp/yscal))

;-----------------------------------------------------------------
; scale error to required accuracy eps and reduce stepsize if 
; necessary
;-----------------------------------------------------------------

errmax = errmax/eps
if errmax gt 1.0 then begin
   h = safety*double(h)*(errmax^pshrnk)
   goto, jump1
endif else begin            
;-----------------------------------------------------------------
; if stepsize ok, compute stepsize for next step
;-----------------------------------------------------------------

   hdid = h                 
   if errmax gt errcon then begin
       hnext = safety*h*(errmax^pgrow)
   endif else begin
       hnext = 4.0*h
   endelse
endelse



;-----------------------------------------------------------------
; reduce truncation error to fifth order
;-----------------------------------------------------------------

y = y + ytemp*fcor

hnext=hnext

end ;_____of RKQC
;*****************************************************************
; ODEINT is the driver routine for the ODE solver
;*****************************************************************

PRO ODEINT, ystart,eps,h1,string1, yp

common somename, tt, FRon, len, an
;common path, kount;, yp
;maxsteps = 1e4-2
maxsteps=500
tiny = 1.0e-20

;-----------------------------------------------------------------
; we always start at arclength(=x) zero; set initial steplength to
; h1; and step counter to zero
;-----------------------------------------------------------------
x=0
h=h1
kount = 0
;-----------------------------------------------------------------
; assign initial values and store them in solution vector
;-----------------------------------------------------------------
y = ystart
yp=dblarr(3)
yp(*,0) = ystart(*)
;-----------------------------------------------------------------
; start integrating ODE
;-----------------------------------------------------------------
for iit=1,maxsteps do begin
	dydx = derivs(x,y)
;-----------------------------------------------------------------
; for monitoring accuracy; see Numerical Recipes for explanation
;-----------------------------------------------------------------
  	yscal = abs(y) + abs(h*dydx) + tiny
;-----------------------------------------------------------------
; make one adaptive step forward and write result to yp
;-----------------------------------------------------------------
  RKQC,y, dydx,x,h,eps,yscal,hdid,hnext,string1
  kount = kount + 1
  yp=[[yp],[dblarr(3)]]
  yp(*,kount) = y(*)
  IF ((ABS(y(0)) GT 10.0*len) OR (ABS(y(1)) GT 10.0*len) OR (ABS(y(2)) GT 7.5*len)) THEN BEGIN
   return
  ENDIF
  ;print, y(*)
  ;STOP
  h = hnext
  
endfor
end;____________of ODEINT
