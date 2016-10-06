PRO FLINE
; FLINE
;  numerically integrates field lines according to field specified in derivs
;  then displays
;  REQUIRES - ODEINT, derivs, RKQC

@ODEINT
common somename, tt, FRon, len, an

 ;thisPalette = Obj_New('IDLgrPalette')
 ;thisPalette->LoadCT, 34
 lscl=1e6
 zscl=1e6
 dims=[1.,1.,1.]
 ;mysymsize=[lscl,lscl,lscl]*0.1*dims
 mysymsize=[1.0,1.0,1.0]*7.*dims  


 oOrb = OBJ_NEW('orb', COLOR=[255, 0 ,0]) 
; oOrb2 = OBJ_NEW('orb', COLOR=[255, 0 ,0]) 
 oOrb->Scale, dims[0]*5., dims[1]*5., dims[2]*5. 
; oOrb2->Scale, .04, .04, .25 
 oSymb = OBJ_NEW('IDLgrSymbol', oOrb)     ;oSymbol is green orb for start point
; oSymbol2 = OBJ_NEW('IDLgrSymbol', oOrb2)   ;osymbol2 is red orb for stop point


 oModel=obj_new('IDLgrModel')
 oModel->add,symbol_obj(2,[0,0,255],[0,0,0])
 oTet = obj_new('IDLgrSymbol', oModel)
 IF Obj_Valid(oTet) THEN oTet->SetProperty, Size=mysymsize 

tt=0.0d
;tim=0.3

n1=[ [[0,0]*lscl],[[0,0]*lscl],[[-50,-50]*zscl]]
n2=[ [[0,0]*lscl],[[0,0]*lscl],[[50,50]*zscl]]
sep=[ [[0,0]*lscl],[[0,0]*lscl],[[-50,50]*zscl]]


;;----- define variables ----------------
eps=0.0000001
h1=0.01     	    ; step size?
string1='derivs'
FRon=1	    	    ;switch flux ring on?

;;----- pick initial conditions ------------
ystart=dblarr(3)
r1=dblarr(3)
r2=dblarr(3)
rsteps=intarr(3)


r1[0]=-4.0e7	;xstart
r2[0]=4.0e7 	;xend
rsteps[0]=5 	;nx
r1[1]=-1.0e7	;ystart
r2[1]=1.0e7 	;yend
rsteps[1]=5	;ny
r1[2]=-8e7  	;zstart
r2[2]=8e7   	;zend
rsteps[2]=5 	;nz

;r1[0]=-1.0e3	;xstart
;r2[0]=1.0e3 	;xend
;rsteps[0]=4 	;nx
;r1[1]=-1.0e3	;ystart
;r2[1]=1.0e3 	;yend
;rsteps[1]=4	;ny
;r1[2]=-3e7  	;zstart
;r2[2]=3e7   	;zend
;rsteps[2]=7 	;nz



IF (rsteps[0] EQ 1) THEN  gdx=1.0d0 ELSE gdx=1.0d0/FLOAT(rsteps[0]-1)
IF (rsteps[1] EQ 1) THEN  gdy=1.0d0 ELSE gdy=1.0d0/FLOAT(rsteps[1]-1)
IF (rsteps[2] EQ 1) THEN  gdz=1.0d0 ELSE gdz=1.0d0/FLOAT(rsteps[2]-1)
gds=[gdx,gdy,gdz]
lbox=[r2(0)-r1(0),r2(1)-r1(1),r2(2)-r1(2)]

pno=0
np=rsteps[0]*rsteps[1]*rsteps[2]

flinecol=[0,0,255]
   ;x_start=0.15
   ;y_start=0.0
   ;z_start=0.0
   ;ystart=[, y_start, z_start]

FOR ix=0,rsteps[0]-1 DO BEGIN
 FOR iy=0,rsteps[1]-1 DO BEGIN
  FOR iz=0,rsteps[2]-1 DO BEGIN
   
   ;gds=[gdx,gdy,gdz]
   i=[ix,iy,iz]
   ;yorig= R1+lbox*(i*1.0D0)*gds
   ystart= R1+lbox*(i*1.0D0)*[gdx,gdy,gdz]
   yp=ystart
   ODEINT, ystart, eps, h1, string1, yp
   
   n_data=size(yp)
      
     ts=[ [yp[0,0],yp[0,0]],[yp[1,0],yp[1,0]],[yp[2,0],yp[2,0]]]
     te=[ [yp[0,n_data[2]-1],yp[0,n_data[2]-1]],[yp[1,n_data[2]-1],yp[1,n_data[2]-1]],[yp[2,n_data[2]-1],yp[2,n_data[2]-1]]]
     oi=[ [yp[0,n_data[2]-1]-yp[0,n_data[2]-2]],[yp[1,n_data[2]-1]-yp[1,n_data[2]-2]],[yp[2,n_data[2]-1]-yp[2,n_data[2]-2]]]
   
   oModel=obj_new('IDLgrModel')
   oModel->add,mk_vector([oi[0],oi[1],oi[2]]/sqrt(oi[0]*oi[0]+oi[1]*oi[1]+oi[2]*oi[2])/dims,color=flinecol)
   oAR3 = obj_new('IDLgrSymbol', oModel)
   IF Obj_Valid(oAR3) THEN oAR3->SetProperty, Size=mysymsize
    
   
    IF pno eq 0 THEN BEGIN
     xplot3d, yp(0,*)/lscl, yp(1,*)/lscl, yp(2,*)/zscl, xrange=[-100,100], yrange=[-100,100], zrange=[-100,100], COLOR=flinecol, $
     ztitle='z (Mm)', xtitle='x (Mm)',ytitle='y (Mm)'
    ENDIF ELSE BEGIN
     xplot3d, yp(0,*)/lscl, yp(1,*)/lscl, yp(2,*)/zscl, COLOR=[0,0,255], /OVERPLOT
    ENDELSE
   ; XPLOT3D, ts[*,0], ts[*,1], ts[*,2], COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol, THICK=2, /OVERPLOT
   
   ; XPLOT3D, te[*,0], te[*,1], te[*,2], COLOR=[0,0,0], NAME='end', SYMBOL=oSymbol2, THICK=2, /OVERPLOT
     xplot3D, te[*,0]/lscl, te[*,1]/lscl, te[*,2]/zscl, COLOR=[0,0,0], NAME='atest', SYMBOL=oAR3, THICK=2, /OVERPLOT
    pno=pno+1
  ENDFOR
 ENDFOR
ENDFOR

xplot3D, n1[*,0]/lscl, n1[*,1]/lscl, n1[*,2]/zscl, COLOR=[0,0,0], NAME='atest', SYMBOL=oSymb, THICK=2, /OVERPLOT
xplot3D, n2[*,0]/lscl, n2[*,1]/lscl, n2[*,2]/zscl, COLOR=[0,0,0], NAME='atest', SYMBOL=oSymb, THICK=2, /OVERPLOT
xplot3D, sep[*,0]/lscl, sep[*,1]/lscl, sep[*,2]/zscl, COLOR=[0,0,0], NAME='atest', THICK=4, /OVERPLOT

;;----3D Plot of Field Lines----------------------------------------------

;xplot3d, yp(0,*), yp(1,*), yp(2,*), xrange=[-1,1], yrange=[-1,1], zrange=[-6,6] 


end
