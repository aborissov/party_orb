function get_finish_data, n_c=n_c, wkdir=wkdir, start=start, exit=exit
on_error, 2

IF n_elements(wkdir) eq 0 THEN wkdir="./"
IF ((n_elements(n_c) eq 0) and (n_elements(start) eq 0)) THEN n_c=10
IF ((n_elements(n_c) eq 0) and (n_elements(start) ne 0)) THEN n_c=9



;files=FINDFILE(wkdir+'RV*.dat',count=count)
;IF count eq 0 THEN BEGIN
; PRINT, "ERROR: no RV*.dat files found"
; return, 0
;ENDIF
;IF (ii lt 1) or (ii gt count) THEN BEGIN
; print, "ERROR: data files from ii=1->"+STRING(COUNT,FORMAT='(i4.4)')$
; +string(ii,format='("; ",i4.4," outside range!")')
; return, 0
;ENDIF

if (n_elements(exit) eq 0) then begin
  ;filename=string(wkdir,format='(a,"finishr.tmp")')
  ;tstring=string(wkdir,format='("wc -l ",a,"finishr.tmp")')
  if n_elements(start) ne 0 then begin
    filename=string(wkdir,format='(a,"start.tmp")')
    tstring=string(wkdir,format='("wc -l ",a,"start.tmp")')
  endif else begin
    filename=string(wkdir,format='(a,"finish.tmp")')
    tstring=string(wkdir,format='("wc -l ",a,"finish.tmp")')
  endelse
  spawn,tstring,res
  n_data=1UL
  dum='blah'
  reads,res,n_data,dum
  A=dblarr(n_c,n_data)
  openr,lun,filename,/get_lun
  readf,lun,A
  Free_lun,lun
  if n_c eq 10 then begin
    index=transpose(A(0,*))
    x=transpose(A(1,*))
    y=transpose(A(2,*))
    z=transpose(A(3,*))
    vpar=transpose(A(4,*))
    vperp2=transpose(A(5,*))
    ke_ev = transpose(A(6,*))
    theta = transpose(A(7,*))
    rho = transpose(A(8,*))
    tfinal = transpose(A(9,*))
  endif else begin
    index=transpose(A(0,*))
    x=transpose(A(1,*))
    y=transpose(A(2,*))
    z=transpose(A(3,*))
    vpar=transpose(A(4,*))
    vperp2=transpose(A(5,*))
    ke_ev = transpose(A(6,*))
    theta = transpose(A(7,*))
    rho = transpose(A(8,*))
  endelse
 
  case 1 of
    n_c eq 10: thedata={name: filename, index:index,x:x,y:y,z:z,vpar:vpar,vperp2:vperp2,ke_ev:ke_ev,theta:theta,rho:rho,tfinal:tfinal}
    n_c eq 9: thedata={name: filename, index:index,x:x,y:y,z:z,vpar:vpar,vperp2:vperp2,ke_ev:ke_ev,theta:theta,rho:rho}
  endcase

endif else begin
  filename=string(wkdir,format='(a,"exit.tmp")')
  tstring=string(wkdir,format='("wc -l ",a,"exit.tmp")')
 
  n_c = 2
  spawn,tstring,res
  n_data=1UL
  dum='blah'
  reads,res,n_data,dum
  A=dblarr(n_c,n_data)
  openr,lun,filename,/get_lun
  readf,lun,A
  Free_lun,lun
 
  pn=transpose(A(0,*))
  status=transpose(A(1,*))
  
  thedata={name: filename, pn:pn,status:status}
endelse

return, thedata

END
