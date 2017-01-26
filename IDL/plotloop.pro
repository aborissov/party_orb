pro plotloop,first,last

	for j = first-1,last-1 do begin
		ds = getrdata(j+1,/all)
		;plot,ds.t,ds.theta,yrange = [0,!pi],title = j+1
		;wait,0.3
		window,0
		plot,ds.x,ds.y,title = j+1,xrange = [-1E2,1E2],yrange = [-2.2E2,2.2E2]
		window,2
		plot,ds.t,ds.Upar,title = 'Vpar'
		;window,3
		;plot,alog10(ds.t),ds.epar,title = 'eparallel'
		;plot,ds.t,ds.epar,title = j+1
		;wait,0.2
		;plot,ds.t,ds.epar,title = j+1
		;oplot,ds.t,(ds.upar - min(ds.upar))/max(ds.upar-min(ds.upar))*(max(ds.epar - min(ds.epar))) + min(ds.epar),color = 100
		wait,0.4
	end

end 
