function AB_hist,x,nb,w

;h = histogram(x,nbins = nb,locations = hbins, reverse_indices = R)
h = histogram(x,binsize = 0.1,locations = hbins, reverse_indices = R)
wh = h
nh = (size(h))(1)

for j = 0,nh-2 do begin
  if (R(j) ne R(j+1)) then begin
    wh(j) = wh(j) + total(w(R(R(j):R(j+1)-1))-1)
  endif
endfor

hist = {whist:wh,whbins:hbins}
return,hist

end
