pro header_write


h = headfits('/Users/mwiesner/DC2_out/woof.fits') 
;sxaddpar,h,'CTYPE1','RA--TAN' 
;modfits,'/Users/mwiesner/DC2_out/woof.fits',0,h 
sxdelpar, h, ['RADECSYS']
;sxaddpar, h, 'RADECSYS','J2000','RA-DEC system '
modfits,'/Users/mwiesner/DC2_out/woof.fits',0,h 

end
