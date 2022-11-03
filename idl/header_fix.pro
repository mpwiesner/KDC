pro header_fix

;This routine reads in a single PhoSim image and fixes the header parameters so 
;it can display properly in DS9.

for i=10, 20 DO BEGIN
	filename = "~/Desktop/woof/lsst_"+strtrim(strcompress(i),2)+"_kn.fits"

	h = headfits(filename) 
	sxdelpar, h, ['RADESYS']
;sxaddpar, h, 'RADESYS ','J2000','RA-DEC system '

;sxdelpar, h, ['CTYPE1']
;sxaddpar,h,'CTYPE1','RA---TAN' 
;sxaddpar,h,'CUNIT1','deg' 
;sxaddpar,h,'CUNIT2','deg' 

;sxdelpar, h, ['RADESYS']

;sxdelpar, h, ['CUNIT1']
;sxdelpar, h, ['CUNIT2']
	sxdelpar, h, ['CTYPE1']
	sxaddpar,h,'CTYPE1','RA---TAN' 

	modfits,filename,0,h 
	endfor

end
