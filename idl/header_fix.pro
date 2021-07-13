pro header_fix


h = headfits('lsst_e_159479_f1_R02_S21_E000.fits') 
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

modfits,'lsst_e_159479_f1_R02_S21_E000.fits',0,h 

end