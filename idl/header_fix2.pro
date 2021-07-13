pro header_fix2
; This edits the headers in PhoSim FITS images so they display with WCS properly in DS9.

name = ['lsst_e_159479_f1_R01_S00_E000.fits','lsst_e_159479_f1_R01_S01_E000.fits', 'lsst_e_159479_f1_R01_S02_E000.fits', 'lsst_e_159479_f1_R01_S10_E000.fits', $
'lsst_e_159479_f1_R01_S11_E000.fits', 'lsst_e_159479_f1_R01_S12_E000.fits', 'lsst_e_159479_f1_R01_S20_E000.fits', 'lsst_e_159479_f1_R01_S21_E000.fits',$
'lsst_e_159479_f1_R01_S22_E000.fits', 'lsst_e_159479_f1_R02_S00_E000.fits', 'lsst_e_159479_f1_R02_S01_E000.fits', 'lsst_e_159479_f1_R02_S02_E000.fits', 'lsst_e_159479_f1_R02_S10_E000.fits',$
'lsst_e_159479_f1_R02_S11_E000.fits', 'lsst_e_159479_f1_R02_S12_E000.fits', 'lsst_e_159479_f1_R02_S20_E000.fits', 'lsst_e_159479_f1_R02_S21_E000.fits']

for i=0, 16 do begin
	print, 'working on ', name[i]
	h = headfits(name[i]) 
	sxdelpar, h, ['RADESYS']
	sxdelpar, h, ['CTYPE1']
	sxaddpar,h,'CTYPE1','RA---TAN' 
	modfits,name[i],0,h 
endfor

end