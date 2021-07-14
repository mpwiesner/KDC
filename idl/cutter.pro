pro cutter

;This routine reads in an image and produces a cutout of it.
 
a=mrdfits('STRIPE82/drC-006580-g3-0577.fits', 0, header)
hextract, a, header, 725, 855, 1299, 1414
writefits, 'IMAGENEW.fits', a ,header
 
 
end
