import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from scipy import interpolate
import math
import sys
from pathlib import Path
import h5py
import bisect
from astropy.cosmology import FlatLambdaCDM
import argparse
import os.path


# Main method for calculating synthetic AB magnitudes.

def calc_abmag(args):

    #  Extract the bandList...
    bandList = args.bandList
    bandList = bandList.split(',')
    if args.verbose > 0:
        print('bandList: ', bandList)

    #  Extract the name of the bandpassFile...
    bandpassFile = args.bandpassFile
    if os.path.isfile(bandpassFile)==False:
        print("""bandpassFile %s does not exist...""" % (bandpassFile))
        print('Returning with error code 1 now...')
        return 1
    if args.verbose > 0:
        print('bandpassFile: ', bandpassFile)

    #  Extract the name of the spectrum file...
    spectrumFile = args.spectrumFile
    if os.path.isfile(spectrumFile)==False:
        print("""spectrumFile %s does not exist...""" % (spectrumFile))
        print('Returning with error code 1 now...')
        return 1
    if args.verbose > 0:
        print('spectrumFile: ', spectrumFile)

    # Try to determine spectrumFile type (FITS file or CSV file)...
    spectrumType = 'Unknown'
    try:
        hdulist = fits.open(spectrumFile)
        hdulist.close()
        spectrumType = 'FITS'
    except IOError:
        if args.verbose > 2:
            print("""spectrumFile %s is not a FITS file...""" % (spectrumFile))
        try:
            df_test = pd.read_csv(spectrumFile)
            spectrumType = 'CSV'
        except IOError:
            if args.verbose > 2:
                print("""spectrumFile %s is not a CSV file...""" % (spectrumFile))

    # Read in spectrumFile and create a SciPy interpolated function of the spectrum...
    if spectrumType == 'FITS':
        flux,wave_lo,wave_hi = getSpectrumSynphot(spectrumFile, fluxFactor=1.0)
    elif spectrumType == 'CSV':
        flux,wave_lo,wave_hi = getSpectrumCSV(spectrumFile, 
                                              colname_wave=args.colname_wave, colname_flam=args.colname_flam, 
                                              fluxFactor=1.0)
    else:
        print("""Spectrum file %s is of unknown type...""" % (spectrumFile))
        print('Returning with error code 1 now...')
        return 1

    # Read the bandpassFile into a Pandas DataFrame...
    df_resp = pd.read_csv(bandpassFile, comment='#')

    # Check to make sure the spectrumFile covers at least the same wavelength range
    #  as the bandpassFile...
    if ( (wave_lo > df_resp['LAMBDA'].min()) or (wave_hi < df_resp['LAMBDA'].max()) ):
        print("""WARNING:  %s does not cover the full wavelength range of %s""" % (spectrumFile, bandpassFile))
        print('Returning with error code 1 now...')
        return 1

    # Create wavelength_array and flux_array...
    delta_wavelength = 1.0 # angstroms
    wavelength_array = np.arange(wave_lo, wave_hi, delta_wavelength)
    flux_array = flux(wavelength_array)
    

    # If needed, convert flux from flam to fnu...
    if args.flux_type == 'Fnu':
        fnu_array = flux_array
    elif args.flux_type == 'Flam':
        c_kms = 299792.5        # speed of light in km/s
        c_ms = 1000.*c_kms      # speed of light in m/s
        c_as = (1.000e10)*c_ms  # speed of light in Angstroms/sec
        fnu_array = flux_array * wavelength_array * wavelength_array / c_as
    else:
        print("""Flux type %s is unknown...""" % (args.flux_type))
        print('Returning with error code 1 now...')
        return 1


    # Print out header...
    outputLine = ''
    for band in bandList:
        outputLine =  """%s,%s""" % (outputLine, band)
  #  print(outputLine[1:])

    outputLine = ''
    for band in bandList:

        response = interpolate.interp1d(df_resp['LAMBDA'], df_resp[band],
                                        bounds_error=False, fill_value=0.,
                                        kind='linear')
        response_array = response(wavelength_array)

        try:
            abmag = calc_abmag_value(wavelength_array, response_array, fnu_array)
        except Exception:
            abmag = -9999.

        outputLine =  """%s,%.4f""" % (outputLine, abmag)

   # print(outputLine[1:])
    f.write(str(outputLine[1:])+","+"{:.2f}".format(t_kn)+"\n")
    return 0

# Create an argparse Namespace and run "calc_abmag(args)"...
#  (e.g., run_calc_abmag('u,g,r,i,z,Y', bandpassFile, tempFile, 'LAMBDA', 'Flam', 'Flam', verbose) )

def run_calc_abmag(bandList, bandpassFile, spectrumFile, colname_wave, colname_flam, flux_type, verbose):
    
    args = argparse.Namespace(bandList = bandList, 
                              bandpassFile = bandpassFile, 
                              spectrumFile = spectrumFile, 
                              colname_wave = colname_wave, 
                              colname_flam = colname_flam, 
                              flux_type = flux_type, 
                              verbose = verbose)
    
    if args.verbose > 0: print(args)

    status = calc_abmag(args)
    
     
    return status

  
#parser.add_argument('--bandList', help='comma-separated list with no spaces', default='g,r,i,z,Y')
#parser.add_argument('--bandpassFile', help='name of the input plan file', default='DES_STD_BANDPASSES_Y3A2_ugrizY.test.csv')
#parser.add_argument('--spectrumFile', help='name of the input plan file (can be CSV file or a synphot-style FITS file')
#parser.add_argument('--colname_wave', help='name of the wavelength column (in case of a CSV spectrumFile)', default='wave')
#parser.add_argument('--colname_flux', help='name of the flux column (in case of a CSV spectrumFile)', default='flux')
#parser.add_argument('--flux_type', help='type of flux (Flam [ergs/sec/cm**2/Angstrom] or Fnu [ergs/sec/cm**2/Hz])? ', default='Flam')
#parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    

# Calculate abmag using the wavelength version of the Fukugita et al. (1996) equation.
def calc_abmag_value(wavelength_array, response_array, fnu_array):

    # Calculate the abmag...
    numerator = np.sum(fnu_array * response_array / wavelength_array)
    denominator = np.sum(response_array / wavelength_array)
    abmag_value = -2.5*math.log10(numerator/denominator) - 48.60

    return abmag_value


# Return a SciPy interpolation function of a Synphot-style FITS spectrum.
#  (Based on code from Keith Bechtol's synthesize_locus.py.)
# Unless otherwise noted, fluxes are assumed to be Flam and wavelengths
#  are assumed to be in Angstroms.

def getSpectrumSynphot(synphotFileName, fluxFactor=1.0):

    try:

        hdulist = fits.open(synphotFileName)
        t = Table.read(hdulist[1])
        hdulist.close()

    except IOError:

        print("""Could not read %s""" % synphotFileName)
        sys.exit(1)


    wave = t['WAVELENGTH'].data.tolist()
    wave_lo = min(wave)
    wave_hi = max(wave)
    t['FLUX'] = fluxFactor*t['FLUX']
    flam = t['FLUX'].data.tolist()
    flam = t['FLUX'].data.tolist()
    data = {'wavelength': wave, 'flux': flam}

    f = interpolate.interp1d(data['wavelength'], data['flux'],
                             bounds_error=True,
                             kind='linear')

    return f,wave_lo,wave_hi


# Return a SciPy interpolation function of a CSV-style spectrum.
#  (Based on code from Keith Bechtol's synthesize_locus.py.)
# Unless otherwise noted, fluxes are assumed to be Flam and wavelengths
#  are assumed to be in Angstroms.

def getSpectrumCSV(csvFileName, colname_wave='wave', colname_flam='flux', fluxFactor=1.0):

    try:

        df = pd.read_csv(csvFileName)

    except IOError:

        print("""Could not read %s""" % csvFileName)
        sys.exit(1)


    columnNameList = df.columns.tolist()

    if colname_wave not in columnNameList:
        print("""Column %s not in %s""" % (colname_wave, csvFileName))
        sys.exit(1)

    if colname_flam not in columnNameList:
        print("""Column %s not in %s""" % (colname_wave, csvFileName))
        sys.exit(1)

    wave = df[colname_wave].tolist()
    wave_lo = min(wave)
    wave_hi = max(wave)
    df[colname_flam] = fluxFactor*df[colname_flam]
    flam = df[colname_flam].tolist()
    data = {'wavelength': wave, 'flux': flam}

    f = interpolate.interp1d(data['wavelength'], data['flux'],
                             bounds_error=True,
                             kind='linear')

    return f,wave_lo,wave_hi


# Convert from redshift to luminosity distance.
#  Default values of H0 and Omega0 are from Bennett et al. (2014).
def zToDlum(z, H0=69.6, Om0=0.286):
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    # comoving distance...
    Dcom = cosmo.comoving_distance(z)
    # luminosity distance...
    Dlum = (1.+z)*Dcom
    return Dlum


# Convert Mpc_to_cm.
def Mpc_to_cm(Dmpc):
    Dcm = Dmpc*1.00e6*3.086e+18
    return Dcm