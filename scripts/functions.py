#### SMALL FUNCTIONS ####

## This function can be use with the python map() function.
# It applies the function bit2power to a bit and checks if a certain bit is checked.
# It simply returns 1 if the bit is checked.
def bit2power_map(x,bit):
    bitlist = bit2power(x) 
    if bit in bitlist:
        return(1)
    else:
        return(0)

## Create a mask with 0 where bit is not set and 1 where the bit is set.
# It also takes a bit list and combines them with AND.
# The output is a matrix with the same shape as the input mask and values 0 or 1.
def mkmask(mask,bitlist=[5]):
    
    mask_out = mask.copy() * 0
    for bit in bitlist:
        a = mask.reshape((1,mask.shape[0]*mask.shape[1]))[0]
        a_unique = np.unique(a)
        r = list(map(bit2power_map , a_unique, itertools.repeat(bit)))
        sel_has_bit = np.where( np.asarray(r) == 1)[0]
        a2 = a.copy() * 0
        for ii in sel_has_bit:
            a2[a == a_unique[ii]] = 1
        mask_out += a2.reshape(mask.shape)
    mask_out[mask_out >= 1] = 1
    return(mask_out)



def bit2power(n):
    sel = np.where( np.asarray(list(bin(n)[:1:-1])).astype("int") == 1 )[0]
    if len(sel) == 0:
        return([0])
    else:
        return(list(sel))
    


## Function to load Sergio's astrometry correction table
def load_astrometry_correction_table(file):
    with open(file, 'r') as fp:
        data = Table(names=["delta_ra_mas","delta_dec_mas"],dtype=["f","f"])
        for cnt, line in enumerate(fp,start=0):
            if cnt > 0:
                data.add_row([ float(line.split()[0]) , float(line.split()[1]) ])
                
    return(data)



## computes a pseudo chi2 measuring the deviation of a residual image from pure noise.
# input the histograms for the residual image and simulated noise from the plt.hist() function
def compute_deviation_noise(hist_res,hist_noise):
    x1 = hist_res[1]
    x1 = x1 - np.abs(x1[1]-x1[0])
    x1 = x1[1:]
    y1 = hist_res[0]

    x2 = hist_noise[1]
    x2 = x2 - np.abs(x2[1]-x2[0])
    x2 = x2[1:]
    y2 = hist_noise[0]
    y2_interp = np.interp(x1,x2,y2)
    
    chi2 = np.nansum( (y1 - y2_interp)**2 ) / len(x1)
    
    return(chi2)
    




# Returns modeled PSFex PSF for a 1D polynomial (e.g., magnitude)
# input loaded psfcube (psfcube) and its header (psfcube_h) together with parameter (e.g., magnitude)
def modelPSF_1D(psfcube,psfcube_h,param):
    X = psfcube[0][0][0].copy() # constant
    if psfcube_h["POLNAXIS"] > 0:
        x_scaled = (param - psfcube_h["POLZERO1"]) / psfcube_h["POLSCAL1"]
        for iii in range(psfcube_h["POLDEG1"]):
            X += psfcube[0][0][iii+1].copy() * x_scaled**(iii+1)
    return(X)


def is_number(n):
    try:
        float(n)   # Type-casting the string to `float`.
                   # If string is not a valid `float`, 
                   # it'll raise `ValueError` exception
    except ValueError:
        return False
    return True

def clip(dat,n,niter):
    dat = dat.ravel()
    
    for jj in range(niter):
        med = np.nanmedian(dat)
        stdev = np.nanstd(dat)
    
        index = np.where( (dat < (med+n*stdev)) & (dat > (med-n*stdev)) )[0]
        
        out = {"med": med,
              "stdev": stdev}

        if len(index) == len(dat):
            return(out)
        
        if len(index) > 0:
            dat = dat[index]
        
    return(out)

def replace_in_file(filename, old_string, new_string):
    # Safely read the input filename using 'with'
    with open(filename) as f:
        s = f.read()
        if old_string not in s:
            print("Could not find " + old_string)
            return

    # Safely write the changed content, if found in the file
    with open(filename, 'w') as f:
        s = s.replace(old_string, new_string)
        f.write(s)

        
        
def get_default_parfile():
    
    PARS_default = {"CATALOG_NAME": "output.cat", # UPDATE
        "CATALOG_TYPE": "ASCII_HEAD",
       "PARAMETERS_NAME": "./input/sex.par", # UPDATE
        "DETECT_TYPE": "CCD",
        "DETECT_MINAREA": 5, #5
        "THRESH_TYPE": "RELATIVE",
        "DETECT_THRESH": 3,
        "ANALYSIS_THRESH": 1.5,
        "FILTER": "Y",
        "FILTER_NAME": "./input/g2.8.conv", # UPDATE
        "FILTER_THRESH": "",
        "DEBLEND_NTHRESH": 32,
        "DEBLEND_MINCONT": 0.01,
        "CLEAN": "Y",
        "CLEAN_PARAM": 1.0,
        "MASK_TYPE": "CORRECT",
        "WEIGHT_TYPE": "NONE", # UPDATE (this is in command line)
        "WEIGHT_IMAGE": "weight.fits", # UPDATE (this is in command line)
        "WEIGHT_GAIN": "Y",
        "WEIGHT_THRESH": "",
        "FLAG_IMAGE":"NONE",
        "FLAG_TYPE":"OR",
        "PHOT_APERTURES": "30", # UPDATE
        "PHOT_AUTOPARAMS": "2.5,3.5",
        "PHOT_PETROPARAMS": "2.0,3.5",
        "PHOT_AUTOAPERS": "0.0,0.0",
        "PHOT_FLUXFRAC": 0.5,
        "SATUR_LEVEL": 50000.0, # UPDATE
        "MAG_ZEROPOINT": 0.0, # UPDATE
        "MAG_GAMMA": 4.0,
        "GAIN": 0.0, # UPDATE
        "PIXEL_SCALE": 0, # UPDATE
        "SEEING_FWHM": 0.1, # UPDATE
        "STARNNW_NAME": "./input/default.nnw", # UPDATE
        "BACK_TYPE": "AUTO",
        "BACK_VALUE": 0.0,
        "BACK_SIZE": 64,
        "BACK_FILTERSIZE": 3,
        "BACKPHOTO_TYPE": "GLOBAL",
        "BACKPHOTO_THICK": 24,
        "BACK_FILTTHRESH": 0.0,
        "CHECKIMAGE_TYPE": "APERTURES,SEGMENTATION,BACKGROUND", # UPDATE
        "CHECKIMAGE_NAME": "aper.fits, seg.fits, back.fits", # UPDATE
        "NTHREADS":1 # UPDATE?
       }
    
    return(PARS_default)


def get_apertures(pixscale,aperturelist_arcsec=[0.5 , 2.0 , 3.0]):
    
    # apertures in arcsec
    apertures_arcsec = np.asarray(aperturelist_arcsec)
    
    # aperture in pixels (using pixel scale)
    apertures_pix = apertures_arcsec / pixscale
    
    return(apertures_pix)


    
def ConvertToMicroJansky(flux,zp):
    
    if zp == -99:
        tmp3 = flux
    else: 
        tmp1 = -2.5*np.log10(flux) + zp # mag
        tmp2 = 10.0**(-0.4*(tmp1 + 48.6) ) # erg/s/cm2/Hz
        tmp3 = tmp2 / 1e-23 * 1e6 # micro-Jy
    
    return(tmp3)

def print_time_summary(dat):
    
    print("\n\n++++++ TIME SUMMARY +++++++")
    for key in dat.keys():
        print("%s: %g seconds" % (dat[key][0],dat[key][1]) )
    print("Total time: %g seconds" % (np.sum(np.asarray( [dat[key][1] for key in dat.keys()] ) )) )
    print("++++++++++++++++++++++++++\n\n")
    
def save_time_summary(dat,file):
    with open(file,"w") as f:
        for key in dat.keys():
            f.write( "%s: %g seconds \n" % (dat[key][0],dat[key][1]) )
        f.write( "Total time: %g seconds" % (np.sum(np.asarray( [dat[key][1] for key in dat.keys()] ) )) )

def save_stats(dat,file):
    with open(file,"w") as f:
        for ii in range(len(dat)):
            f.write( "%s\n" % (dat[ii]) )



## This function saves a FITS file with original image, model, residual, segmentation map from the HR fitting.
def save_hr_fits(img,mod,res,seg,header,outfile):
    
    h1 = header.copy()
    h2 = header.copy()
    try:
        np.isfinite(h2["SIMPLE"])
    except Exception as e:
        h2.rename_keyword("SIMPLE","XTENSION")
    h2["PCOUNT"] = 0
    h2["GCOUNT"] = 1
    h2["XTENSION"] = "IMAGE"

    hdu_1 = fits.PrimaryHDU()
    hdu_1.data = img.copy()
    hdu_1.header = h1.copy()
    hdu_1.header["EXTNAME"] = "ORIGINAL"

    hdu_2 = fits.ImageHDU()
    hdu_2.data = mod.copy()
    hdu_2.header = h2.copy()
    hdu_2.header["EXTNAME"] = "COMPL_MODEL"

    hdu_3 = fits.ImageHDU()
    hdu_3.data = res.copy()
    hdu_3.header = h2.copy()
    hdu_3.header["EXTNAME"] = "COMPL_RES"

    hdu_4 = fits.ImageHDU()
    hdu_4.data = seg.copy() # already loaded in the begining.
    hdu_4.header = h2.copy()
    hdu_4.header["EXTNAME"] = "SEGMENTATION"

    hdul_new = fits.HDUList([hdu_1,hdu_2,hdu_3,hdu_4]) # original, compl:[model, residual], segmentation
    hdul_new.verify("silentfix") # fix everything that is broken
    hdul_new.writeto(outfile, overwrite=True)
    
    return(True)
    

## This function saves a FITS file with original image, model, residual, chi2 map from the LR fitting
def save_lr_fits(img,mod,res,seg,seg_sextractor,mask_bad,mask,header,mask_header,outfile):

    h1 = header.copy()
    h2 = header.copy()
    h3 = mask_header.copy()
    try:
        np.isfinite(h2["SIMPLE"])
    except Exception as e:
        h2.rename_keyword("SIMPLE","XTENSION")
    h2["PCOUNT"] = 0
    h2["GCOUNT"] = 1
    h2["XTENSION"] = "IMAGE"

    hdu_0 = fits.PrimaryHDU()
    hdu_0.data = img.copy()
    hdu_0.header = h1.copy()
    hdu_0.header["EXTNAME"] = "ORIGINAL"

    hdu_1 = fits.ImageHDU()
    hdu_1.data = mod.copy()
    hdu_1.header = h2.copy()
    hdu_1.header["EXTNAME"] = "COMPL_MODEL"

    hdu_2 = fits.ImageHDU()
    hdu_2.data = res.copy()
    hdu_2.header = h2.copy()
    hdu_2.header["EXTNAME"] = "COMPL_RES"
    
    hdu_3 = fits.ImageHDU()
    hdu_3.data = seg.copy()
    hdu_3.header = h2.copy()
    hdu_3.header["EXTNAME"] = "SEGMAP"
    
    hdu_4 = fits.ImageHDU()
    hdu_4.data = seg_sextractor.copy()
    hdu_4.header = h2.copy()
    hdu_4.header["EXTNAME"] = "SEGMAP_SEXTRACTOR"

    hdu_5 = fits.ImageHDU()
    hdu_5.data = mask.copy()
    hdu_5.header = h3.copy()
    hdu_5.header["EXTNAME"] = "MASK"

    hdu_6 = fits.ImageHDU()
    hdu_6.data = mask_bad.copy()
    hdu_6.header = h2.copy()
    hdu_6.header["EXTNAME"] = "MASK_BAD"
    

    hdul_new = fits.HDUList([hdu_0,hdu_1,hdu_2,hdu_3,hdu_4,hdu_5,hdu_6]) # original, compl:[model, residual, ], segmap, segmap_sextractor , mask, mask_bad
    hdul_new.verify("silentfix") # fix everything that is broken
    hdul_new.writeto(outfile, overwrite=True)
    
    return(True)
    
    
    
##################


#### MAIN FUNCTIONS #####

