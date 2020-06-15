import os, sys

import numpy as np

import subprocess
import random
import time
import multiprocessing as multiproc

import itertools

import json

import sh

from astropy.io import fits, ascii
from astropy.table import Table, Column, MaskedColumn, hstack, vstack
import astropy.wcs as wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.nddata import Cutout2D

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects as path_effects

from scipy import stats
from scipy import ndimage

import tractor
from tractor import *
from tractor.galaxy import *
from tractor.sersic import *

### FORCED PHOTOMETRY FITTING OF LR IMAGE WITH APPLIED POSITION OFFSET AND OBJECT BY OBJECT (new) #####

## Fit LR image using the HR model
# only works for one tile ID
# seg_map_type can be "sextractor" or "hrmodel"
def fit_LR_image(userinput , tileid , seg_map_type): 
    '''
    This function creates uses the model created on the high-resolution image to perform
    forced photometry on the low-resolution image.

    USAGE: get_LR_model(userinput , tileids , seg_map_type)
    where
        - userinput: is a user input dictionary
        - tileids: the IDs of the tiles on which the model should be created
        - seg_map_type: either "sextractor" or "hrmodel"
    '''



    ## 1. Create process ID ===========================
    this_tile_id = tileid
    #this_tile_id = 3

    TIMESUMMARY = dict()
    STATS = []

    process_id = "%04.0f_%s" % (this_tile_id , userinput["lr_large_image_name"].split("/")[-1].split(".fits")[0] ) # set the process ID to the name of the large image array + the tile number
    print("This Process is: " + str(process_id))
    STATS.append("This Process is: " + str(process_id))


    ## 2. create (temporary) working directory ===========================
    dir_this_process = os.path.join(userinput["workdir"], "%s_%s" % (userinput["out_prefix"],process_id) )
    if not os.path.exists(dir_this_process):
        sh.mkdir('-p', dir_this_process)
    else:
        print("Directory %s already exists" % dir_this_process)


    ## 3. Load images ===========
    print("reading images . . .",end="")
    start_time = time.time()

    ## low resolution segmentation map and extent file
    ## choose here either segmap from HR image model or from SExtractor
    if seg_map_type == "sextractor":
        lr_segmap_sextractor_path = os.path.join(dir_this_process,"lr_seg.fits")
        lr_segmap_path = os.path.join(dir_this_process,"lr_segmap_from_sextractor.fits")
        lr_segmap_extentfile_path = os.path.join(dir_this_process,"lr_segmap_extensions_sextractor.csv")
    elif seg_map_type == "hrmodel":
        lr_segmap_sextractor_path = os.path.join(dir_this_process,"lr_seg.fits") # still need SExtractor segmentation map for flagging in the end.
        lr_segmap_path = os.path.join(dir_this_process,"lr_segmap_from_hr.fits")
        lr_segmap_extentfile_path = os.path.join(dir_this_process,"lr_segmap_extensions.csv")
    else:
        print("Segmentation map type is not known")
        #return(False)

    with fits.open(lr_segmap_sextractor_path) as hdul:
        lr_segmap_sextractor = hdul[0].data
    with fits.open(lr_segmap_path) as hdul:
        lr_segmap = hdul[0].data
    lr_segmap_extent_table = ascii.read(lr_segmap_extentfile_path , format="csv")



    ## low resolution
    lr_image_path = os.path.join(userinput["cutout_output_dir"],userinput["lr_large_image_name"].split("/")[-1].split(".fits")[0],"%04.0f_%s.fits" % (this_tile_id,userinput["lr_large_image_name"].split("/")[-1].split(".fits")[0]))
    with fits.open(lr_image_path) as hdul:
        lr_img = hdul[0].data
        lr_img_h = hdul[0].header
        lr_mask = hdul[1].data
        lr_mask_h = hdul[1].header

    lr_pixscale = np.abs(lr_img_h["CD1_1"]*3600) # pixel scale in arcsec/px
    lr_img_wcs = wcs.WCS(lr_img_h)


    ## high resolution (need this for pixel scale)
    hr_image_path = os.path.join(userinput["cutout_output_dir"],userinput["lr_large_image_name"].split("/")[-1].split(".fits")[0],"%04.0f_%s.fits" % (this_tile_id,userinput["hr_large_image_name"].split("/")[-1].split(".fits")[0]))
    with fits.open(hr_image_path) as hdul:
        hr_img_h = hdul[0].header

    hr_pixscale = np.abs(hr_img_h["CD1_1"]*3600) # pixel scale in arcsec/px
    hr_img_wcs = wcs.WCS(hr_img_h)


    print(" done (in %g seconds)" % (round((time.time()-start_time),2)) )
    TIMESUMMARY[0] = ["Reading images",round((time.time()-start_time),2)]


    ## 3.5 Get LR image properties ============

    ## Check some stuff
    if "CD1_2" not in lr_img_h.keys():
        lr_img_h["CD1_2"] = 0
    if "CD2_1" not in lr_img_h.keys():
        lr_img_h["CD2_1"] = 0

    ## Get Pixel noise
    tmp = clip(lr_img[lr_segmap == 0], n=3, niter=10)
    lr_img_pixnoise = tmp["stdev"]
    lr_img_medbkg = tmp["med"]
    print("Pixelnoise of low-res image: %g" % lr_img_pixnoise )
    print("Median background of low-res image: %g" % lr_img_medbkg )
    STATS.append("Pixelnoise of low-res image: %g" % lr_img_pixnoise )
    STATS.append("Median background of low-res image: %g" % lr_img_medbkg)


    ## Get Astrometry offset in RA and DEC in pixels
    # Be careful here with direction!!!!
    astro_offset = load_astrometry_correction_table(userinput["lr_astrometry_correction_name"])
    hsc_gaia_median_ra = np.nanmedian(astro_offset["delta_ra_mas"]) # in mas
    hsc_gaia_median_dec = np.nanmedian(astro_offset["delta_dec_mas"]) # in mas

    hr_astro_offset = load_astrometry_correction_table(userinput["hr_astrometry_correction_name"])
    acs_gaia_median_ra = np.nanmedian(hr_astro_offset["delta_ra_mas"]) # in mas
    acs_gaia_median_dec = np.nanmedian(hr_astro_offset["delta_dec_mas"]) # in mas
    #acs_gaia_median_ra = 15.8
    #acs_gaia_median_dec = 2.3

    hsc_acs_median_ra = hsc_gaia_median_ra - acs_gaia_median_ra # in mas
    hsc_acs_median_dec = hsc_gaia_median_dec - acs_gaia_median_dec # in mas

    delta_X_median_px = (-1)*hsc_acs_median_dec/1000 / lr_pixscale
    delta_Y_median_px = hsc_acs_median_ra/1000 / lr_pixscale

    #delta_X_median_px = (-1)*np.nanmedian(astro_offset["delta_dec_mas"])/1000 / lr_pixscale
    #delta_Y_median_px = np.nanmedian(astro_offset["delta_ra_mas"])/1000 / lr_pixscale
    print("Astrometric offset in mas of LR image (ra,dec) = (%5.3g,%5.3g)" % (hsc_gaia_median_ra,hsc_gaia_median_dec))
    print("Astrometric offset in mas of HR image (ra,dec) = (%5.3g,%5.3g)" % (acs_gaia_median_ra,acs_gaia_median_dec))
    print("Astrometric offset in mas of image (ra,dec) = (%5.3g,%5.3g)" % (hsc_acs_median_ra,hsc_acs_median_dec))
    print("Astrometric offset in pixels (X,Y) = (%5.3g,%5.3g)" % (delta_X_median_px,delta_Y_median_px))

    STATS.append("Astrometric offset in mas of LR image (ra,dec) = (%5.3g,%5.3g)" % (hsc_gaia_median_ra,hsc_gaia_median_dec))
    STATS.append("Astrometric offset in mas of HR image (ra,dec) = (%5.3g,%5.3g)" % (acs_gaia_median_ra,acs_gaia_median_dec))
    STATS.append("Astrometric offset in mas (ra,dec) = (%5.3g,%5.3g)" % (hsc_acs_median_ra,hsc_acs_median_dec))
    STATS.append("Astrometric offset in pixels (X,Y) = (%5.3g,%5.3g)" % (delta_X_median_px,delta_Y_median_px))


    ## 4. Read High-resolution table with model fits ===========
    print("Reading high-resolution model table . . .",end="")
    start_time = time.time()

    outname = os.path.join(dir_this_process,"hr_compl_table_final.fits")
    if not os.path.exists(outname):
        print("The models for this high-resolution image (%s, tile %g) do not exist, yet. Please run get_HR_model() first." % (userinput["lr_large_image_name"],this_tile_id) )
        #continue # UNCOMMENT THIS ONCE EVERYTHING IS IN A LOOP
    cat_models_hr = Table.read(outname)
    cat_models_hr.sort("NUMBER.hr")

    print(" done (in %g seconds)" % (round((time.time()-start_time),2)) )
    TIMESUMMARY[len(TIMESUMMARY)+1] = ["Reading high-resolution model table",round((time.time()-start_time),2)]



    ## 5. Load PSF (LR image) ==============
    if is_number(userinput["lr_image_psf"]):
        PSF_FWHM_ARCSEC_LR = userinput["lr_image_psf"]
        PSF_LR_TYPE = "fwhm"
        print("Using given FWHM and gaussian for PSF of low-res image.")
        STATS.append("Using given FWHM and gaussian for PSF of low-res image.")
    else:
        PSF_LR_TYPE = "pixel"
        with fits.open(os.path.join(userinput["lr_image_psf"])) as hdul:
            PSF_LR_PIXEL = modelPSF_1D(psfcube=hdul[1].data,psfcube_h=hdul[1].header,param=21) #22
            #tmp = Cutout2D(data=PSF_LR_PIXEL,position=ndimage.measurements.maximum_position(PSF_LR_PIXEL),size=(21,21) , copy=True , mode="trim")
            #PSF_LR_PIXEL = tmp.data.copy()
            PSF_LR_PIXEL = PSF_LR_PIXEL / np.nansum(PSF_LR_PIXEL)
        print("Using pixel PSF for low-res image.")
        STATS.append("Using pixel PSF for low-res image.")


    ## 6. Select galaxies ===========

    # add local x/y to the table
    cat_models_hr["X_IMAGE_hr.lr"] = np.zeros(len(cat_models_hr))
    cat_models_hr["Y_IMAGE_hr.lr"] = np.zeros(len(cat_models_hr))
    cat_models_hr["X_IMAGE_hr_withoffset.lr"] = np.zeros(len(cat_models_hr))
    cat_models_hr["Y_IMAGE_hr_withoffset.lr"] = np.zeros(len(cat_models_hr))
    for ii in range(len(cat_models_hr)):
        tmp = lr_img_wcs.all_world2pix([[cat_models_hr["RA_tractor.hr"][ii],cat_models_hr["DEC_tractor.hr"][ii]]],1)
        cat_models_hr["X_IMAGE_hr.lr"][ii] = tmp[0][0]
        cat_models_hr["Y_IMAGE_hr.lr"][ii] = tmp[0][1]
        cat_models_hr["X_IMAGE_hr_withoffset.lr"][ii] = tmp[0][0] + delta_X_median_px
        cat_models_hr["Y_IMAGE_hr_withoffset.lr"][ii] = tmp[0][1] + delta_Y_median_px

    # select galaxies on image
    sel_good = np.where( (cat_models_hr["X_IMAGE_hr.lr"] >= 1)
                        & (cat_models_hr["X_IMAGE_hr.lr"] <= lr_img.shape[1])
                        & (cat_models_hr["Y_IMAGE_hr.lr"] >= 1)
                        & (cat_models_hr["Y_IMAGE_hr.lr"] <= lr_img.shape[0])
                        )[0]
    cat_models_hr_use = cat_models_hr[sel_good].copy()

    print("%g galaxies on the whole image" % len(cat_models_hr_use))
    STATS.append("%g galaxies on the whole image" % len(cat_models_hr_use))



    ## 7. Create object lists ==============
    OBJECTS_TO_FIT = Table(data=np.asarray(cat_models_hr_use["NUMBER.hr"]).reshape((len(cat_models_hr_use),1)),names=["NUMBER.hr"]) # store here the NUMBERs of the objects that still have to be fitted
    OBJECTS_FITTED = [] # store there NUMBERs of objects that have been fitted.
    OBJECTS_FAILED = [] # store there NUMBERs of objects for which the fit failed.


    ## 8. Run Tractor object by object =========
    print("Run Tractor object by object . . .")
    start_time = time.time()

    TABLE = Table()

    group_counter = 0
    while len(OBJECTS_TO_FIT) > 0:



        ## 8.1 get object number and see where it is in the *_use model catalog ####
        main_obj_number = OBJECTS_TO_FIT["NUMBER.hr"][0] # always take next galaxy in line
        print("\nMain Object: %g (%g objects left)" % (main_obj_number , len(OBJECTS_TO_FIT)) )

        id_in_hr_cat = np.where( cat_models_hr_use["NUMBER.hr"] == main_obj_number )[0][0]

        ## check in which segmap area it falls.
        seg_area_this_object = lr_segmap[int(cat_models_hr_use["Y_IMAGE_hr.lr"][id_in_hr_cat]-1),int(cat_models_hr_use["X_IMAGE_hr.lr"][id_in_hr_cat]-1)] # remember to flip X and Y
        print("Object belongs to segmap area %g" % seg_area_this_object)

        ## create box in which to choose the objects and make cutout
        obj_x = cat_models_hr_use["X_IMAGE_hr.lr"][id_in_hr_cat]-1
        obj_y = cat_models_hr_use["Y_IMAGE_hr.lr"][id_in_hr_cat]-1

        if seg_area_this_object != 0: # in case segmentation area exists
            xmin = float(lr_segmap_extent_table[lr_segmap_extent_table["label"] == seg_area_this_object]["XMIN"])
            xmax = float(lr_segmap_extent_table[lr_segmap_extent_table["label"] == seg_area_this_object]["XMAX"])
            ymin = float(lr_segmap_extent_table[lr_segmap_extent_table["label"] == seg_area_this_object]["YMIN"])
            ymax = float(lr_segmap_extent_table[lr_segmap_extent_table["label"] == seg_area_this_object]["YMAX"])

            xmin = xmin - (xmax-xmin)*0.1
            xmax = xmax + (xmax-xmin)*0.1
            ymin = ymin - (ymax-ymin)*0.1
            ymax = ymax + (ymax-ymin)*0.1

        else: # if no segmap area, cut a region 4" x 4"
            xmin = obj_x - 2 / lr_pixscale # 2"
            xmax = obj_x + 2 / lr_pixscale # 2"
            ymin = obj_y - 2 / lr_pixscale # 2"
            ymax = obj_y + 2 / lr_pixscale # 2"


        ## 8.2 cut out #####
        tmp = Cutout2D(data=lr_img,position=( xmin+(xmax-xmin)/2-1,ymin+(ymax-ymin)/2-1),size=(ymax-ymin , xmax-xmin ) , copy=True , mode="trim",wcs=lr_img_wcs)
        img_cutout = tmp.data.copy()
        img_cutout_wcs = tmp.wcs.copy()

        tmp = Cutout2D(data=lr_segmap,position=( xmin+(xmax-xmin)/2-1,ymin+(ymax-ymin)/2-1),size=(ymax-ymin , xmax-xmin ) , copy=True , mode="trim",wcs=lr_img_wcs)
        lr_segmap_cutout = tmp.data.copy()


        ## 8.3 Get all objects that are in the cutout, get their new X/Y values, and flag those that do not belong to the same segmap area ####
        sel_objects_in_cutout = np.where( (cat_models_hr_use["X_IMAGE_hr.lr"]-1 >= xmin)
                        & (cat_models_hr_use["X_IMAGE_hr.lr"]-1 < xmax)
                        & (cat_models_hr_use["Y_IMAGE_hr.lr"]-1 >= ymin)
                        & (cat_models_hr_use["Y_IMAGE_hr.lr"]-1 < ymax)
                        )[0]
        cat_models_hr_use_incutout = cat_models_hr_use[sel_objects_in_cutout].copy()
        id_objects_in_cutout = cat_models_hr_use_incutout["NUMBER.hr"]
        print("Number of objects in cutout: %g" % len(sel_objects_in_cutout))

        # get X/Y of these objects
        cat_models_hr_use_incutout["X_IMAGE_hr_cutout.lr"] = np.zeros(len(id_objects_in_cutout))
        cat_models_hr_use_incutout["Y_IMAGE_hr_cutout.lr"] = np.zeros(len(id_objects_in_cutout))
        for ii in range(len(sel_objects_in_cutout)):
            tmp = img_cutout_wcs.all_world2pix([[cat_models_hr_use_incutout["RA_tractor.hr"][ii],cat_models_hr_use_incutout["DEC_tractor.hr"][ii]]],1)
            cat_models_hr_use_incutout["X_IMAGE_hr_cutout.lr"][ii] = tmp[0][0]
            cat_models_hr_use_incutout["Y_IMAGE_hr_cutout.lr"][ii] = tmp[0][1]


        # flag the ones that do not belong to the same segmap area   
        # However, if current segmap area does not exist (i.e., 0), then refit all the other galaxies automatically
        cat_models_hr_use_incutout["in_same_segmap_area"] = (np.zeros(len(id_objects_in_cutout))).astype("int") # 1 if in same area, start all with 0
        if seg_area_this_object == 0:
            cat_models_hr_use_incutout["in_same_segmap_area"][cat_models_hr_use_incutout["NUMBER.hr"] == main_obj_number] = 1 # means: fit all galaxies in box, but only keep main object
        else:
            for iii in range(len(sel_objects_in_cutout)):
                i_this = int(cat_models_hr_use_incutout["Y_IMAGE_hr_cutout.lr"][iii]-1)
                j_this = int(cat_models_hr_use_incutout["X_IMAGE_hr_cutout.lr"][iii]-1)
                if (i_this >= 0) & (i_this < lr_segmap_cutout.shape[0]) & (j_this >= 0) & (j_this < lr_segmap_cutout.shape[1]):
                    pix_val_this = lr_segmap_cutout[i_this,j_this] # remember to flip X and Y
                    if pix_val_this == seg_area_this_object:
                        cat_models_hr_use_incutout["in_same_segmap_area"][iii] = int(1)
                else:
                    if len(sel_objects_in_cutout) == 1:
                        cat_models_hr_use_incutout["in_same_segmap_area"][iii] = int(1)
                    else:
                        cat_models_hr_use_incutout["in_same_segmap_area"][iii] = int(0)


        ########

        ## 8.4 Create Tractor image  ###########
        # Create a complex image with pixelized PSF as well as a simple image with simple gaussian PSF

        # create Tractor image
        if PSF_LR_TYPE == "fwhm":
            tim_compl_lr = Image(data=img_cutout, invvar=np.ones_like(img_cutout) / (lr_img_pixnoise**2),
                        psf=NCircularGaussianPSF([userinput["lr_image_psf"]/2.35/lr_pixscale], [1.]),
                        wcs=NullWCS(), photocal=NullPhotoCal(),
                        sky=ConstantSky(lr_img_medbkg))
            tim_simple_lr = Image(data=img_cutout, invvar=np.ones_like(img_cutout) / (lr_img_pixnoise**2),
                        psf=NCircularGaussianPSF([0.7/2.35/lr_pixscale], [1.]),
                        wcs=NullWCS(), photocal=NullPhotoCal(),
                        sky=ConstantSky(lr_img_medbkg))
        if PSF_LR_TYPE == "pixel":
            tim_compl_lr = Image(data=img_cutout, invvar=np.ones_like(img_cutout) / (lr_img_pixnoise**2),
                        psf=PixelizedPSF(PSF_LR_PIXEL),
                        wcs=NullWCS(), photocal=NullPhotoCal(),
                        sky=ConstantSky(lr_img_medbkg))
            tim_simple_lr = Image(data=img_cutout, invvar=np.ones_like(img_cutout) / (lr_img_pixnoise**2),
                        psf=NCircularGaussianPSF([0.7/2.35/lr_pixscale], [1.]),
                        wcs=NullWCS(), photocal=NullPhotoCal(),
                        sky=ConstantSky(lr_img_medbkg))


        ## 8.5 Run TRACTOR in simple mode to get position offset (in case astrometry is not perfect) #######

        # create source list (create point sources at the position of the HR models)
        # This can be change to non-pointsource fits if needed.
        src_simple_lr = []
        for ii in range(len(sel_objects_in_cutout)):    
            tmp = img_cutout_wcs.all_world2pix([[cat_models_hr_use_incutout["RA_tractor.hr"][ii],cat_models_hr_use_incutout["DEC_tractor.hr"][ii]]],0)
            #pos_init = PixPos(tmp[0][0],tmp[0][1])
            pos_init = PixPos(tmp[0][0]+delta_X_median_px , tmp[0][1]+delta_Y_median_px)
            pos_init.lowers = [pos_init[0]-0.5,pos_init[1]-0.5] # vary in 1 pixel
            pos_init.uppers = [pos_init[0]+0.5,pos_init[1]+0.5] # vary in 1 pixel

            flux_init = Flux(cat_models_hr_use_incutout["brightness.Flux.hr"][ii] * 10**(0.4*( userinput["lr_zp"] - userinput["hr_zp"] ))) 

            if cat_models_hr_use_incutout["is_pointsource"][ii] == 2: # if point source
                thissource = PointSource(pos_init,flux_init)
            else: # if extended 
                n_init = SersicIndex(cat_models_hr_use_incutout["sersicindex.SersicIndex.hr"][ii])
                ba_init = cat_models_hr_use_incutout["shape.ab.hr"][ii]
                re_init = cat_models_hr_use_incutout["shape.re.hr"][ii] * (hr_pixscale / lr_pixscale)
                phi_init = cat_models_hr_use_incutout["shape.phi.hr"][ii]
                thissource = SersicGalaxy(pos_init,
                                            flux_init,
                                            GalaxyShape(re_init, ba_init, phi_init), n_init)

            src_simple_lr.append(thissource)

        print("    Number of galaxies for fitting (simple models): " + str(len(src_simple_lr)))

        # create tractor object (use simple image - not pixelized PSF for faster run)
        tractor_simple_lr = Tractor([tim_compl_lr] , src_simple_lr , optimizer=ConstrainedOptimizer())

        # freeze image calibration params
        tractor_simple_lr.freezeParam('images')

        # save initial positions
        posx0 = [src_simple_lr[iii].pos[0] for iii in range(len(src_simple_lr))]
        posy0 = [src_simple_lr[iii].pos[1] for iii in range(len(src_simple_lr))]

        # Freese all parameters except position
        # None to freeze for pointsource fit.
        for ii in range(len(src_simple_lr)):
            if cat_models_hr_use_incutout["is_pointsource"][ii] == 2: # if point source
                pass
            else: # if extended
                src_simple_lr[ii].shape.freezeParams("phi")
                src_simple_lr[ii].shape.freezeParams("ab")
                src_simple_lr[ii].shape.freezeParams("re")
                src_simple_lr[ii].freezeParams("sersicindex")
                #src_simple_lr[ii].freezeParams("pos")
                pass


        # run it and extract models and chi2
        print("    Forced photometry pointsource fitting to LR image to get position offset ",end="")
        start_time2 = time.time()
        for ii in range(userinput["lr_compl_niter"]):
            dlnp,X,alpha,variances_simple = tractor_simple_lr.optimize(variance=True)
            print(str(ii)+" ("+str(dlnp)+") ",end=" ")
            if dlnp < 1e-3:
                break
        print(" done (in %g seconds)" % (round((time.time()-start_time2),2)) )

        ## Create table with variances for position (NEW)
        # DO THIS BEFORE THAWNING PARAMETERS!
        # Note: Variances array comes in packages for each parameter (e.g., variances = x,x,x , y,y,y , z,z,z). These 
        # packages might not have same length as the parameters fitted for each sources may vary.
        # Need to use pop() to remove number of parameters.
        # Here it's much easier because only flux is fitted. The length of the variance array should therefore
        # be the number of sources.
        variances_simple = list(variances_simple) # make list
        variance_table_position = Table(data=np.zeros((len(src_simple_lr),3))-99, names=["pos.x.var","pos.y.var","brightness.Flux.var"])
        for sss,src in enumerate(src_simple_lr):
            params_fitted = src.getParamNames()
            params_variances = [variances_simple.pop(0) for vvv in range(len(params_fitted))]
            for ppp,param in enumerate(params_fitted):
                variance_table_position["%s.var" % param][sss] = params_variances[ppp]
        variance_table_position = variance_table_position[["pos.x.var","pos.y.var"]].copy()

        # save final positions
        posx1 = [src_simple_lr[iii].pos[0] for iii in range(len(src_simple_lr))]
        posy1 = [src_simple_lr[iii].pos[1] for iii in range(len(src_simple_lr))]

        # compute offset between final and initial positions in pixels
        Dx = np.asarray(posx1) - np.asarray(posx0)
        Dy = np.asarray(posy1) - np.asarray(posy0)

        print("Median offset in x: %g px" % np.median(Dx))
        print("Median offset in y: %g px" % np.median(Dy))

        # get the model, chi2, and residual image (don't need this)      
        #lr_mod_simple = tractor_simple_lr.getModelImage(0)
        #lr_chi2_simple = tractor_simple_lr.getChiImage(0)
        #lr_residual_simple = img_cutout - lr_mod_simple



        ## 8.6 Run TRACTOR on complex models with applied position offset ############

        # create source list (fixed position from the pointsource fit)
        src_compl_lr = []
        for ii in range(len(sel_objects_in_cutout)):    
            #tmp = img_cutout_wcs.all_world2pix([[cat_models_hr_use_incutout["RA_tractor.hr"][ii],cat_models_hr_use_incutout["DEC_tractor.hr"][ii]]],0)
            #pos_init = PixPos(tmp[0][0]+delta_X_median_px+Dx[ii],tmp[0][1]+delta_Y_median_px+Dy[ii])
            pos_init = PixPos( src_simple_lr[ii].pos[0] , src_simple_lr[ii].pos[1] )
            flux_init = Flux(cat_models_hr_use_incutout["brightness.Flux.hr"][ii] * 10**(0.4*( userinput["lr_zp"] - userinput["hr_zp"] ))) 

            if cat_models_hr_use_incutout["is_pointsource"][ii] == 2: # if point source
                thissource = PointSource(pos_init,flux_init)
            else: # if extended 
                n_init = SersicIndex(cat_models_hr_use_incutout["sersicindex.SersicIndex.hr"][ii])
                ba_init = cat_models_hr_use_incutout["shape.ab.hr"][ii]
                re_init = cat_models_hr_use_incutout["shape.re.hr"][ii] * (hr_pixscale / lr_pixscale)
                phi_init = cat_models_hr_use_incutout["shape.phi.hr"][ii]
                thissource = SersicGalaxy(pos_init,
                                            flux_init,
                                            GalaxyShape(re_init, ba_init, phi_init), n_init)
            src_compl_lr.append(thissource)


        print("    Number of galaxies for fitting (complex models): " + str(len(src_compl_lr)))


        # create tractor object
        tractor_compl_lr = Tractor([tim_compl_lr] , src_compl_lr)

        # freeze image calibration params
        tractor_compl_lr.freezeParam('images')

        # Freese all parameters except brightness
        for ii in range(len(src_compl_lr)):

            if cat_models_hr_use_incutout["is_pointsource"][ii] == 2:
                src_compl_lr[ii].freezeParams("pos")
                pass
            else:
                src_compl_lr[ii].shape.freezeParams("phi")
                src_compl_lr[ii].shape.freezeParams("ab")
                src_compl_lr[ii].shape.freezeParams("re")
                src_compl_lr[ii].freezeParams("sersicindex")
                src_compl_lr[ii].freezeParams("pos")
                pass


        # run it and extract models and chi2
        print("    Forced photometry fitting to LR image ",end="")
        start_time2 = time.time()
        for ii in range(userinput["lr_compl_niter"]):
            dlnp,X,alpha,variances = tractor_compl_lr.optimize(variance=True)
            print(str(ii)+" ("+str(dlnp)+") ",end=" ")
            if dlnp < 1e-3:
                break
        print(" done (in %g seconds)" % (round((time.time()-start_time2),2)) )

        # get the model, chi2, and residual image (comment after testing, will create final image later)   
        #lr_mod_compl = tractor_compl_lr.getModelImage(0)
        #lr_chi2_compl = tractor_compl_lr.getChiImage(0)
        #lr_residual_compl = img_cutout - lr_mod_compl

        ## Create table with variances (NEW)
        # DO THIS BEFORE THAWNING PARAMETERS!
        # Note: Variances array comes in packages for each parameter (e.g., variances = x,x,x , y,y,y , z,z,z). These 
        # packages might not have same length as the parameters fitted for each sources may vary.
        # Need to use pop() to remove number of parameters.
        # Here it's much easier because only flux is fitted. The length of the variance array should therefore
        # be the number of sources.
        variances = list(variances) # make list
        variance_table = Table(data=np.zeros((len(src_compl_lr),1))-99, names=["brightness.Flux.var"])
        for sss,src in enumerate(src_compl_lr):
            params_fitted = src.getParamNames()
            params_variances = [variances.pop(0) for vvv in range(len(params_fitted))]
            for ppp,param in enumerate(params_fitted):
                variance_table["%s.var" % param][sss] = params_variances[ppp]

        variance_table = hstack([variance_table_position , variance_table])


        ## 8.7. Create Table ##############
        for ii in range(len(src_compl_lr)):
            src_compl_lr[ii].thawAllRecursive()

        param_keys = ['pos.x',
                        'pos.y',
                        'brightness.Flux',
                        'shape.re','shape.ab',
                        'shape.phi',
                        'sersicindex.SersicIndex']
        param_keys = param_keys + ["pos.x.var","pos.y.var","brightness.Flux.var"]
        param_keys = [key + ".lr" for key in param_keys]
        tab_this_cutout = Table(data=np.zeros((len(src_compl_lr),len(param_keys)))-1 , names=param_keys , dtype=["f"]*len(param_keys))
        for jj in range(len(src_compl_lr)):
            for jjj,key in enumerate(src_compl_lr[jj].getParamNames(),start=0):
                tab_this_cutout[key+".lr"][jj] = src_compl_lr[jj].getParams()[jjj]
                if key in ["pos.x","pos.y","brightness.Flux"]:
                    tab_this_cutout[key+".var.lr"][jj] = variance_table[key+".var"][jj]

        # Add new RA and DEC that is used for final Tractor fit.
        tab_this_cutout["RA_tractor.lr"] = np.zeros(len(tab_this_cutout))
        tab_this_cutout["DEC_tractor.lr"] = np.zeros(len(tab_this_cutout))
        for ii in range(len(tab_this_cutout)):
            tmp = img_cutout_wcs.all_pix2world([[tab_this_cutout["pos.x.lr"][ii],tab_this_cutout["pos.y.lr"][ii]]],0)
            tab_this_cutout["RA_tractor.lr"][ii] = tmp[0][0]
            tab_this_cutout["DEC_tractor.lr"][ii] = tmp[0][1]

        # Add small x/y offset to position
        tab_this_cutout["Dx_tractorfit.lr"] = Dx
        tab_this_cutout["Dy_tractorfit.lr"] = Dy

        # Add astrometry offsets in pixels
        tab_this_cutout["astro_offset_px_x.lr"] = np.repeat(delta_X_median_px,len(tab_this_cutout) )
        tab_this_cutout["astro_offset_px_y.lr"] = np.repeat(delta_Y_median_px,len(tab_this_cutout) )

        # combine LR and HR table
        TABLE_this_cutout = hstack([cat_models_hr_use_incutout,tab_this_cutout])

        # compute flux in muJy for LR and HR
        TABLE_this_cutout["brightness.Flux.mujy.hr"] = np.zeros(len(TABLE_this_cutout))
        TABLE_this_cutout["brightness.Flux.mujy.lr"] = np.zeros(len(TABLE_this_cutout))

        TABLE_this_cutout["brightness.Flux.mujy.hr"] = ConvertToMicroJansky(flux=TABLE_this_cutout["brightness.Flux.hr"], zp=userinput["hr_zp"] )
        TABLE_this_cutout["brightness.Flux.mujy.lr"] = ConvertToMicroJansky(flux=TABLE_this_cutout["brightness.Flux.lr"], zp=userinput["lr_zp"] )

        TABLE_this_cutout["brightness.Flux.var.mujy.hr"] = TABLE_this_cutout["brightness.Flux.mujy.hr"] * TABLE_this_cutout["brightness.Flux.var.hr"]/TABLE_this_cutout["brightness.Flux.hr"]
        TABLE_this_cutout["brightness.Flux.var.mujy.lr"] = TABLE_this_cutout["brightness.Flux.mujy.lr"] * TABLE_this_cutout["brightness.Flux.var.lr"]/TABLE_this_cutout["brightness.Flux.lr"]

        TABLE_this_cutout["brightness.Flux.mujy.hr"][TABLE_this_cutout["brightness.Flux.hr"] < 0] = -1
        TABLE_this_cutout["brightness.Flux.mujy.lr"][TABLE_this_cutout["brightness.Flux.lr"] < 0] = -1
        TABLE_this_cutout["brightness.Flux.var.mujy.hr"][TABLE_this_cutout["brightness.Flux.hr"] < 0] = -1
        TABLE_this_cutout["brightness.Flux.var.mujy.lr"][TABLE_this_cutout["brightness.Flux.lr"] < 0] = -1


        # add other stuff
        TABLE_this_cutout["fit_group.lr"] = np.repeat(group_counter,len(TABLE_this_cutout))
        TABLE_this_cutout["nbr_group_members.lr"] = np.repeat(len( np.where(TABLE_this_cutout["in_same_segmap_area"] == 1)[0] ),len(TABLE_this_cutout))

        # Add to large table
        sel_to_add = np.where( TABLE_this_cutout["in_same_segmap_area"] == 1 )[0]
        TABLE_this_cutout_add = TABLE_this_cutout.copy()[sel_to_add]
        TABLE = vstack([TABLE,TABLE_this_cutout_add])
        


        ## 8.8. Do the Bookkeeping ###########

        ids_to_add = list(TABLE_this_cutout["NUMBER.hr"][TABLE_this_cutout["in_same_segmap_area"] == 1])

        # add fitted objects to list and remove from OBJECTS_TO_FIT
        for idd in ids_to_add:

            OBJECTS_FITTED.append( idd )

            sel_tmp = np.where( OBJECTS_TO_FIT["NUMBER.hr"] == idd )[0]
            if len(sel_tmp) > 0:
                OBJECTS_TO_FIT.remove_row( int(sel_tmp) )

        #print(OBJECTS_FITTED)
        #print(list(OBJECTS_TO_FIT["NUMBER.hr"]))

        ## Advance group counter
        group_counter += 1


        ### DONE WITH ALL GALAXIES =========



    print(" done (in %g seconds)" % (round((time.time()-start_time),2)) )
    TIMESUMMARY[len(TIMESUMMARY)+1] = ["Running Tractor on low-res image object by object",round((time.time()-start_time),5)]




    ## 9. Add some additional stuff ==============

    # add global X/Y (do this outside of the loop to save time)
    xy_world = [ [TABLE["RA_tractor.lr"][jj] , TABLE["DEC_tractor.lr"][jj] ] for jj in range(len(TABLE))  ]
    xy = lr_img_wcs.all_world2pix(xy_world,1)
    TABLE["X_IMAGE_tractor.lr"] = xy[:,0]
    TABLE["Y_IMAGE_tractor.lr"] = xy[:,1]



    # Add flags
    # three types of flags:
    # 1. masked area on the LR image
    # 2. mask also objects that are not detected on the LR image.
    # 3. Close to the edge (within 1")

    lr_maskbad = mkmask(lr_mask,bitlist=[lr_mask_h["MP_NO_DATA"],lr_mask_h["MP_INTRP"],lr_mask_h["MP_SUSPECT"]]) # get areas without data and turn into NaN
    #lr_segmap_sextractor

    TABLE["flag_bad.lr"] = np.zeros(len(TABLE))
    TABLE["flag_detected.lr"] = np.zeros(len(TABLE))
    TABLE["flag_edge.lr"] = np.zeros(len(TABLE))
    for ii in range(len(TABLE)):
        if lr_maskbad[ int(TABLE["Y_IMAGE_tractor.lr"][ii]-1) , int(TABLE["X_IMAGE_tractor.lr"][ii]-1) ] != 0:
            TABLE["flag_bad.lr"][ii] = 1
        if lr_segmap_sextractor[ int(TABLE["Y_IMAGE_tractor.lr"][ii]-1) , int(TABLE["X_IMAGE_tractor.lr"][ii]-1) ] != 0:
            TABLE["flag_detected.lr"][ii] = 1
        if ( ( (TABLE["X_IMAGE_tractor.lr"][ii]-1) < (1/lr_pixscale) )
            | ( (TABLE["X_IMAGE_tractor.lr"][ii]-1) > (lr_img.shape[1]-1-(1/lr_pixscale)))
            | ( (TABLE["Y_IMAGE_tractor.lr"][ii]-1) < (1/lr_pixscale) )
            | ( (TABLE["Y_IMAGE_tractor.lr"][ii]-1) > (lr_img.shape[0]-1-(1/lr_pixscale)))
            ):
            TABLE["flag_edge.lr"][ii] = 1


    ## 10. Save Table  ==========
    TABLE.sort("NUMBER.hr") # sort again by number
    outname = os.path.join(dir_this_process,"lr_compl_table_final.fits")
    TABLE.write(outname, format='fits' , overwrite=True)



    ## 11. Create residual image ==========
    print("Creating residual image . . . ",end="")
    start_time = time.time()

    # create Tractor image
    if PSF_LR_TYPE == "fwhm":
        tim_model_lr = Image(data=np.zeros(lr_img.shape), invvar=np.ones_like(lr_img) / (lr_img_pixnoise**2),
                    psf=NCircularGaussianPSF([userinput["lr_image_psf"]/2.35/lr_pixscale], [1.]),
                    wcs=NullWCS(), photocal=NullPhotoCal(),
                    sky=ConstantSky(lr_img_medbkg))
    if PSF_LR_TYPE == "pixel":
        tim_model_lr = Image(data=np.zeros(lr_img.shape), invvar=np.ones_like(lr_img) / (lr_img_pixnoise**2),
                    psf=PixelizedPSF(PSF_LR_PIXEL),
                    wcs=NullWCS(), photocal=NullPhotoCal(),
                    sky=ConstantSky(lr_img_medbkg))


    # add sources
    try:
        src_model_lr.clear()
    except Exception as e:
        pass
    src_model_lr = []
    for ii in range(len(TABLE)):
        pos_init = PixPos(TABLE["X_IMAGE_tractor.lr"][ii]-1,TABLE["Y_IMAGE_tractor.lr"][ii]-1) # DONT FORGET THAT PYTHON STARTS AT 0
        flux_init = Flux(TABLE["brightness.Flux.lr"][ii])
        if TABLE["is_pointsource"][ii] == 2:
            thissource = PointSource(pos_init,flux_init)
        else:
            n_init = SersicIndex(TABLE["sersicindex.SersicIndex.lr"][ii])
            ba_init = TABLE["shape.ab.lr"][ii]
            re_init = TABLE["shape.re.lr"][ii]
            phi_init = TABLE["shape.phi.lr"][ii]
            thissource = SersicGalaxy(pos_init,
                                        flux_init,
                                        GalaxyShape(re_init, ba_init, phi_init), n_init)

        src_model_lr.append(thissource)


    # create tractor object
    tractor_model = Tractor([tim_model_lr] , src_model_lr)

    lr_mod = tractor_model.getModelImage(0)
    lr_res = lr_img - lr_mod

    print(" done (in %g seconds)" % (round((time.time()-start_time),2)) )
    TIMESUMMARY[len(TIMESUMMARY)+1] = ["Creating residual image",round((time.time()-start_time),2)]

    ## Save the images ====
    print("Writing LR FITS files . . . ",end="")
    start_time = time.time()

    save_lr_fits(img = lr_img,
                    mod = lr_mod,
                    res = lr_res,
                    seg = lr_segmap,
                    seg_sextractor = lr_segmap_sextractor,
                    mask_bad = lr_maskbad,
                    mask = lr_mask,
                    header=lr_img_h,
                    mask_header = lr_mask_h,
                    outfile=os.path.join(dir_this_process,"lr_tractor_results.fits")
                )

    print(" done (in %g seconds)" % (round((time.time()-start_time),2)) )
    TIMESUMMARY[len(TIMESUMMARY)+1] = ["Writing LR FITS files",round((time.time()-start_time),2)]


    ## Normalized sum of residuals
    v = lr_res.ravel()
    o = lr_img.ravel()
    normalized_sum_of_residuals = np.nansum(v**2)/len(o)
    print("Normalized sum of residuals = %5.5f" % (normalized_sum_of_residuals))
    STATS.append("Normalized sum of residuals = %5.5f" % (normalized_sum_of_residuals))


    print_time_summary(TIMESUMMARY)
    save_time_summary(TIMESUMMARY,file=os.path.join(dir_this_process,"lr_timesummary.txt") )
    save_stats(STATS,file=os.path.join(dir_this_process,"lr_stats.txt"))



    return(True)

                                

    '''##### Check #####
    fig = plt.figure(figsize=(10,5))
    ima = dict(cmap=plt.get_cmap("Greys"),vmin=-2*lr_img_pixnoise,vmax=10*lr_img_pixnoise,origin="lower",interpolation="nearest")

    ax1 = fig.add_subplot(1,2,1)
    ax1.imshow(img_cutout, **ima)
    ax1.plot(cat_models_hr_use_incutout["X_IMAGE_hr_cutout.lr"]-1 , cat_models_hr_use_incutout["Y_IMAGE_hr_cutout.lr"]-1,"x",color="red")
    ax1.plot(cat_models_hr_use_incutout["X_IMAGE_hr_cutout.lr"]-1 + delta_X_median_px , cat_models_hr_use_incutout["Y_IMAGE_hr_cutout.lr"]-1 + delta_Y_median_px,"d",color="red")
    for ii in range(len(src_compl_lr)):
        ax1.plot(src_compl_lr[ii].pos[0],src_compl_lr[ii].pos[1],"o",color="red",fillstyle="none")
        ax1.text(src_compl_lr[ii].pos[0],src_compl_lr[ii].pos[1]+2,cat_models_hr_use_incutout["NUMBER.hr"][ii],fontsize=10,color="red",va="bottom",ha="center")
    ax1.set_title("Original")

    ax2 = fig.add_subplot(1,2,2)
    ax2.imshow(lr_residual_compl,**ima)
    ax2.set_title("Residual")

    plt.show()

    ##################
'''
                             
