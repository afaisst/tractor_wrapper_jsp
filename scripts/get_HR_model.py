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

### GET HR MODEL ######
def get_HR_model(userinput,tileids):
    '''
    This function creates a model from the high-resolution image.

    USAGE: get_HR_model(userinput , tileids)
    where
        - userinput: is a user input dictionary
        - tileids: the IDs of the tiles on which the model should be created
    '''

    for tt in range(len(tileids)):
        this_tile_id = int(tileids[tt])


        ## 1. Create process ID ===========================
        #this_tile_id = 2

        TIMESUMMARY = dict()
        STATS = []
        
        process_id = "%04.0f_%s" % (this_tile_id , userinput["lr_large_image_name"].split("/")[-1].split(".fits")[0] ) # set the process ID to the name of the large image array + the tile number
        print("This Process is: " + str(process_id))
        STATS.append("This Process is: " + str(process_id))

        ## 2. create (temporary) working directory ===========================
        dir_this_process = os.path.join(userinput["workdir"], "%s_%s" % (userinput["out_prefix"],process_id) )
        sh.mkdir('-p', dir_this_process)


        ## 3. High resolution image ===========
        print("reading images . . .",end="")
        start_time = time.time()

        ## high resolution
        hr_image_path = os.path.join(userinput["cutout_output_dir"],userinput["lr_large_image_name"].split("/")[-1].split(".fits")[0],"%04.0f_%s.fits" % (this_tile_id,userinput["hr_large_image_name"].split("/")[-1].split(".fits")[0]))
        with fits.open(hr_image_path) as hdul:
            hr_img = hdul[0].data
            hr_img_h = hdul[0].header

        hr_pixscale = np.abs(hr_img_h["CD1_1"]*3600) # pixel scale in arcsec/px
        hr_img_wcs = wcs.WCS(hr_img_h)


        print(" done (in %g seconds)" % (round((time.time()-start_time),2)) )
        TIMESUMMARY[0] = ["Reading images",round((time.time()-start_time),2)]


        ## 3.5 Get the PSFs.

        # HR image
        if userinput["hr_psf_type"] == "fwhm":
            PSF_FWHM_ARCSEC_HR = userinput["hr_image_psf"]
            PSF_HR_TYPE = "fwhm"
            print("Using given FWHM and gaussian for PSF of high-res image.")
            STATS.append("Using given FWHM and gaussian for PSF of high-res image.")
        elif userinput["hr_psf_type"] == "fits":
            PSF_HR_TYPE = "pixel"
            with fits.open(os.path.join(userinput["hr_image_psf"])) as hdul:
                PSF_HR_PIXEL = hdul[0].data
                PSF_HR_PIXEL = PSF_HR_PIXEL / np.nansum(PSF_HR_PIXEL)
            print("Using pixel PSF for high-res image in FITS format.")
            STATS.append("Using pixel PSF for high-res image in FITS format.")
        elif userinput["hr_psf_type"] == "psfex":
            PSF_HR_TYPE = "pixel"
            with fits.open(os.path.join(userinput["hr_image_psf"])) as hdul:
                PSF_HR_PIXEL = modelPSF_1D(psfcube=hdul[1].data,psfcube_h=hdul[1].header,param=21) #22
                PSF_HR_PIXEL = PSF_HR_PIXEL / np.nansum(PSF_HR_PIXEL)
            print("Using pixel PSF for high-res image in PSFex format.")
            STATS.append("Using pixel PSF for high-res image in PSFex format.")
        else:
            print("PSF type not understood.")
            STATS.append("PSF type not understood.")
            quit()

        ## 4. Run SExtractor ===========================


        ## 4.1 Copy SExtractor config file ##
        config_default = os.path.join(userinput["sexinputdir"],"default.conf.interactive")
        config_this_process = os.path.join(dir_this_process,"hr_default.conf")
        cmd = "cp " + config_default + " " + config_this_process
        subprocess.run(cmd, shell=True)


        # 4.2 Adjust configuration file ##
        sex_output_cat_hr_file = os.path.join(dir_this_process,"hr_sex_output.cat")
        PARS_this_process = get_default_parfile()
        PARS_this_process["CATALOG_NAME"] = sex_output_cat_hr_file
        PARS_this_process["PARAMETERS_NAME"] = os.path.join(userinput["sexinputdir"],"sex.par")
        PARS_this_process["FILTER_NAME"] = os.path.join(userinput["sexinputdir"],"g2.8.conv")
        PARS_this_process["STARNNW_NAME"] = os.path.join(userinput["sexinputdir"],"default.nnw")
        PARS_this_process["CHECKIMAGE_NAME"] = "%s" % os.path.join(dir_this_process,"hr_seg.fits")
        PARS_this_process["CHECKIMAGE_TYPE"] = "SEGMENTATION"
        PARS_this_process["PHOT_APERTURES"] = ','.join(map(str, get_apertures(hr_pixscale,aperturelist_arcsec=[3.0]).tolist())) # Note: if number of apertures are changes, also change sex.par file!!!
        PARS_this_process["PIXEL_SCALE"] = hr_pixscale
        PARS_this_process["DEBLEND_MINCONT"] = 0.1 #0.01
        PARS_this_process["DEBLEND_NTHRESH"] = 32 #0.01
        PARS_this_process["DETECT_MINAREA"] = 20 # 5
        PARS_this_process["DETECT_THRESH"] = 2.0 #2
        PARS_this_process["MAG_ZEROPOINT"] = userinput["hr_zp"]
        PARS_this_process["ANALYSIS_THRESH"] = 1.5#1.5
        PARS_this_process["BACKPHOTO_TYPE"] = "LOCAL"  #GLOBAL
        if PSF_HR_TYPE == "fwhm": #seeing has to be in arcsec
            PARS_this_process["SEEING_FWHM"] = userinput["hr_image_psf"]
        else:
            PARS_this_process["SEEING_FWHM"] = 0.15 # Doesn't really matter for SExtractor

        for key in PARS_this_process.keys():
            replace_in_file(filename = config_this_process,
                           old_string = "*" + key + "*",
                           new_string = str(PARS_this_process[key]))


        # 4.3 Run SExtractor and load catalog ##
        start_time = time.time()
        print("running SExtractor . . .",end="")
        cmd = "%s %s -c %s " % (userinput["sex_command"],
                                                                      hr_image_path,
                                                                   config_this_process
                                                                                  )
        subprocess.run(cmd , shell=True)
        print(" done (in %g minutes)" % (round((time.time()-start_time)/60,2)) )
        TIMESUMMARY[len(TIMESUMMARY)+1] = ["Running SExtractor on high-resolution image",round((time.time()-start_time),2)]


        # 4.4 Read SExtractor catalog and change column names to .hr and add "stellarity flag"
        sexcat_hr = ascii.read(sex_output_cat_hr_file)
        for key in sexcat_hr.keys():
            sexcat_hr.rename_column(key, key + ".hr" )
        sexcat_hr.sort("FLUX_AUTO.hr") # sort by flux (faint first).
        sel_stars = np.where( (sexcat_hr["CLASS_STAR.hr"] > 0.8)
                             | ((-2.5*np.log10(sexcat_hr["FLUX_APER.hr"])+userinput["hr_zp"] < 19)
                                & (sexcat_hr["B_IMAGE.hr"]/sexcat_hr["A_IMAGE.hr"]>0.8)
                               )
                             | (sexcat_hr["FLUX_RADIUS.hr"] < userinput["hr_psf_fwhm_arcsec"]/hr_pixscale)
                            )[0]
        sexcat_hr["is_pointsource"] = np.zeros(len(sexcat_hr))
        sexcat_hr["is_pointsource"][sel_stars] = 2
        print("Number of galaxies in high-resolution SExtractor catalog: %g" % len(sexcat_hr))
        STATS.append("Number of galaxies in high-resolution SExtractor catalog: %g" % len(sexcat_hr))

        # 4.45 Read SExtractor Segmentation map
        with fits.open(os.path.join(dir_this_process,"hr_seg.fits") ) as hdul:
            hr_seg = hdul[0].data

        # 4.5 Create DS9 file with extractions ##
        hr_ds9radec_name = dir_this_process + "/hr_radec.reg"
        with open(hr_ds9radec_name,"w+") as f:
            f.write('global color=magenta dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
            f.write('fk5\n')
            for jj in range(len(sexcat_hr)):
                if sexcat_hr["is_pointsource"][jj] > 0:
                    col_this = "cyan"
                else:
                    col_this = "magenta"
                f.write('circle(%s,%s,0.5") # color=%s width=2 text={%s}\n' % (sexcat_hr["ALPHA_J2000.hr"][jj] , sexcat_hr["DELTA_J2000.hr"][jj], col_this, sexcat_hr["NUMBER.hr"][jj]) )

        

        ## Get some properties of the images ================

        ## Check some stuff
        if "CD1_2" not in hr_img_h.keys():
            hr_img_h["CD1_2"] = 0
        if "CD2_1" not in hr_img_h.keys():
            hr_img_h["CD2_1"] = 0

        ## Get Pixel noise ##
        tmp = clip(hr_img[hr_seg == 0], n=3, niter=10)
        hr_img_pixnoise = tmp["stdev"]
        hr_img_medbkg = tmp["med"]
        print("Pixelnoise of high-res image: %g" % hr_img_pixnoise )
        print("Median background of high-res image: %g" % hr_img_medbkg )
        STATS.append("Pixelnoise of high-res image: %g" % hr_img_pixnoise )
        STATS.append("Median background of high-res image: %g" % hr_img_medbkg)


        ## Now, go through each object.
        OBJECTS_TO_FIT = Table(data=np.asarray(sexcat_hr["NUMBER.hr"]).reshape((len(sexcat_hr),1)),names=["NUMBER.hr"]) # store here the NUMBERs of the objects that still have to be fitted
        OBJECTS_FITTED = [] # store there NUMBERs of objects that have been fitted.
        OBJECTS_FAILED = [] # store there NUMBERs of objects for which the fit failed.
        
        ## DEBUGGING
        #sel_debug = np.where(OBJECTS_TO_FIT["NUMBER.hr"] == 98)[0] ######## ADDED
        #OBJECTS_TO_FIT = OBJECTS_TO_FIT[sel_debug] ######## ADDED
        
        TABLE = Table()

        group_counter = 0
        start_time = time.time()
        #for cc in range(1): # change this to while loop later
        while len(OBJECTS_TO_FIT) > 0:

            main_obj_number = OBJECTS_TO_FIT["NUMBER.hr"][0] # always take next galaxy in line
            print("\nMain Object: %g (%g objects left)" % (main_obj_number,len(OBJECTS_TO_FIT)) )

            id_in_sexcat = np.where( sexcat_hr["NUMBER.hr"] == main_obj_number )[0][0]
            
            ## extract image part that should be fitted. ==============
            # Don't forget new WCS!!!

            # get extension (note that SExtractor starts with 1 while Python starts with 0)
            xmin = sexcat_hr["XMIN_IMAGE.hr"][id_in_sexcat]-1
            ymin = sexcat_hr["YMIN_IMAGE.hr"][id_in_sexcat]-1
            xmax = sexcat_hr["XMAX_IMAGE.hr"][id_in_sexcat]-1
            ymax = sexcat_hr["YMAX_IMAGE.hr"][id_in_sexcat]-1

            
            # get some extra space
            xmin = xmin - (xmax-xmin)*0.5
            xmax = xmax + (xmax-xmin)*0.5
            ymin = ymin - (ymax-ymin)*0.5
            ymax = ymax + (ymax-ymin)*0.5
            
            # get final numbers (new: Aug 4, 2020)
            x_extension = np.float(xmax - xmin)
            y_extension = np.float(ymax - ymin)
            #x_extension = np.nanmax( np.asarray([np.float(xmax - xmin) , 1.0/hr_pixscale]) ) # make cutouts at least 1 arcsec across
            #y_extension = np.nanmax( np.asarray([np.float(ymax - ymin) , 1.0/hr_pixscale]) ) # make cutouts at least 1 arcsec across
            x_center = xmin+(xmax-xmin)/2 # this is the same as before (only center counts)
            y_center = ymin+(ymax-ymin)/2 # this is the same as befor (only center counts)


            # for cutout: X and Y reversed and Python starts with 0 (therefore X_IMAGE-1 and Y_IMAGE-1)!!
            #tmp = Cutout2D(data=hr_img,position=( xmin+(xmax-xmin)/2,ymin+(ymax-ymin)/2),size=(ymax-ymin , xmax-xmin ) , copy=True , mode="trim",wcs=hr_img_wcs)
            tmp = Cutout2D(data=hr_img,position=( x_center , y_center) , size=(y_extension , x_extension ) , copy=True , mode="trim",wcs=hr_img_wcs)
            img_cutout = tmp.data.copy()
            img_cutout_wcs = tmp.wcs.copy()

            #tmp = Cutout2D(data=hr_seg,position=( xmin+(xmax-xmin)/2,ymin+(ymax-ymin)/2),size=(ymax-ymin , xmax-xmin ) , copy=True , mode="trim",wcs=hr_img_wcs)
            tmp = Cutout2D(data=hr_seg,position=( x_center , y_center) , size=(y_extension , x_extension ) , copy=True , mode="trim",wcs=hr_img_wcs)
            hr_seg_cutout = tmp.data.copy()
            

            ## Fit that small image with Tractor ============
            print(np.unique(hr_seg_cutout.ravel()))
            
            ## get objects that are in cutout
            id_objects_in_cutout = list(np.unique(hr_seg_cutout.ravel()))# id_objects_in_cutout.copy()
            if [0] in id_objects_in_cutout:
                id_objects_in_cutout.remove(0) # remove 0
            sel_objects_in_cutout = np.asarray([ np.where( sexcat_hr["NUMBER.hr"] == idd )[0][0] for idd in id_objects_in_cutout])
            sexcat_hr_in_cutout = sexcat_hr[sel_objects_in_cutout].copy()
            print("Number of sources in cutout: %g" % len(sexcat_hr_in_cutout) )
            print("Objects to be fitted: %s" % str(list(sexcat_hr_in_cutout["NUMBER.hr"])) )
            
            ## Check some things before running ------
            # This was added on July 20 2020 when a (probably?) bug in SExtractor was discovered
            # when running job_1008 (9813, 6_7, tile 0008). Specifically, the NUMBER 1043 (although in the
            # SExtractor catalog) was *not* on the Segmentation map!. I added the following check
            # to prevent the script to crash. If it doesn't find a NUMBER on the segmentation map, it
            # simply skips the galaxy. This happened in SExtractor version 2.19.1
            FATAL_FAIL = 0
            
            ## 1. Check if there is an object at all.
            # If not, skip all below and add the object to OBJECTS_FAILED (remove from OBJECTS_TO_FIT).
            if len(sexcat_hr_in_cutout) == 0: # check if there is an object in the cutout at all.
                print("NO SOURCES IN CUTOUTS - SKIP")
                FATAL_FAIL = 1
                OBJECTS_FAILED.append(main_obj_number)
                sel_tmp = np.where( OBJECTS_TO_FIT["NUMBER.hr"] == main_obj_number )[0]
                if len(sel_tmp) > 0:
                    OBJECTS_TO_FIT.remove_row( int(sel_tmp) )
               
            ## 2. Check if main object is in cutout.
            # It happened at least ones that there is no segmentation for an object NUMBER.
            # This might be a big in SExtractor. Therefore, check here if the main object
            # is in the cutout. If not, go straight to failed and skip the rest.
            if FATAL_FAIL == 0: # only do this if there is an object, else don't even look at it.
                if main_obj_number not in np.unique(hr_seg_cutout):
                    print("MAIN OBJECT %g IS NOT IN THIS SEGMENTATION CUTOUT - this might be bug in SExtractor: skip this object" % main_obj_number)
                    FATAL_FAIL = 1
                    OBJECTS_FAILED.append(main_obj_number)
                    sel_tmp = np.where( OBJECTS_TO_FIT["NUMBER.hr"] == main_obj_number )[0]
                    if len(sel_tmp) > 0:
                        OBJECTS_TO_FIT.remove_row( int(sel_tmp) )

            
            ## If no FATAL_FAIL, continue
            if FATAL_FAIL == 0: # continue here with the rest for this object of FATAL_FAIL == 0

                ## Get percentage of the objects in cutout (if not fully included, fit again and don't add to table at the end!)
                coverage = Table(np.zeros((len(id_objects_in_cutout),4)) , names=["id","nbr_pix_tot","nbr_pix_cutout","ratio"] , dtype=["int","int","int","f"] )
                for iii,idd in enumerate(id_objects_in_cutout,start=0):
                    coverage["id"][iii] = idd
                    coverage["nbr_pix_tot"][iii] = len( np.where( hr_seg.ravel() == idd )[0] )
                    coverage["nbr_pix_cutout"][iii] = len( np.where( hr_seg_cutout.ravel() == idd )[0] )
                    coverage["ratio"][iii] = coverage["nbr_pix_cutout"][iii]/coverage["nbr_pix_tot"][iii]


                # create source list
                try:
                    src_compl_hr.clear()
                except Exception as e:
                    pass
                src_compl_hr = []
                for ii in range(len(sexcat_hr_in_cutout)):

                    #pos_init = PixPos(sexcat_hr_in_cutout["X_IMAGE.hr"][ii],sexcat_hr_in_cutout["Y_IMAGE.hr"][ii])
                    tmp = img_cutout_wcs.all_world2pix([[sexcat_hr_in_cutout["ALPHA_J2000.hr"][ii],sexcat_hr_in_cutout["DELTA_J2000.hr"][ii]]],0)
                    pos_init = PixPos(tmp[0][0],tmp[0][1])
                    pos_init.lowers = [tmp[0][0]-1,tmp[0][1]-1]
                    pos_init.uppers = [tmp[0][0]+1,tmp[0][1]+1]
                    flux_init = Flux(sexcat_hr_in_cutout["FLUX_AUTO.hr"][ii])
                    #tmp = hr_img_wcs.all_world2pix([[sexcat_hr["ALPHA_J2000.hr"][ii],sexcat_hr["DELTA_J2000.hr"][ii]]],0)

                    #if sexcat_hr_in_cutout["is_pointsource"][ii] == 20: # Do not use pointsource fitting in HR image
                    if sexcat_hr_in_cutout["is_pointsource"][ii] == 2: # Use pointsource fitting in HR image
                        print("Fitting as point source!")
                        thissource = PointSource(pos_init , flux_init )
                    else:
                        n_init = SersicIndex(2.0)
                        n_init.lower = 0
                        n_init.upper = 20
                        ba_init = sexcat_hr_in_cutout["B_IMAGE.hr"][ii]/sexcat_hr_in_cutout["A_IMAGE.hr"][ii]
                        re_init = sexcat_hr_in_cutout["FLUX_RADIUS.hr"][ii]/np.sqrt(ba_init)  # semimajor axis!
                        phi_init = sexcat_hr_in_cutout["THETA_IMAGE.hr"][ii]-90.0
                        thissource = SersicGalaxy(pos_init,
                                                    flux_init,
                                                    GalaxyShape(re_init, ba_init, phi_init), n_init)

                    src_compl_hr.append(thissource)


                # create Tractor image
                if PSF_HR_TYPE == "fwhm":
                    tim_compl_hr = Image(data=img_cutout.astype(np.float32), invvar=np.ones_like(img_cutout) / (hr_img_pixnoise**2),
                                psf=NCircularGaussianPSF([userinput["hr_image_psf"]/2.35/hr_pixscale], [1.]),
                                wcs=NullWCS(), photocal=NullPhotoCal(),
                                sky=ConstantSky(hr_img_medbkg))
                if PSF_HR_TYPE == "pixel":
                    tim_compl_hr = Image(data=img_cutout.astype(np.float32), invvar=np.ones_like(img_cutout) / (hr_img_pixnoise**2),
                                psf=PixelizedPSF(PSF_HR_PIXEL.astype(np.float32)),
                                wcs=NullWCS(), photocal=NullPhotoCal(),
                                sky=ConstantSky(hr_img_medbkg))


                # create tractor object
                tractor_compl_hr = Tractor([tim_compl_hr] , src_compl_hr , optimizer=ConstrainedOptimizer())
                #tractor_compl_hr = Tractor([tim_compl_hr] , src_compl_hr)

                # freeze image calibration params
                tractor_compl_hr.freezeParam('images')

                for ii in range(len(src_compl_hr)):
                    if src_compl_hr[ii].getSourceType().lower() == "pointsource":
                        #src_compl_hr[ii].freezeParams('pos')
                        pass
                    else:
                        #src_compl_hr[ii].freezeParams('pos')
                        #src_compl_hr[ii].freezeParams("sersicindex")
                        src_compl_hr[ii].shape.freezeParams("phi")

                print("    Fitting complex models ",end="")
                for ii in range(userinput["hr_compl_niter"]):
                    try:
                        dlnp,X,alpha,variances = tractor_compl_hr.optimize(variance=True)  # NEW: ADDED VARIANCE
                        print(str(ii)+" ("+str(dlnp)+") ",end=" ")
                        FIT_SUCCESS_1 = True
                        FIT_SUCCESS_2 = True
                        iter_max = ii
                        if dlnp < 1e-3:
                            break
                    except Exception as e: # first fit failed -> set FIT_SUCCESS_1 to False
                        FIT_SUCCESS_1 = False
                        iter_max = ii
                        break

                #print("\n")
                #print(src_compl_hr)


                ## If not fitted, try to fit up to the iteration where it failed. Not perfec (indicate that somewhere) but maybe good enough.
                # Note that if number of iterations = 0, than this step won't work and abort.
                # Have to create a new src list because copy does not work. 
                if (not FIT_SUCCESS_1):
                    if iter_max > 0:
                        print("\nFIT FAILED AT iter_max = %g (try again to that iteration)" % iter_max)
                        src_compl_hr_2 = []
                        for ii in range(len(sexcat_hr_in_cutout)):

                            tmp = img_cutout_wcs.all_world2pix([[sexcat_hr_in_cutout["ALPHA_J2000.hr"][ii],sexcat_hr_in_cutout["DELTA_J2000.hr"][ii]]],0)
                            pos_init = PixPos(tmp[0][0],tmp[0][1])
                            pos_init.lowers = [tmp[0][0]-1,tmp[0][1]-1]
                            pos_init.uppers = [tmp[0][0]+1,tmp[0][1]+1]
                            #pos_init = PixPos(sexcat_hr_in_cutout["X_IMAGE.hr"][ii],sexcat_hr_in_cutout["Y_IMAGE.hr"][ii])
                            flux_init = Flux(sexcat_hr_in_cutout["FLUX_AUTO.hr"][ii])

                            #if sexcat_hr_in_cutout["is_pointsource"][ii] == 20: # Do not fit point source in HR image
                            if sexcat_hr_in_cutout["is_pointsource"][ii] == 2: # Fit point source in HR image
                                thissource = PointSource(pos_init , flux_init )
                            else:
                                n_init = SersicIndex(2.0)
                                n_init.lower = 0
                                n_init.upper = 20
                                ba_init = sexcat_hr_in_cutout["B_IMAGE.hr"][ii]/sexcat_hr_in_cutout["A_IMAGE.hr"][ii]
                                re_init = sexcat_hr_in_cutout["FLUX_RADIUS.hr"][ii]/np.sqrt(ba_init)  # semimajor axis!
                                phi_init = sexcat_hr_in_cutout["THETA_IMAGE.hr"][ii]-90.0
                                thissource = SersicGalaxy(pos_init,
                                                            flux_init,
                                                            GalaxyShape(re_init, ba_init, phi_init), n_init)
                            src_compl_hr_2.append(thissource)


                        tractor_compl_hr = Tractor([tim_compl_hr] , src_compl_hr_2 , optimizer=ConstrainedOptimizer())
                        tractor_compl_hr.freezeParam('images')

                        for ii in range(len(src_compl_hr_2)):
                            if src_compl_hr_2[ii].getSourceType().lower() == "pointsource":
                                #src_compl_hr_2[ii].freezeParams('pos')
                                pass
                            else:
                                #src_compl_hr_2[ii].freezeParams('pos')
                                #src_compl_hr_2[ii].freezeParams("sersicindex")
                                src_compl_hr_2[ii].shape.freezeParams("phi")

                        for ii in range(iter_max): # changed to iter_max (removed "-1" because range(N) is going up to N-1)
                            try:
                                dlnp,X,alpha,variances = tractor_compl_hr.optimize(variance=True)
                                print(str(ii)+" ("+str(dlnp)+") ",end=" ")
                                FIT_SUCCESS_2 = True
                                iter_max = ii
                                if dlnp < 1e-3:
                                    break
                            except Exception as e: # second fit failes -> set FIT_SUCCESS_2 to False
                                FIT_SUCCESS_2 = False
                                iter_max = -99
                    else:
                        FIT_SUCCESS_2 = False
                        iter_max = -99
                
                
                
                if not FIT_SUCCESS_2:
                    print("FIT FAILED. MOVING ON")

                    ## Remove main object from the list (bc that might have cause the trouble)
                    OBJECTS_FAILED.append( main_obj_number )
                    sel_tmp = np.where( OBJECTS_TO_FIT["NUMBER.hr"] == main_obj_number )[0]

                    if len(sel_tmp) > 0:
                        OBJECTS_TO_FIT.remove_row( int(sel_tmp) )
                    continue

                elif (not FIT_SUCCESS_1) & (FIT_SUCCESS_2): # failed first time, success second time
                    src_compl_hr.clear() # clear
                    src_compl_hr = src_compl_hr_2.copy() # now copy final back to original name
                else:
                    pass
                
                

	            ## IF ALL IS GOOD (i.e., fit success in first or second round): =========

                ## Create table with variances (NEW)
                # DO THIS BEFORE THAWNING PARAMETERS!
                # Note: Variances array comes in packages for each parameter (e.g., variances = x,x,x , y,y,y , z,z,z). These 
                # packages might not have same length as the parameters fitted for each sources may vary.
                # Need to use pop() to remove number of parameters.
                variances = list(variances) # make list
                variance_table = Table(data=np.zeros((len(src_compl_hr),7))-99, names=["pos.x.var","pos.y.var",
                                                "brightness.Flux.var","shape.re.var","shape.ab.var",
                                                "shape.phi.var","sersicindex.SersicIndex.var"])
                for sss,src in enumerate(src_compl_hr):
                    params_fitted = src.getParamNames()
                    params_variances = [variances.pop(0) for vvv in range(len(params_fitted))]
                    for ppp,param in enumerate(params_fitted):
                        variance_table["%s.var" % param][sss] = params_variances[ppp]




                ## Create table for Tractor output values
                # NEW: also add variances
                for ii in range(len(src_compl_hr)):
                    src_compl_hr[ii].thawAllRecursive() # Thawn all parameters so that we can see them.

                param_keys = ['pos.x',
                                'pos.y',
                                'brightness.Flux',
                                'shape.re','shape.ab',
                                'shape.phi',
                                'sersicindex.SersicIndex']
                param_keys = param_keys + [key + ".var" for key in param_keys]
                param_keys = [key + ".hr" for key in param_keys]
                this_table_tractor = Table(data=np.zeros((len(src_compl_hr),len(param_keys)))-1 , names=param_keys , dtype=["f"]*len(param_keys))
                this_table = Table()
                for jj in range(len(src_compl_hr)):
                    for jjj,key in enumerate(src_compl_hr[jj].getParamNames(),start=0):
                        this_table_tractor[key+".hr"][jj] = src_compl_hr[jj].getParams()[jjj]
                        this_table_tractor[key+".var.hr"][jj] = variance_table[key+".var"][jj]
                this_table = hstack([sexcat_hr_in_cutout,this_table_tractor])

                ## Add Success flags
                this_table["fit_success_1"] = np.repeat(int(FIT_SUCCESS_1),len(this_table))
                this_table["fit_success_2"] = np.repeat(int(FIT_SUCCESS_2),len(this_table))

                ## Sum up residual in segmentation map.
                # first create residual quickly
                mod_this_cutout = tractor_compl_hr.getModelImage(0)
                res_this_cutout = img_cutout - mod_this_cutout
                res_sums = []
                for ccc in range(len(this_table)):
                    mask_tmp = np.zeros(hr_seg_cutout.shape)
                    mask_tmp[hr_seg_cutout == this_table["NUMBER.hr"][ccc]] = 1
                    res_sums.append( 1 * np.sqrt( np.nansum( mask_tmp * res_this_cutout**2 ) ) )
                this_table["res_sum_square.hr"] = res_sums



                ## from this table, remove the objects that were in the FoV but not fully contained
                # request that objects must be 100% contained. Else fit again.
                ids_not_fully_contained = coverage["id"][coverage["ratio"] < 1]
                if len(ids_not_fully_contained) > 0:
                    for idd in ids_not_fully_contained:
                        this_table.remove_rows( np.where( this_table["NUMBER.hr"] == idd )[0] )


                ## from this table from this run, remove the objects that are already fitted (do something smarter later.)
                sel_remove = []
                for ii in range(len(this_table)):
                    try: # have to try here because in the first row TABLE is undefined.
                        if this_table["NUMBER.hr"][ii] in TABLE["NUMBER.hr"]:
                            sel_tmp = np.where( this_table["NUMBER.hr"][ii] == TABLE["NUMBER.hr"] )[0]

                            # check if residual looks better/worse than what already exists in table.
                            if this_table["res_sum_square.hr"][ii] > TABLE["res_sum_square.hr"][sel_tmp]:
                                sel_remove.append(ii)
                            else:
                                TABLE.remove_rows(sel_tmp) # if better in this fit, remove from the large table and add later value from new table.

                    except Exception as e:
                        pass

                if len(sel_remove) < len(this_table): # only do the following if there is still something in the table.

                    this_table.remove_rows(sel_remove)

                    # add RA/DEC
                    xy = [ [this_table["pos.x.hr"][jj] , this_table["pos.y.hr"][jj] ] for jj in range(len(this_table))  ]
                    xy_world = img_cutout_wcs.all_pix2world(xy,0)
                    this_table["RA_tractor.hr"] = xy_world[:,0]
                    this_table["DEC_tractor.hr"] = xy_world[:,1]

                    # add other stuff
                    this_table["iter_max.hr"] = np.repeat(iter_max,len(this_table))
                    this_table["fit_group.hr"] = np.repeat(group_counter,len(this_table))
                    this_table["nbr_group_members.hr"] = np.repeat(len(this_table),len(this_table))

                    # Add to large table
                    TABLE = vstack([TABLE,this_table])

                    # add fitted objects to list and remove from OBJECTS_TO_FIT
                    for ii in range(len(this_table)):
                        OBJECTS_FITTED.append( this_table["NUMBER.hr"][ii] )
                        sel_tmp = np.where( OBJECTS_TO_FIT["NUMBER.hr"] == this_table["NUMBER.hr"][ii] )[0]
                        if len(sel_tmp) > 0:
                            OBJECTS_TO_FIT.remove_row( int(sel_tmp) )

                else:
                    print("All objects are already in the table.")



                ## Advance group counter
                group_counter += 1

                ## not really necessary
                #hr_mod_compl = tractor_compl_hr.getModelImage(0)
                #hr_chi2_compl = tractor_compl_hr.getChiImage(0)
                #hr_residual_compl = img_cutout - hr_mod_compl



        # add global X/Y (do this outside of the loop to save time)
        xy_world = [ [TABLE["RA_tractor.hr"][jj] , TABLE["DEC_tractor.hr"][jj] ] for jj in range(len(TABLE))  ]
        xy = hr_img_wcs.all_world2pix(xy_world,1)
        TABLE["X_IMAGE_tractor.hr"] = xy[:,0]
        TABLE["Y_IMAGE_tractor.hr"] = xy[:,1]

        ## ADD FAILED OBJECTS TO CATALOG
        for obj in OBJECTS_FAILED:
            tmp = TABLE[0:1].copy() # just pick one to copy
            id_failed_in_sexcat = np.where( (sexcat_hr["NUMBER.hr"] == obj) )
            for key in tmp.colnames:
                if key in sexcat_hr.keys():
                    tmp[key] = sexcat_hr[key][id_failed_in_sexcat]
                else:
                    tmp[key] = -99
            TABLE = vstack([TABLE,tmp])


        print(" DONE (in %g minutes)" % (round((time.time()-start_time)/60,2)) )
        TIMESUMMARY[len(TIMESUMMARY)+1] = ["Running Tractor on high-resolution image",round((time.time()-start_time),5)]


        ## SAVE TABLE  ==========
        TABLE.sort("NUMBER.hr") # sort again by number
        outname = dir_this_process + "/" + "hr_compl_table_final.fits"
        TABLE.write(outname, format='fits' , overwrite=True)


        ## CREATE RESIDUAL IMAGE ==========
        print("Creating residual image . . . ",end="")
        start_time = time.time()
        
        # create Tractor image
        if PSF_HR_TYPE == "fwhm":
            tim_model_hr = Image(data=np.zeros(hr_img.shape), invvar=np.ones_like(hr_img) / (hr_img_pixnoise**2), # CHANGE: img_cutout to hr_img
                        psf=NCircularGaussianPSF([userinput["hr_image_psf"]/2.35/hr_pixscale], [1.]),
                        wcs=NullWCS(), photocal=NullPhotoCal(),
                        sky=ConstantSky(0.))
        if PSF_HR_TYPE == "pixel":
            tim_model_hr = Image(data=np.zeros(hr_img.shape), invvar=np.ones_like(hr_img) / (hr_img_pixnoise**2), # CHANGE: img_cutout to hr_img
                        psf=PixelizedPSF(PSF_HR_PIXEL.astype(np.float32)),
                        wcs=NullWCS(), photocal=NullPhotoCal(),
                        sky=ConstantSky(0.))

        # add sources
        try:
            src_model_hr.clear()
        except Exception as e:
            pass
        src_model_hr = []
        for ii in range(len(TABLE)):
            pos_init = PixPos(TABLE["X_IMAGE_tractor.hr"][ii]-1,TABLE["Y_IMAGE_tractor.hr"][ii]-1) # DONT FORGET THAT PYTHON STARTS AT 0
            flux_init = Flux(TABLE["brightness.Flux.hr"][ii])
            #if TABLE["is_pointsource"][ii] == 20:
            if TABLE["is_pointsource"][ii] == 2:
                thissource = PointSource(pos_init,flux_init)
            else:
                n_init = SersicIndex(TABLE["sersicindex.SersicIndex.hr"][ii])
                ba_init = TABLE["shape.ab.hr"][ii]
                re_init = TABLE["shape.re.hr"][ii]
                phi_init = TABLE["shape.phi.hr"][ii]
                thissource = SersicGalaxy(pos_init,
                                          flux_init,
                                          GalaxyShape(re_init, ba_init, phi_init), n_init)

            src_model_hr.append(thissource)


        # create tractor object
        tractor_model = Tractor([tim_model_hr] , src_model_hr)

        hr_mod = tractor_model.getModelImage(0)
        hr_res = hr_img - hr_mod
        
        print(" done (in %g seconds)" % (round((time.time()-start_time),2)) )
        TIMESUMMARY[len(TIMESUMMARY)+1] = ["Creating residual image",round((time.time()-start_time),2)]
        
        v = hr_res.ravel()
        o = hr_img.ravel()
        normalized_sum_of_residuals = np.nansum(v**2)/len(o)
        STATS.append("Normalized sum of residuals = %5.5f" % (normalized_sum_of_residuals))
        
        
        ## SAVE ALL STUFF ====
        print("Writing HR FITS files . . . ",end="")
        start_time = time.time()

        save_hr_fits(img=hr_img,
                     mod=hr_mod,
                     res=hr_res,
                     seg=hr_seg,
                     header=hr_img_h,
                     outfile=os.path.join(dir_this_process,"hr_tractor_results.fits")
                    )

        print(" done (in %g seconds)" % (round((time.time()-start_time),2)) )
        TIMESUMMARY[len(TIMESUMMARY)+1] = ["Writing HR FITS files",round((time.time()-start_time),2)]



        
        
        print_time_summary(TIMESUMMARY)
        save_time_summary(TIMESUMMARY,file=os.path.join(dir_this_process,"hr_timesummary.txt") )
        save_stats(STATS,file = os.path.join(dir_this_process,"hr_stats.txt"))



        print("\n\n SUMMARY")
        print("Objects still to fit:")
        print(list(OBJECTS_TO_FIT["NUMBER.hr"]))

        print("Objects sucessfully fitted:")
        print(OBJECTS_FITTED)

        print("Objects failed to fit:")
        print(OBJECTS_FAILED)
        
    ## End for each tile id -------
    
## End for function ----
