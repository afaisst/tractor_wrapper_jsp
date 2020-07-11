## This scripts goes through all the input directories and creates job JSON files that can be run separately.

import os, sys
import numpy as np
import json
import glob
import copy

import sh

#from astropy.io import fits, ascii
#from astropy.table import Table, Column, MaskedColumn, hstack, vstack



#### USER INPUT ####

## 1. Directories where to find the files and output directory
directories = {"hr_dir":"/stage/irsa-jointproc-data03/ACS_COSMOS/from_irsa/downloads/", # containing HR image (patch)
			"lr_dir":"/stage/irsa-jointproc-data03/ACS_COSMOS/grizli-stacks/HSC_corresponding_calexp_oct152019/", # containing LR image (patch)
			"hr_psf_dir":"/stage/irsa-jointproc-data03/ACS_COSMOS/from_irsa/sextractor/psf/", # containing HR PSF (patch)
			"lr_psf_dir":"/stage/irsa-jointproc-data02/PSFex/allcalexp_psfex_pipeline/output_variableMag/psfs/", # containing LR PSF (patch)
			"hr_astrometry_dir":"/stage/irsa-jointproc-data03/ACS_COSMOS/from_irsa/catalogs_analys/", # containing LR astrometry correction (patch)
			"lr_astrometry_dir":"/stage/irsa-jointproc-data02/pdr1_udeep/warp_calexp_output/allcalexp_sextractor_output/catalogs_analys/", # containing LR astrometry correction (patch)
			"job_output_dir":"/stage/irsa-jointproc-data03/TRACTOR/run_acs_jsp_full_SB4_Jul20/jobs/" # output directory where jobs are saved to disk
             }

## 2. Template file (containing all the other information)
'''template_userinput = {"mainpath":"/stage/irsa-jointproc-data03/TRACTOR/run_acs_jsp_full_SB4_Jan20/SANDBOX4/", # main path (all the rest is relative to this path)
             "userinputdir":"userinput/", # user input (includes ACS and HSC image)
             "workdir":"work/", # working directory
             "sexinputdir":"sextractor_input/", # SExtractor input directory (includes config files)
             "out_prefix":"SB4", # output prefix (directory will be created in "outputdir")
             
             "cutoutsize_arcsec":[3*60,3*60], # cutout size of tiles in arcseconds (length of side of box)
             "cutoutoverlap_arcsec":[10,10], # overlap region (is added to the cutoutsize_arcsec)
             "cutout_output_dir": "/stage/irsa-jointproc-data03/TRACTOR/run_acs_jsp_full_SB4_Jan20/SANDBOX4/cutouts/", # where to save the cutouts (absolute path)             
             
             "hr_wht_image_type":"MAP_WEIGHT", # in SExtractor lingo (variance image, weight image, etc)
             
             "hr_compl_niter":15, # number of tractor iterations for complex models (high-resolution image)
             "hr_ps_niter":5, # number of tractor interations for pointsource fitting (high-resolution image)
             "lr_compl_niter":10, # number of tractor iterations for complex models (low-resolution image)
             "hr_zp":25.94734, # photometric zeropoint of high-resolution image
             "lr_zp":27,
             
             "run_on_mock": False, # if true, it uses a mock LR image created with "create_mock_image_from_HR_image()" (NAME_OF_LR_IMAGE_MOCK.fits)

             "sex_command":"/usr/bin/sextractor" # SExtractor command (changes depending on system)
             #"sex_command":"/Users/afaisst/Work/Tools/SExtractor/sextractor-2.19.5/bin/sex"
            }'''
template_userinput = {
             "workdir":"../work/", # working directory
             "sexinputdir":"../sextractor_input/", # SExtractor input directory (includes config files)
             "out_prefix":"SB4", # output prefix (directory will be created in "outputdir")
             "cutoutsize_arcsec":[3*60,3*60], # cutout size of tiles in arcseconds (length of side of box)
             "cutoutoverlap_arcsec":[10,10], # overlap region (is added to the cutoutsize_arcsec)
             "cutout_output_dir": "../cutouts/", # where to save the cutouts (absolute path)                          
             "hr_compl_niter":15, # number of tractor iterations for complex models (high-resolution image)
             "lr_compl_niter":10, # number of tractor iterations for complex models (low-resolution image)
             "hr_zp":25.94734, # photometric zeropoint of high-resolution image
             "lr_zp":27, # photometric zeropoint of low-resolution image
             "sex_command":"/usr/bin/sextractor" # SExtractor command (changes depending on system)
			 "hr_psf_type":"fits", # PSF type for high-resolution image ("fwhm", "fits", or "psfex")
			 "lr_psf_type":"psfex", # PSF type for low-resolution image ("fwhm", "fits", or "psfex")
			 "hr_large_image_name":"tbd", # high-resolution large image (to be changed)
			 "lr_large_image_name":"tbd", # low-resolution large image (to be changed)
			 "hr_image_psf":"tbd", # high-resolution PSF (to be changed)
			 "lr_image_psf":"tbd", # low-resolution PSF (to be changed)
			 "hr_astrometry_correction_name":"tbd", # astrometry offsets HR - Gaia (to be changed)
			 "lr_astrometry_correction_name":"tbd", # astrometry offsets LR - Gaia (to be changed)
			 "tile_id":0 # tile ID (to be changed)
            }


## 3. Other input
patch_size = [12,12] # in arc minutes (needed to compute maximal number of tiles)


####################


##### RUN ########

## get a list of files

# all HR images
hr_sci = glob.glob("%s/*_sci.fits" % directories["hr_dir"])
#print(hr_sci)

# all LR image
lr_img = glob.glob("%s/*.fits" % directories["lr_dir"])
#print(lr_img)

# check
if len(hr_sci) == len(lr_img):
	print("Number of patches found: %g" % len(hr_sci))
else:
	print("Something is wrong: not the same number of LR and HR files. Exit program!")
	quit()


# expected number of tiles.
number_of_tiles = (patch_size[0]*patch_size[1]) / (template_userinput["cutoutsize_arcsec"][0]/60*template_userinput["cutoutsize_arcsec"][1]/60)
print("Total number of tiles: %g" % number_of_tiles ) 


# go through each patch, collect astrometry and PSFs then write job file for each tile.
running_job_number = 1
for ii in range(len(hr_sci)):

	# start with lr image
	lr_large_image_name = lr_img[ii]
	lr_name = lr_large_image_name.split("/")[-1].split(".fits")[0]
	
	# get hr image
	hr_large_image_name = glob.glob("%s/%s*_sci.fits" % (directories["hr_dir"],lr_name) )
	if len(hr_large_image_name) > 1:
		print("HR image name is not unique! - Abort!")
		quit()
	else:
		hr_large_image_name = hr_large_image_name[0]



	# get PSF
	#hr_psf = "f814w_flat_tinytim_psf_60mas.fits"
	#hr_psf = "psf_F814W_undistorted_repixto0.03_center.fits"
	hr_psf = "HSC-I-9813-5_4-2812_psf.fits"
	hr_image_psf = os.path.join(directories["hr_psf_dir"],hr_psf)

	lr_psf = glob.glob("%s/*%s*.psf" % (directories["lr_psf_dir"],lr_name))
	if len(lr_psf) > 1:
		print("More than 1 PSF found. Abort!")
		quit()
	elif len(lr_psf) == 0:
		print("No PSF found. Abort!")
		quit()
	elif len(lr_psf) == 1:
		lr_image_psf = lr_psf[0]
		#print(lr_image_psf)

	# get astrometry (for LR image)
	lr_astro = glob.glob("%s/rdelt-%s*.txt" % (directories["lr_astrometry_dir"],lr_name))
	if len(lr_astro) > 1:
		print("More than 1 astrometry file found. Abort!")
		quit()
	elif len(lr_astro) == 0:
		print("No LR astrometry file found. Abort!")
		quit()
	elif len(lr_astro) == 1:
		lr_astrometry_correction_name = lr_astro[0]
		#print(lr_astrometry_correction_name)

	# get astrometry (for HR image)
	hr_astro = glob.glob("%s/delta-%s*.txt" % (directories["hr_astrometry_dir"],lr_name.replace("calexp-","").replace("_",","))) # the last replacement is a bit risky but works here because no other "_" in the name.
	if len(hr_astro) > 1:
		print("More than 1 astrometry file found. Abort!")
		quit()
	elif len(hr_astro) == 0:
		print("No HR astrometry file found. Abort!")
		quit()
	elif len(lr_astro) == 1:
		hr_astrometry_correction_name = hr_astro[0]
		#print(hr_astrometry_correction_name)


	# put everything together (for each tile)
	for tt in range(int(number_of_tiles)): # CHANGE BACK
	#for tt in range(1):

		userinput_this = copy.deepcopy(template_userinput)
		userinput_this["hr_large_image_name"] = hr_large_image_name
		userinput_this["lr_large_image_name"] = lr_large_image_name	
		userinput_this["hr_image_psf"] = hr_image_psf
		userinput_this["lr_image_psf"] = lr_image_psf
		userinput_this["lr_astrometry_correction_name"] = lr_astrometry_correction_name
		userinput_this["hr_astrometry_correction_name"] = hr_astrometry_correction_name
		userinput_this["tile_id"] =  tt ## CHANGE BACK

		# and save
		with open( os.path.join(directories["job_output_dir"],"job_%g.json" % running_job_number), 'w') as fp:
			json.dump(userinput_this, fp)

		# advance running job number
		running_job_number += 1 # CHANGE BACK






## also add tile ID