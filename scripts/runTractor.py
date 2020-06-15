## Run as runTractor.py jobX.json
#
# A. Faisst
#
###############################################

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

# import ConstrainedOptimizer()
import_file_list = ["constrained_optimizer.py",
                   "functions.py",
                   "run_sextractor_residual.py",
                   "get_HR_model.py",
                   "fit_LR_image.py",
                   "create_LR_segmap_from_SExtractor.py",
                   "create_LR_segmap_from_HR_model.py",
                   "create_cutouts.py",
                   "plot_lr_images.py",
                   "plot_hr_images.py",
                   "run_MP_functions.py"]
for file in import_file_list:
    exec(compile(open(file, "rb").read(), file, 'exec'))


## Plotting stuff
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
#mpl.rc('text', usetex=True)
mpl.rc('font', family='serif')
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['hatch.linewidth'] = 4

def_cols = plt.rcParams['axes.prop_cycle'].by_key()['color']




########### RUN ###########
LOG = []

## Read user input file
userinput_file = sys.argv[1]
with open(userinput_file) as json_file:
    userinput = json.load(json_file)


## get tile ID
#tile_id = int(sys.argv[1])
tile_id = userinput["tile_id"]
print("tile ID: %g" % tile_id)

## Get file name for log file
#output_logfile_name = sys.argv[2]
output_logfile_name = "%s.txt" % ".".join(userinput_file.split("/")[-1].split(".")[0:-1])
userinput["output_logfile_name"] = output_logfile_name # add just in case it will be needed in the future
print("output log file name: %s" % output_logfile_name)



## Create cutouts
cutout_results = create_cutouts(userinput,overwrite_cutouts=True,tile_nbr=tile_id) # generates N cutouts (starting at 0)
print(cutout_results)


if cutout_results["success"]:
	## Run Tractor
	start_time = time.time() # get time
	runTractor_helper(userinput=userinput , usemultiproc=False , tile_id_list=tile_id , doplots=False)
	elapsed_time = time.time() - start_time
	print("\n\n ALL FINISHED IN: " + str(round(elapsed_time/60,2)) + " minutes")
	LOG.append("ALL FINISHED IN: " + str(round(elapsed_time/60,2)) + " minutes")

	## Make plots
	print("MAKING PLOTS FOR TILE ID %g" % tile_id)
	LOG.append("MAKING PLOTS FOR TILE ID %g" % tile_id)

	plot_hr_images(userinput=userinput , tileid=tile_id , plot_ids=False)
	plot_lr_images(userinput=userinput , tileid=tile_id , plot_ids=False)

	print("ALL DONE!")
	LOG.append("ALL DONE!")
else:
	LOG.append("Failed. Invalid tile ID.")

## save log file
with open(output_logfile_name,"w") as f:
	for ii in range(len(LOG)):
		f.write( "%s\n" % (LOG[ii]) )


