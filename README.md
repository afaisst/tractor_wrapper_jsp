# Tractor Wrapper for Joint Survey Processing

### Introduction

Wrapper for Tractor to run it on ACS and HSC images (or other combinations, however, not tested). The code is structured so it can be run on computer clusters such as at NERSC. 


### Usage

#### General Setup

It is very straight forward to run the code. All the information about the input data is captured in job-files, which are dictionaries in JSON format. These job files are discussed in a section below.
To run a job file, simply run the following line in the "scripts/" directory:
```
python runTractor.py JOBFILE
```
where "JOBFILE" is the path to the job file that you want to run.

#### Running things in parallel

Several job files can be run in parallel easily using the GNU parallel package.
The script "run.sh" does exactly this. The script contains only one line (and a bunch of explanations):

```
#!/bin/bash
# Run all N jobs on C CPUs (2*C vCPUs)
# Do not forget to time it
# parallel -j C ./task_runner.sh {%} {} ::: {1..N}

parallel -j 1 ./task_runner.sh {%} {} ::: {1..1}
```

The scripts runs a series of job files (called job_X.json, where X = 1,2, ... , N and N the total number of job files). Each of these jobs is run on C CPUs (or 2 * C virtual CPUs). The example above would run 1 job on 1 CPU. Alternatively, you can say 
```
parallel -j 2 ./task_runner.sh {%} {} ::: {1..4}
```
which would run 4 jobs on 2 CPUs each. If you have 4 CPUs in total, this would run two jobs at the same time in parallel.

The task runner file (task_runner.sh) looks something like this:
```
#! /bin/bash
# task_runner.sh
 
slot=$1
 
jobnbr=$2
 
# This function determines which virtual CPUs correspond to the specified slot:
cpu_threads() {
    minus=$((slot - 1))
    first_cpu_thread=$((2 * minus))
    second_cpu_thread=$((first_cpu_thread + 1))
    echo $first_cpu_thread,$second_cpu_thread
}

taskset -c $(cpu_threads) python runTractor.py ../example_data/jobs/job_${jobnbr}.json
```


Note that there are other ways to handle parallel computing. These can be easily implemented given the very simple syntax to run the main code (python runTractor.py JOBFILE). One way is to directly feed a list of JSON files directly to "parallel", for example:

```
ls my/path/*.json | parallel -j C ./task_runner_list.sh {%} {}
```
Note that the file "task_runner_list.sh" is not the same as above ("task_runner_jsp.sh"). It needs to be slightly modified to accept a list instead of a number:

```
#! /bin/bash
# task_runner_list.sh
 
slot=$1
 
jobfile=$2
 
# This function determines which virtual CPUs correspond to the specified slot:
cpu_threads() {
    minus=$((slot - 1))
    first_cpu_thread=$((2 * minus))
    second_cpu_thread=$((first_cpu_thread + 1))
    echo $first_cpu_thread,$second_cpu_thread
}

taskset -c $(cpu_threads) python runTractor.py ${jobfile}
```


#### Things to addjust and to create

If you run the wrapper on a single job file, then the only thing you need to change is the job file itself (see below).
If you run the wrapper in parallel, you need to change
- the job file(s) (see below)
- the path to the job files in the "task_runner_jsp.sh" (or "task_runner_list.sh") file


#### The job file

A job file is a JSON file containing all information for the code to run.
It looks something like this (the file names are taken from the example data, see below):
```
{"workdir": "../work/",
"sexinputdir": "../sextractor_input/",
"out_prefix": "S1",
"cutoutsize_arcsec": [60, 60],
"cutoutoverlap_arcsec": [10, 10],
"cutout_output_dir": "../cutouts/",
"hr_compl_niter": 15,
"lr_compl_niter": 10,
"hr_zp": 25.94734,
"lr_zp": 27,
"sex_command": "/usr/bin/sextractor",
"hr_large_image_name": "../example_data/images/3sqarcmin/hr.fits",
"lr_large_image_name": "../example_data/images/3sqarcmin/lr.fits",
"hr_image_psf": "../example_data/psf/HSC-I-9813-5_4-2812_psf.fits",
"hr_psf_type":"fits",
"lr_image_psf": "../example_data/psf/calexp-HSC-I-9813-5_4-9813_2812.psf",
"lr_psf_type":"psfex",
"lr_astrometry_correction_name": "../example_data/astrometry_offsets/lr_astro.txt",
"hr_astrometry_correction_name": "../example_data/astrometry_offsets/hr_astro.txt",
"hr_psf_fwhm_arcsec":0.1,
"tile_id": 0
}
```

with the following options:
- workdir: Path (absolute or relative) to working directory. All the output will be in there.
- sexinputdir: Path (absolute or relative) to SExtractor input directory. It contains the configuration file template as well as convolution filters.
- out_prefix: Output prefix that should be added to the name. The output name is then [prefix]_[image name]_[tile number]
- cutoutsize_arcsec: The size of the cutout in arcseconds (the script tiles the input image in tiles on which Tractor is run)
- cutoutoverlap_arcsec: Overlap in arcseconds of the cutout tiles
- cutout_output_dir: Path (absolute or relative) where to save the cutouts
- hr_compl_niter: Maximal number of iteration to fit a model to a source in the high-resolution image. Fit is aborted (but result still saved and might not be bad!) if the number of iterations exceeds this limit and no good fit (delta < 1e-3) has been found.
- lr_compl_niter: Maximal number of iteration to fit a model to a source in the low-resolution image. Fit is aborted (but result still saved and might not be bad!) if the number of iterations exceeds this limit and no good fit (delta < 1e-3) has been found.
- hr_zp: Photometric zeropoint of high-resolution image
- lr_zp: Photometric zeropoitn of low-resolution image
- sex_command: Command-line command to run SExtractor (depends on the installation of SExtractor)
- hr_large_image_name: Large (un-tiled) name of the high-resolution image. This image gets tiles using the cutoutsize_arcsec keyword).
- lr_large_image_name: Large (un-tiled) name of the low-resolution image. This image gets tiles using the cutoutsize_arcsec keyword). Note that a mask image should be included in one of the fits extensions using the same format as the HSC images (they have an "IMAGE" and "MASK" and "VARIANCE" extension that are used by this code. Make sure to include them if you, for example, simulated data).
- hr_image_psf: PSF of the high-resolution image (see "hr_psf_type").
- hr_psf_type: PSF type, can be "fits" (FITS image), "fwhm" (the Gaussian FWHM in arcsec, float), "psfex" (PSFex gemerated PSF). 
- lr_image_psf: PSF of the low-resolution image (see "lr_psf_type").
- lr_psf_type: PSF type, can be "fits" (FITS image), "fwhm" (the Gaussian FWHM in arcsec, float), "psfex" (PSFex gemerated PSF). 
- lr_astrometry_correction_name: Path (absolute or relative) to an astrometry offset file for the low-resolution image. The file lists offsets in milli-arcseconds between stars on the low-resolution image and a reference catalog (e.g., Gaia). The sense is LRI-RefCat. There can be multiple stars (in this case the median is taken) or just one value for RA and DEC.
- hr_astrometry_correction_name: Path (absolute or relative) to an astrometry offset file for the high-resolution image. The file lists offsets in milli-arcseconds between stars on the high-resolution image and a reference catalog (e.g., Gaia). The sense is HRI-RefCat. There can be multiple stars (in this case the median is taken) or just one value for RA and DEC.
- hr_psf_fwhm_arcsec: give here an approximate PSF FWHM of the high-resolution PSF. This is used to flag point sources, which are then fit with the Tractor PointSource() class. Simulations showed that fitting point sources as point sources increase the accuracy of the photometry.
-tile_id: This integer number defines the tile on which the code is run. This is the essential number linked to the job file. See below for a detailed explanation.


Note 1: you can define relative paths. For example, use "../FOLDER" if script is run in "scripts/". 
Note 2: PSF can be either a FITS file, a number (given the FWHM of a Gaussian PSF), or in PSFex format. For the latter case: so far, the code only accepts a magnitude-dependent PSFex PSF. It evaluates the PSF at a hard-coded value of mag = 21. Yes, this needs to be changed.

#### Tiling and jobs
The code takes an input image ("lr_large_image_name" and "hr_large_image_name") and first creates cutout tiles of that image according the the size defined by the user (using the keyword "cutoutsize_arcsec").
For example, if the image is of size 3' x 3', and the cutout size is chosen to be [60,60] (meaning 1' x 1'), then 9 tiles are created. **Make sure that the cutout size covers the whole image!**. The keyword "tile_id" lets you tell the program on which tile it should run. This defines a job.
Reasonably, the tile (cutout) size should not exceed 3' x 3'. Else Tractor will run for quite a while and running multiple jobs of smaller tiles is beneficial.

As a practical example, if you have a large image with the size 3' x 3' and you choose a cutout size of 1' x 1' (60" x 60"), you need 9 jobs to finish the image.
You can create 9 job files which are exact copies, except you change the "tile_id" keyword so it is 0 in job_1.json, 1 in job_2.json, .... , 8 in job_9.json. **Not that the tile number starts with 0 but the job files start with 1 (at least in the bash shell scripts that are included here)**


#### Example Data

Example data (image, PSF, astrometry offsets, one job file) can be downloaded here:
https://caltech.box.com/s/dekvuoyf1wzpqi5d6isajksfrt5duh9i

Download the .tar file and unzip it in the same directory as the files in this repository.
You should be able to run the code directly by running (in the "scripts/" directory):

```
python runTractor.py ../example_data/jobs/job_1.json
```

The script should finish in between 1 and 2 minutes depending on your computer.

All the above scripts to run things in parallel work directly with this example data. If you have GNU parallel installed, you can run the one job in parallel (which is non sense, but never mind) by executing:
```
sh run.sh
```
in the "scripts/" directory. This runs 1 job with 1 CPU.

The example data is a 3' x 3' simulation of ACS (hr.fits) and HSC (lr.fits) data. The cutout size is chosen such that there are 9 tiles (see section above). The example data includes one job file ("jobs/job_1.json") that runs tile number 0.

Add more jobs files (changing tile_id to 1, 2, 3, 4, 5, 6, 7, or 8) and adjust the "run.sh" file to run all of them in parallel.


