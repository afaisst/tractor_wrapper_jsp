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
# parallel -j C ./task_runner_jsp.sh {%} {} ::: {1..N}

parallel -j 1 ./task_runner_jsp.sh {%} {} ::: {1..1}
```

The scripts runs a series of job files (called job_X.json, where X = 1,2, ... , N and N the total number of job files). Each of these jobs is run on C CPUs (or 2 * C virtual CPUs). The example above would run 1 job on 1 CPU. Alternatively, you can say 
```
parallel -j 2 ./task_runner_jsp.sh {%} {} ::: {1..4}
```
which would run 4 jobs on 2 CPUs each. If you have 4 CPUs in total, this would run two jobs at the same time in parallel.

The task runner file (task_runner_jsp.sh) looks something like this:
```
#! /bin/bash
# task_runner_jsp.sh
 
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
It looks something like this:
```
{"workdir": "../work/",
"sexinputdir": "../sextractor_input/",
"out_prefix": "S1",
"cutoutsize_arcsec": [60, 60],
"cutoutoverlap_arcsec": [10, 10],
"cutout_output_dir": "../cutouts/",
"hr_wht_image_type": "MAP_WEIGHT",
"hr_compl_niter": 15,
"hr_ps_niter": 5,
"lr_compl_niter": 10,
"hr_zp": 25.94734,
"lr_zp": 27,
"sex_command": "/usr/bin/sextractor",
"hr_large_image_name": "../example_data/images/3sqarcmin/hr.fits",
"hr_large_image_wht_name": "none",
"lr_large_image_name": "../example_data/images/3sqarcmin/lr.fits",
"hr_image_psf": "../example_data/psf/HSC-I-9813-5_4-2812_psf.fits",
"lr_image_psf": "../example_data/psf/calexp-HSC-I-9813-5_4-9813_2812.psf",
"lr_astrometry_correction_name": "../example_data/astrometry_offsets/lr_astro.txt",
"hr_astrometry_correction_name": "../example_data/astrometry_offsets/hr_astro.txt",
"tile_id": 0
}
```

### Example
