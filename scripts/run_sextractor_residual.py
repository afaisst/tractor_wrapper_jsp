## RUN SEXTRACTOR (AND TRACTOR) ON RESIDUAL IMAGE #####
def run_sextractor_residual(userinput,tileid,doplot):

    # 1. Mask the detected/fitted sources.

    ## 1. Create process ID ===========================
    this_tile_id = tileid
    #this_tile_id = 3

    TIMESUMMARY = dict()

    process_id = "%04.0f_%s" % (this_tile_id , userinput["lr_large_image_name"].split("/")[-1].split(".fits")[0] ) # set the process ID to the name of the large image array + the tile number
    print("This Process is: " + str(process_id))


    ## 2. create (temporary) working directory ===========================
    dir_this_process = os.path.join(userinput["workdir"], "%s_%s" % (userinput["out_prefix"],process_id) )
    if not os.path.exists(dir_this_process):
        sh.mkdir('-p', dir_this_process)
    else:
        print("Directory %s already exists" % dir_this_process)


    ## 3. Load images ===========
    print("reading images . . .",end="")
    start_time = time.time()
    lr_results_image_path = os.path.join(dir_this_process,"lr_tractor_results.fits")
    with fits.open(lr_results_image_path) as hdul:
        lr_img = hdul["ORIGINAL"].data
        lr_img_h = hdul["ORIGINAL"].header
        lr_model = hdul["COMPL_MODEL"].data
        lr_res = hdul["COMPL_RES"].data
        lr_seg = hdul["SEGMAP_SEXTRACTOR"].data
        lr_seg_broad = hdul["SEGMAP"].data
        lr_mask = hdul["MASK"].data    
        #hdul.info()

    lr_pixscale = np.abs(lr_img_h["CD1_1"]*3600) # pixel scale in arcsec/px
    lr_img_wcs = wcs.WCS(lr_img_h)

    print(" done (in %g seconds)" % (round((time.time()-start_time),2)) )
    TIMESUMMARY[0] = ["Reading images",round((time.time()-start_time),2)]


    ## 3.5 Get LR image properties ============

    ## Check some stuff
    if "CD1_2" not in lr_img_h.keys():
        lr_img_h["CD1_2"] = 0
    if "CD2_1" not in lr_img_h.keys():
        lr_img_h["CD2_1"] = 0

    ## Get Pixel noise
    tmp = clip(lr_img[lr_seg_broad == 0], n=3, niter=15)
    lr_img_pixnoise = tmp["stdev"]
    lr_img_medbkg = tmp["med"]
    print("Pixelnoise of low-res image: %g" % lr_img_pixnoise )
    print("Median background of low-res image: %g" % lr_img_medbkg )


    ## 4. Run SExtractor =====================

    ## Copy SExtractor config file
    config_default = os.path.join(userinput["sexinputdir"],"default.conf.interactive")
    config_this_process = os.path.join(dir_this_process,"lr_residual_default.conf")
    cmd = "cp " + config_default + " " + config_this_process
    subprocess.run(cmd, shell=True)


    ## Adjust configuration file
    sex_output_cat_lr_file = os.path.join(dir_this_process,"lr_residual_sex_output.cat")
    PARS_this_process = get_default_parfile()
    PARS_this_process["CATALOG_NAME"] = sex_output_cat_lr_file
    PARS_this_process["PARAMETERS_NAME"] = os.path.join(userinput["sexinputdir"],"sex.par")
    PARS_this_process["FILTER_NAME"] = os.path.join(userinput["sexinputdir"],"g2.8.conv")
    PARS_this_process["STARNNW_NAME"] = os.path.join(userinput["sexinputdir"],"default.nnw")
    PARS_this_process["CHECKIMAGE_NAME"] = "%s" % (os.path.join(dir_this_process,"lr_residual_seg.fits") )
    PARS_this_process["CHECKIMAGE_TYPE"] = "SEGMENTATION"
    PARS_this_process["PHOT_APERTURES"] = ','.join(map(str, get_apertures(lr_pixscale,aperturelist_arcsec=[3.0]).tolist())) # Note: if number of apertures are changes, also change sex.par file!!!
    PARS_this_process["PIXEL_SCALE"] = lr_pixscale
    PARS_this_process["DEBLEND_MINCONT"] = 0.01
    PARS_this_process["DETECT_MINAREA"] = 3
    PARS_this_process["DETECT_THRESH"] = 1.5#2
    PARS_this_process["MAG_ZEROPOINT"] = userinput["lr_zp"]
    PARS_this_process["ANALYSIS_THRESH"] = 1.5#1.5
    PARS_this_process["SEEING_FWHM"] = 0.6 # Doesn't really matter for SExtractor

    for key in PARS_this_process.keys():
        replace_in_file(filename = config_this_process,
                       old_string = "*" + key + "*",
                       new_string = str(PARS_this_process[key]))


    ## Run SExtractor
    start_time = time.time()
    print("running SExtractor on residual image . . .",end="")
    cmd = "%s %s[2] -c %s" % (userinput["sex_command"],
                           lr_results_image_path,  # 0 = original, 1=model, 2=residual
                           config_this_process
                          )
    subprocess.run(cmd , shell=True)
    #print(cmd)
    print(" done (in %g minutes)" % (round((time.time()-start_time)/60,2)) )
    TIMESUMMARY[len(TIMESUMMARY)+1] = ["Running SExtractor on residual image",round((time.time()-start_time),2)]


    ## Load catalog
    sexcat_lr_residual = ascii.read(sex_output_cat_lr_file)
    print("%g objects detected!" % len(sexcat_lr_residual))


    ## Load segmentation map
    with fits.open(os.path.join(dir_this_process,"lr_residual_seg.fits") ) as hdul:
        lr_residual_seg = hdul[0].data # 0 is segmentation


    
    ## 5. Clean catalog (remove detections that are already fitted) ==========================

    ## first have to create a masks so we know where we fitted objects
    mask_fitted = lr_model.copy()
    mask_fitted[lr_model > 3*lr_img_pixnoise] = 1
    mask_fitted[lr_model <= 3*lr_img_pixnoise] = 0
    mask_fitted_broad = ndimage.maximum_filter(mask_fitted, size=1.0 / lr_pixscale)

    ## remove objects that are fitted already (including those who are close to the residuals, we don't want junk.)
    ## Also change the segmentation map (remove areas of galaxies that were fitted already or are close to residuals)
    sexcat_lr_residual["flag_fitted"] = np.zeros(len(sexcat_lr_residual))
    lr_residual_seg_final = lr_residual_seg.copy()
    for ii in range(len(sexcat_lr_residual)):
        if (mask_fitted_broad[ int(sexcat_lr_residual["Y_IMAGE"][ii]-1) , int(sexcat_lr_residual["X_IMAGE"][ii]-1) ] != 0) & (lr_mask[ int(sexcat_lr_residual["Y_IMAGE"][ii]-1) , int(sexcat_lr_residual["X_IMAGE"][ii]-1) ] != 0):
            sexcat_lr_residual["flag_fitted"][ii] = 1
            lr_residual_seg_final[ lr_residual_seg == sexcat_lr_residual["NUMBER"][ii] ] = 0

    sexcat_lr_residual_final = sexcat_lr_residual.copy()[ sexcat_lr_residual["flag_fitted"] == 0 ]
    

    ## Save catalog
    sexcat_lr_residual_final.write(os.path.join(dir_this_process,"lr_residual_catalog.fits") , format="fits", overwrite=True)
    print("%g new objects!" % len(sexcat_lr_residual_final))

    ## save final segmentation map
    h1 = lr_img_h.copy()

    hdu_1 = fits.PrimaryHDU()
    hdu_1.data = lr_residual_seg_final.copy()
    hdu_1.header = h1.copy()
    hdu_1.header["EXTNAME"] = "RESIDUAL_SEGMAP"
    hdul_new = fits.HDUList([hdu_1]) # segmap
    hdul_new.verify("silentfix") # fix everything that is broken
    hdul_new.writeto(os.path.join(dir_this_process,"lr_residual_segmap_clean.fits"), overwrite=True)




    # Create DS9 file with extractions
    lr_residual_ds9radec_name = os.path.join(dir_this_process,"lr_residual_radec.reg")
    with open(lr_residual_ds9radec_name,"w+") as f:
        f.write('global color=magenta dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        f.write('fk5\n')
        for jj in range(len(sexcat_lr_residual_final)):
            col_this = "magenta"
            f.write('circle(%s,%s,0.5") # color=%s width=2 text={%s}\n' % (sexcat_lr_residual_final["ALPHA_J2000"][jj] , sexcat_lr_residual_final["DELTA_J2000"][jj], col_this, sexcat_lr_residual_final["NUMBER"][jj]) )
    
    
    
    
    ## 6. Make a plot ======================
    if doplot:
	    ima = dict(cmap=plt.get_cmap("Greys"),vmin=-2*lr_img_pixnoise,vmax=10*lr_img_pixnoise,origin="lower",interpolation="nearest")
	    ima_sqrt = dict(cmap=plt.get_cmap("Greys"),vmin=-2*np.sqrt(lr_img_pixnoise),vmax=10*np.sqrt(lr_img_pixnoise),origin="lower",interpolation="nearest")


	    fig = plt.figure(figsize=(10,10))
	    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)

	    ## Original image
	    ax1 = fig.add_subplot(1,1,1,projection = lr_img_wcs)
	    ax1.imshow(lr_img,**ima)
	    
	    ax1.plot(sexcat_lr_residual_final["X_IMAGE"]-1,sexcat_lr_residual_final["Y_IMAGE"]-1,"s",markersize=10,markeredgewidth=1,fillstyle="none",alpha=1,color=def_cols[0])
	    for iii in range(len(sexcat_lr_residual_final)):
	        text = ax1.text(sexcat_lr_residual_final["ALPHA_J2000"][iii],sexcat_lr_residual_final["DELTA_J2000"][iii]+0.7/3600,sexcat_lr_residual_final["NUMBER"][iii],fontsize=11,color="black",va="bottom",ha="left",transform=ax1.get_transform('world'))
	        text.set_path_effects([path_effects.Stroke(linewidth=1.0, foreground='white'),
	                           path_effects.Normal()])

	    ax1.set_title("Original low-resolution image with residual sources",fontsize=15)
	    ax1.set_xlabel("Right-ascension (J2000)",fontsize=12,labelpad=2)
	    ax1.set_ylabel("Declination (J2000)",fontsize=12,labelpad=2)
	    ax1.xaxis.set_tick_params(labelsize=13)
	    ax1.yaxis.set_tick_params(labelsize=13)
	    ax1.tick_params(axis="both",which="minor",length=2)
	    ax1.tick_params(axis="both",which="major",length=4)
	    ax1.coords[0].set_major_formatter('hh:mm:ss.s')
	    ax1.coords[1].set_major_formatter('dd:mm:ss.s')
	    ax1.coords[0].set_ticklabel(size=11)
	    ax1.coords[1].set_ticklabel(size=11)
	    ax1.minorticks_on()





	    outfile = os.path.join(dir_this_process,"lr_residual_fig.pdf")
	    SHOWPLOT = userinput["show_plots"]
	    plt.savefig(outfile,bbox_inches='tight')
	    if SHOWPLOT:
	        plt.show()
	    else:
	        plt.close()
    
    
    
    
    return(True)
