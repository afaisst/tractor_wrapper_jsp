### CREATE SEGMENTATION MAP FROM SEXTRATOR SEGEMENTATION MAP ####
def create_LR_segmap_from_SExtractor(userinput,tileid,doplot):

    ## 1. Create process ID ===========================
    this_tile_id = tileid

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

    ## low resolution
    lr_image_path = os.path.join(userinput["cutout_output_dir"],userinput["lr_large_image_name"].split("/")[-1].split(".fits")[0],"%04.0f_%s.fits" % (this_tile_id,userinput["lr_large_image_name"].split("/")[-1].split(".fits")[0]))
    with fits.open(lr_image_path) as hdul:
        lr_img = hdul[0].data
        lr_img_h = hdul[0].header

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
    tmp = clip(lr_img, n=5, niter=10)
    lr_img_pixnoise = tmp["stdev"]
    lr_img_medbkg = tmp["med"]
    print("Pixelnoise of low-res image: %g" % lr_img_pixnoise )
    print("Median background of low-res image: %g" % lr_img_medbkg )


    ## 4. Run SExtractor =====================

    ## Copy SExtractor config file
    config_default = os.path.join(userinput["sexinputdir"],"default.conf.interactive")
    config_this_process = os.path.join(dir_this_process,"lr_default.conf")
    cmd = "cp " + config_default + " " + config_this_process
    subprocess.run(cmd, shell=True)


    ## Adjust configuration file
    sex_output_cat_lr_file = os.path.join(dir_this_process,"lr_sex_output.cat")
    PARS_this_process = get_default_parfile()
    PARS_this_process["CATALOG_NAME"] = sex_output_cat_lr_file
    PARS_this_process["PARAMETERS_NAME"] = os.path.join(userinput["sexinputdir"],"sex.par")
    PARS_this_process["FILTER_NAME"] = os.path.join(userinput["sexinputdir"],"g2.8.conv")
    PARS_this_process["STARNNW_NAME"] = os.path.join(userinput["sexinputdir"],"default.nnw")
    PARS_this_process["CHECKIMAGE_NAME"] = "%s" % ( os.path.join(dir_this_process,"lr_seg.fits") )
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


    ## Run SExtractor and load catalog
    start_time = time.time()
    print("running SExtractor . . .",end="")
    cmd = "%s %s[0] -c %s" % (userinput["sex_command"],
                           lr_image_path,  # 0 = image , 1 = mask
                           config_this_process
                          )
    subprocess.run(cmd , shell=True)
    print(" done (in %g minutes)" % (round((time.time()-start_time)/60,2)) )
    TIMESUMMARY[len(TIMESUMMARY)+1] = ["Running SExtractor on low-resolution image",round((time.time()-start_time),2)]


    ## Load segmentation map
    with fits.open(os.path.join(dir_this_process,"lr_seg.fits") ) as hdul:
        lr_seg = hdul[0].data # 0 is segmentation




    ## 4. Find connected area ==========


    # create binary map
    lr_seg_binary = lr_seg.copy()
    lr_seg_binary[lr_seg != 0] = 1

    # apply maximum filter
    imgf = ndimage.maximum_filter(lr_seg_binary, size=1.0 / lr_pixscale)

    # label connected structures
    img_labeled, nr_objects = ndimage.label(imgf) 
    img_labeled_binary = img_labeled.copy()
    img_labeled_binary[img_labeled >= 1] = 1

    print("Number of objects is {}".format(nr_objects))



    ## 5. Get extension of labeled connected area ============

    # find the extension of the labeled objects
    # Extinsion in Python pixel format (starting with (0,0)!!!)
    extensions = ndimage.find_objects(img_labeled)
    EXT_TABLE = Table(names=["label","XMIN","XMAX","YMIN","YMAX"] , dtype=["int"]+4*["f"])
    for ii in range(len(extensions)):
        EXT_TABLE.add_row([ii+1,
                           extensions[ii][1].indices(img_labeled.shape[1])[0],
                           extensions[ii][1].indices(img_labeled.shape[1])[1],
                           extensions[ii][0].indices(img_labeled.shape[0])[0],
                           extensions[ii][0].indices(img_labeled.shape[0])[1]
                          ])
    ascii.write(EXT_TABLE,output=os.path.join(dir_this_process,"lr_segmap_extensions_sextractor.csv") , format="csv", overwrite=True)

    ## 6. Save segmentation map ==========

    # save the "segmentation map"
    h1 = lr_img_h.copy()
    hdu_1 = fits.PrimaryHDU()
    hdu_1.data = img_labeled.copy()
    hdu_1.header = h1.copy()
    hdu_1.header["EXTNAME"] = "SEGMAP"

    hdul_new = fits.HDUList([hdu_1]) # segmap
    hdul_new.verify("silentfix") # fix everything that is broken
    hdul_new.writeto(os.path.join(dir_this_process,"lr_segmap_from_sextractor.fits"), overwrite=True)




    ## 7. Plot them ================
    if doplot:
	    fig = plt.figure(figsize=(30,10))
	    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.01, hspace=0.0)

	    ima = dict(cmap=plt.get_cmap("Greys"),vmin=-2*lr_img_pixnoise,vmax=10*lr_img_pixnoise,origin="lower",interpolation="nearest")
	    ima_labeled = dict(cmap=plt.get_cmap("Greys"),vmin=0,vmax=1,origin="lower",interpolation="nearest")

	    ax1 = fig.add_subplot(1,3,1)
	    im1 = ax1.imshow(lr_seg,**ima_labeled)
	    fig.colorbar(im1,shrink=0.5)
	    ax1.set_xlim(0,lr_img.shape[1])
	    ax1.set_ylim(0,lr_img.shape[0])
	    ax1.set_title("Segmentation map produced by SExtractor",fontsize=13)

	    ax2 = fig.add_subplot(1,3,2)
	    im2 = ax2.imshow(img_labeled, **ima_labeled)
	    fig.colorbar(im2,shrink=0.5)
	    ax2.set_xlim(0,lr_img.shape[1])
	    ax2.set_ylim(0,lr_img.shape[0])
	    ax2.set_title("Segmentation map convolved with box-kernel (maximum filter)",fontsize=13) 


	    ax3 = fig.add_subplot(1,3,3)
	    im3 = ax3.imshow(lr_img, **ima)
	    fig.colorbar(im3,shrink=0.5)
	    ax3.contour(img_labeled_binary,colors="red",linewidths=0.3)

	    for ii in range(len(EXT_TABLE)):
	        ax3.plot([EXT_TABLE["XMIN"][ii],EXT_TABLE["XMIN"][ii],EXT_TABLE["XMAX"][ii],EXT_TABLE["XMAX"][ii],EXT_TABLE["XMIN"][ii]],
	                 [EXT_TABLE["YMIN"][ii],EXT_TABLE["YMAX"][ii],EXT_TABLE["YMAX"][ii],EXT_TABLE["YMIN"][ii],EXT_TABLE["YMIN"][ii]],
	                 "--",
	                  color="blue",
	                  linewidth=0.5,
	                dashes=(5,5))

	    ax3.set_xlim(0,lr_img.shape[1])
	    ax3.set_ylim(0,lr_img.shape[0])
	    ax3.set_title("LR image with segmentation",fontsize=13) 


	    # Save
	    outfile = os.path.join(dir_this_process,"lr_segmap_sextractor.pdf")
	    SHOWPLOT = userinput["show_plots"]
	    plt.savefig(outfile,bbox_inches='tight')
	    if SHOWPLOT:
	        plt.show()
	    else:
	        plt.close()

    return(True)
