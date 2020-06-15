##### CREATE SEGMENTATION MAP FOR LR IMAGE #####
# This will then be used to run the forced photometry on an "object by object" basis
# Steps:
# 1) create model image on LR scale
# 2) convolve model image with simple gaussian kernel
# 3) find connected area

def create_LR_segmap_from_HR_model(userinput , tileid):

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
    tmp = clip(lr_img, n=5, niter=10)
    lr_img_pixnoise = tmp["stdev"]
    lr_img_medbkg = tmp["med"]
    print("Pixelnoise of low-res image: %g" % lr_img_pixnoise )
    print("Median background of low-res image: %g" % lr_img_medbkg )


    ## 4. Read High-resolution table with model fits ===========
    print("Reading high-resolution model table . . .",end="")
    start_time = time.time()

    outname = dir_this_process + "/" + "hr_compl_table_final.fits"
    if not os.path.exists(outname):
        print("The models for this high-resolution image (%s, tile %g) do not exist, yet. Please run get_HR_model() first." % (userinput["lr_large_image_name"],this_tile_id) )
        #continue # UNCOMMENT THIS ONCE EVERYTHING IS IN A LOOP
    cat_models_hr = Table.read(outname)

    print(" done (in %g seconds)" % (round((time.time()-start_time),2)) )
    TIMESUMMARY[len(TIMESUMMARY)+1] = ["Reading high-resolution model table",round((time.time()-start_time),2)]


    
    ## 5. Select galaxies ===========

    # add local x/y to the table
    cat_models_hr["X_IMAGE_hr.lr"] = np.zeros(len(cat_models_hr))
    cat_models_hr["Y_IMAGE_hr.lr"] = np.zeros(len(cat_models_hr))
    for ii in range(len(cat_models_hr)):
        tmp = lr_img_wcs.all_world2pix([[cat_models_hr["RA_tractor.hr"][ii],cat_models_hr["DEC_tractor.hr"][ii]]],1)
        cat_models_hr["X_IMAGE_hr.lr"][ii] = tmp[0][0]
        cat_models_hr["Y_IMAGE_hr.lr"][ii] = tmp[0][1]

    # select galaxies on image
    sel_good = np.where( (cat_models_hr["X_IMAGE_hr.lr"]-1 >= 0)
                       & (cat_models_hr["X_IMAGE_hr.lr"]-1 < lr_img.shape[1])
                       & (cat_models_hr["Y_IMAGE_hr.lr"]-1 >= 0)
                       & (cat_models_hr["Y_IMAGE_hr.lr"]-1 < lr_img.shape[0])
                       )[0]
    cat_models_hr_use = cat_models_hr[sel_good].copy()

    
    
    ## 6. Create model image on LR scale =================

    # create Tractor image
    tim_model_lr = Image(data=np.zeros(lr_img.shape), invvar=np.ones_like(lr_img) / (lr_img_pixnoise**2),
                psf=NCircularGaussianPSF([0.3/2.35/lr_pixscale], [1.]),
                wcs=NullWCS(), photocal=NullPhotoCal(),
                sky=ConstantSky(0.))

    # add sources
    try:
        src_model_hr_to_lr.clear()
    except Exception as e:
        pass
    src_model_hr_to_lr = []
    for ii in range(len(cat_models_hr_use)):

        tmp = lr_img_wcs.all_world2pix([[cat_models_hr_use["RA_tractor.hr"][ii],cat_models_hr_use["DEC_tractor.hr"][ii]]],0)
        pos_init = PixPos(tmp[0][0],tmp[0][1]) # apply reversed astrometric offset??    
        flux_init = Flux(cat_models_hr_use["brightness.Flux.hr"][ii] * 10**(0.4*( userinput["lr_zp"] - userinput["hr_zp"] ))) 
        if cat_models_hr_use["is_pointsource"][ii] == 2:
            thissource = PointSource(pos_init,flux_init)
        else:
            n_init = SersicIndex(cat_models_hr_use["sersicindex.SersicIndex.hr"][ii])
            ba_init = cat_models_hr_use["shape.ab.hr"][ii]
            re_init = cat_models_hr_use["shape.re.hr"][ii] * (hr_pixscale / lr_pixscale)
            phi_init = cat_models_hr_use["shape.phi.hr"][ii]
            thissource = SersicGalaxy(pos_init,
                                      flux_init,
                                      GalaxyShape(re_init, ba_init, phi_init), n_init)

        src_model_hr_to_lr.append(thissource)


    # create tractor object
    tractor_model = Tractor([tim_model_lr] , src_model_hr_to_lr)
        
    hr_to_lr_model = tractor_model.getModelImage(0) + lr_img_medbkg
    #hr_res = hr_img - hr_mod


    ## 7. Find connected area ================

    # first pass: create segmentation map for hr objects
    img_labeled1, nr_objects1 = ndimage.label(hr_to_lr_model > lr_pixscale/3) 
    img_labeled1_binary = img_labeled1.copy()
    img_labeled1_binary[img_labeled1 >= 1] = 1

    # second pass: smooth that segmentation map to adjust to GB PSF
    imgf = ndimage.maximum_filter(img_labeled1_binary, size=1.0 / lr_pixscale)
    img_labeled, nr_objects = ndimage.label(imgf > lr_pixscale/3) 
    img_labeled_binary = img_labeled.copy()
    img_labeled_binary[img_labeled >= 1] = 1

    print("Number of objects is {}".format(nr_objects))


    ## 8. find the extension of the labeled objects =================
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
    ascii.write(EXT_TABLE,output=os.path.join(dir_this_process,"lr_segmap_extensions.csv") , format="csv", overwrite=True)


    ## 9. save the "segmentation map" =======================
    h1 = lr_img_h.copy()

    hdu_1 = fits.PrimaryHDU()
    hdu_1.data = img_labeled.copy()
    hdu_1.header = h1.copy()
    hdu_1.header["EXTNAME"] = "SEGMAP"

    hdul_new = fits.HDUList([hdu_1]) # segmap
    hdul_new.verify("silentfix") # fix everything that is broken
    hdul_new.writeto(os.path.join(dir_this_process,"lr_segmap_from_hr.fits"), overwrite=True)



    ## 10. plot them =======================
    fig = plt.figure(figsize=(30,10))
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.01, hspace=0.0)

    ima = dict(cmap=plt.get_cmap("Greys"),vmin=-2*lr_img_pixnoise,vmax=10*lr_img_pixnoise,origin="lower",interpolation="nearest")
    ima_labeled = dict(cmap=plt.get_cmap("Greys"),vmin=0,vmax=1,origin="lower",interpolation="nearest")

    # models
    ax1 = fig.add_subplot(1,3,1)
    im1 = ax1.imshow(hr_to_lr_model,**ima)
    
    for ii in range(len(EXT_TABLE)):
        ax1.plot([EXT_TABLE["XMIN"][ii],EXT_TABLE["XMIN"][ii],EXT_TABLE["XMAX"][ii],EXT_TABLE["XMAX"][ii],EXT_TABLE["XMIN"][ii]],
                 [EXT_TABLE["YMIN"][ii],EXT_TABLE["YMAX"][ii],EXT_TABLE["YMAX"][ii],EXT_TABLE["YMIN"][ii],EXT_TABLE["YMIN"][ii]],
                 "--",
                  color="blue",
                  linewidth=0.5,
                dashes=(5,5))
    
    fig.colorbar(im1,shrink=0.5)
    ax1.set_xlim(0,lr_img.shape[1])
    ax1.set_ylim(0,lr_img.shape[0])
    ax1.set_title("Model from HR image",fontsize=13)

    # segmentation map
    ax2 = fig.add_subplot(1,3,2)
    im2 = ax2.imshow(img_labeled, **ima_labeled)
    fig.colorbar(im2,shrink=0.5)
    ax2.set_xlim(0,lr_img.shape[1])
    ax2.set_ylim(0,lr_img.shape[0])
    ax2.set_title("Segmentation map produced from HR image",fontsize=13) 

    # original LR with segmap
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
    outfile = os.path.join(dir_this_process,"lr_segmap.pdf")
    SHOWPLOT = userinput["show_plots"]
    plt.savefig(outfile,bbox_inches='tight')
    if SHOWPLOT:
        plt.show()
    else:
        plt.close()
        
    return(True)
