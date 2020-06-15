#### CREATE FAKE IMAGE FROM HR MODEL ######
def create_mock_image_from_HR_image(userinput, tileid , astrometric_offset_sigma_mas):

    this_tile_id = tileid
    #this_tile_id = 3
    #astrometric_offset_sigma_mas = [60,60] # in mas



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
    tmp = clip(lr_img, n=3, niter=10)
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


    ## 5. Get Astrometry offset in RA and DEC in pixels ==========
    # Be careful here with direction!!!!
    astro_offset = load_astrometry_correction_table(userinput["lr_astrometry_correction_name"])
    hsc_gaia_median_ra = np.nanmedian(astro_offset["delta_ra_mas"]) # in mas
    hsc_gaia_median_dec = np.nanmedian(astro_offset["delta_dec_mas"]) # in mas 
    acs_gaia_median_ra = 15.8
    acs_gaia_median_dec = 2.3
    hsc_acs_median_ra = hsc_gaia_median_ra - acs_gaia_median_ra
    hsc_acs_median_dec = hsc_gaia_median_dec - acs_gaia_median_dec

    delta_X_median_px = (-1)*hsc_acs_median_dec/1000 / lr_pixscale
    delta_Y_median_px = hsc_acs_median_ra/1000 / lr_pixscale

    print("Astrometric offset in mas (ra,dec) = (%5.3g,%5.3g)" % (hsc_acs_median_ra,hsc_acs_median_dec))
    print("Astrometric offset in pixels (X,Y) = (%5.3g,%5.3g)" % (delta_X_median_px,delta_Y_median_px))



    ## 5. Load PSF (LR image) ==============
    if is_number(userinput["lr_image_psf"]):
        PSF_FWHM_ARCSEC_LR = userinput["lr_image_psf"]
        PSF_LR_TYPE = "fwhm"
        print("Using given FWHM and gaussian for PSF of low-res image.")
    else:
        PSF_LR_TYPE = "pixel"
        with fits.open(os.path.join(userinput["lr_image_psf"])) as hdul:
            PSF_LR_PIXEL = modelPSF_1D(psfcube=hdul[1].data,psfcube_h=hdul[1].header,param=22)
            PSF_LR_PIXEL = PSF_LR_PIXEL / np.nansum(PSF_LR_PIXEL)
        print("Using pixel PSF for low-res image.")



    ## 6. Select galaxies ===========

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
        src_model_hr_to_lr.clear()
    except Exception as e:
        pass
    src_model_hr_to_lr = []
    for ii in range(len(cat_models_hr_use)):

        tmp = lr_img_wcs.all_world2pix([[cat_models_hr_use["RA_tractor.hr"][ii],cat_models_hr_use["DEC_tractor.hr"][ii]]],0)
        pos_init = PixPos(tmp[0][0],tmp[0][1]) # apply reversed astrometric offset??
        x_random_offset = 1e9
        y_random_offset = 1e9
        while np.abs(x_random_offset) > 2*astrometric_offset_sigma_mas[1]:
            x_random_offset = np.random.normal(loc=0,scale=astrometric_offset_sigma_mas[1]*1e-3 / lr_pixscale , size=1)
        while np.abs(y_random_offset) > 2*astrometric_offset_sigma_mas[0]:
            y_random_offset = np.random.normal(loc=0,scale=astrometric_offset_sigma_mas[0]*1e-3 / lr_pixscale , size=1)
        pos_init = PixPos(tmp[0][0]+x_random_offset + delta_X_median_px,
                          tmp[0][1]+y_random_offset + delta_Y_median_px)
        #pos_init = PixPos(tmp[0][0] + delta_X_median_px , tmp[0][1] + delta_Y_median_px) # apply reversed astrometric offset
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

    hr_to_lr_model = tractor_model.getModelImage(0)
    hr_to_lr_img = hr_to_lr_model + np.random.normal(loc=0 , scale=1*lr_img_pixnoise , size=hr_to_lr_model.shape)
    #hr_res = hr_img - hr_mod


    ## 7. Save mock image ========
    h1 = lr_img_h.copy()

    hdu_0 = fits.PrimaryHDU()
    hdu_0.data = hr_to_lr_img.copy()
    hdu_0.header = lr_img_h.copy()
    hdu_0.header["EXTNAME"] = "IMAGE"

    hdu_1 = fits.ImageHDU(lr_mask.copy())
    hdu_1.header = lr_mask_h.copy()
    hdu_1.header["EXTNAME"] = "MASK"

    hdu_0.verify('silentfix')
    hdu_1.verify('silentfix')

    hdul_new = fits.HDUList([hdu_0,hdu_1]) # IMAGE , MASK

    hdul_new.verify("silentfix") # fix everything that is broken
    mock_output_path = lr_image_path.replace(".fits","_MOCK.fits")
    hdul_new.writeto(mock_output_path, overwrite=True)

    print("Mock image created, it's here: %s!" % (mock_output_path) )
    
    return(True)