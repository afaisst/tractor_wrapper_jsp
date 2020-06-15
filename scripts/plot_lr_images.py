### CREATE LR IMAGES FOR CHECKING ####
def plot_lr_images(userinput , tileid , plot_ids):
    
    print("PLOTTING LR IMAGES")
    
    ## 1. Create process ID ===========================
    this_tile_id = tileid
    #this_tile_id = 4


    ## Process ID
    process_id = "%04.0f_%s" % (this_tile_id , userinput["lr_large_image_name"].split("/")[-1].split(".fits")[0] ) # set the process ID to the name of the large image array + the tile number
    print("This Process is: " + str(process_id))


    ## Working directory
    dir_this_process = os.path.join(userinput["workdir"], "%s_%s" % (userinput["out_prefix"],process_id) )
    if not os.path.exists(dir_this_process):
        print("Looks like there is no data to plot")
        return(False)


    ### Things we need


    # 1. images (all in the final FITS product)
    # extensions: original, compl:[model, residual, ], segmap, segmap_sextractor , mask, mask_bad
    with fits.open(os.path.join(dir_this_process,"lr_tractor_results.fits") ) as hdul:
        lr_img = hdul[0].data
        lr_img_h = hdul[0].header
        lr_res = hdul[2].data
        lr_segmap = hdul[3].data
        lr_mask_bad = hdul[6].data

    lr_pixscale = np.abs(lr_img_h["CD1_1"]*3600) # pixel scale in arcsec/px
    lr_img_wcs = wcs.WCS(lr_img_h)

    # Check some stuff
    if "CD1_2" not in lr_img_h.keys():
        lr_img_h["CD1_2"] = 0
    if "CD2_1" not in lr_img_h.keys():
        lr_img_h["CD2_1"] = 0

    # Get Pixel noise
    tmp = clip(lr_img[lr_segmap == 0], n=3, niter=10)
    lr_img_pixnoise = tmp["stdev"]
    lr_img_medbkg = tmp["med"]
    print("Pixelnoise of low-res image: %g" % lr_img_pixnoise )
    print("Median background of low-res image: %g" % lr_img_medbkg )



    # 2. Table
    TABLE = Table.read( os.path.join(dir_this_process,"lr_compl_table_final.fits") )


    # 3. astrometric offset
    ## Get Astrometry offset in RA and DEC in pixels
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



    ### PLOT  #####


    ## 1) Plot spatial offset and histogram -----------------

    fig = plt.figure(figsize=(15,7))
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=0.025)


    ## spatial offset in tractor rest-frame
    ax1 = fig.add_subplot(1,2,1)

    ax1.grid(linestyle="--",zorder=1)

    Dx_before = TABLE["Dx_tractorfit.lr"] - delta_X_median_px
    Dy_before = TABLE["Dy_tractorfit.lr"] - delta_Y_median_px
    Dx_after = TABLE["Dx_tractorfit.lr"]
    Dy_after = TABLE["Dy_tractorfit.lr"]
    ax1.plot(Dx_before,Dy_before,"s",markersize=1.5,color=def_cols[1],alpha=0.7,label="Before astrometry offset applied")
    ax1.plot(Dx_after , Dy_after ,"o",markersize=1.5,color=def_cols[0],alpha=0.7,label="After astrometry offset applied")

    #ax1.plot(Dx - delta_X_median_px,Dy - delta_Y_median_px,"s",markersize=1.5,color=def_cols[1],alpha=0.7,label="Before astrometry offset applied")

    perc_x = np.percentile(Dx_after,q=(50-68/2,50,50+68/2))
    perc_y = np.percentile(Dy_after,q=(50-68/2,50,50+68/2))
    ax1.errorbar(perc_x[1],
                 perc_y[1],
                 xerr=[[ perc_x[1]-perc_x[0] ] , [ perc_x[2]-perc_x[1] ]],
                 yerr=[[ perc_y[1]-perc_y[0] ] , [ perc_y[2]-perc_y[1] ]],
                 fmt="o",markersize=10,capsize=4,color="red",zorder=1000)

    ax1.legend(loc="best",fontsize=10)

    ax1.set_title("Offset in pixels between ACS and best-fit HSC",fontsize=15)
    ax1.set_xlabel("X [pixels]",fontsize=15,labelpad=10)
    ax1.set_ylabel("Y [pixels]",fontsize=15,labelpad=10)
    ax1.xaxis.set_tick_params(labelsize=13)
    ax1.yaxis.set_tick_params(labelsize=13)
    ax1.tick_params(axis="both",which="minor",length=2)
    ax1.tick_params(axis="both",which="major",length=4)
    ax1.minorticks_on()

    ## histogram
    ax2 = fig.add_subplot(1,2,2)

    stats_bkg = clip(dat=lr_img[lr_segmap == 0],n=3,niter=15)
    stats_res = clip(dat=lr_res[lr_segmap == 0],n=3,niter=10)

    hist_bkg = lr_img.ravel()
    hist_res = lr_res.ravel()
    hist_noise = np.random.normal(loc=stats_res["med"],scale=stats_bkg["stdev"],size=len(hist_res))

    h_out_noise = ax2.hist(hist_noise,cumulative=False,density=False,histtype="stepfilled",
                 bins=np.linspace( stats_bkg["med"]-15*stats_bkg["stdev"] , stats_bkg["med"]+15*stats_bkg["stdev"] , num=250 ),
                 color="lightgray",
                 label="Simulated noise")

    h_out_bkg = ax2.hist(hist_bkg,cumulative=False,density=False,histtype="step",
                 bins=np.linspace( stats_bkg["med"]-15*stats_bkg["stdev"] , stats_bkg["med"]+15*stats_bkg["stdev"] , num=250 ),
                 color=def_cols[0],
                linewidth=2,
                label="Original image")

    h_out_res = ax2.hist(hist_res,cumulative=False,density=False,histtype="step",
                 bins=np.linspace( stats_bkg["med"]-15*stats_bkg["stdev"] , stats_bkg["med"]+15*stats_res["stdev"] , num=250 ),
                 color=def_cols[1],
                linewidth=2,
                label="Residual image")

    v = lr_res.ravel()
    o = lr_img.ravel()
    normalized_sum_of_residuals = np.nansum(v**2)/len(o)


    txt2 = r"$\frac{1}{N}\,\sum\,{\rm res^{2}}$ = %5.5f" % ( normalized_sum_of_residuals )
    ax2.text(0.95,0.95,txt2,transform=ax2.transAxes,va="center",ha="right",fontsize=12,color="black")

    ax2.legend(loc="upper left",fontsize=13)
    ax2.set_xlim( (stats_bkg["med"]-6*stats_bkg["stdev"],stats_bkg["med"]+6*stats_bkg["stdev"]))

    ax2.set_title("Comparison of residual to noise",fontsize=15)
    ax2.set_xlabel("Flux (arbitrary units)",fontsize=15,labelpad=10)
    ax2.set_ylabel("Number",fontsize=15,labelpad=10)
    ax2.xaxis.set_tick_params(labelsize=13)
    ax2.yaxis.set_tick_params(labelsize=13)
    ax2.tick_params(axis="both",which="minor",length=2)
    ax2.tick_params(axis="both",which="major",length=4)
    ax2.minorticks_on()

    outfile = os.path.join(dir_this_process,"lr_fig1.pdf")
    plt.savefig(outfile,bbox_inches='tight')
    plt.close()



    ## 2) Plot flux difference between lr and hr image

    #TAB_use = TABLE.copy()[TABLE["brightness.Flux.lr"] > 0]
    TAB_use = TABLE.copy()[(TABLE["flag_bad.lr"] == 0) & (TABLE["flag_detected.lr"] == 1) & (TABLE["flag_edge.lr"] == 0) & (TABLE["brightness.Flux.lr"] > 0) ]
    TAB_not_use = TABLE.copy()[(TABLE["flag_bad.lr"] != 0) | (TABLE["flag_detected.lr"] != 1) | (TABLE["flag_edge.lr"] != 0) | (TABLE["brightness.Flux.lr"] <= 0) ]


    #TAB_use = TABLE.copy()[(TABLE["brightness.Flux.lr"] > 0) & (TABLE["FLUX_AUTO.hr"] > 0)]

    TAB_notuse = TABLE.copy()[ (TABLE["flag_bad.lr"] == 1) | ( TABLE["flag_detected.lr"] == 0) | ( TABLE["flag_edge.lr"] == 1) ]

    ## let's make some selection

    # select compact objects (on ACS)
    sel_compact = np.where( (TAB_use["is_pointsource"] == 2) )[0]

    # select blended objects
    sel_blended = []
    for ii in range(len(TAB_use)):
        dist = np.sqrt( ( TAB_use["RA_tractor.hr"]-TAB_use["RA_tractor.hr"][ii] )**2 + ( TAB_use["DEC_tractor.hr"]-TAB_use["DEC_tractor.hr"][ii] )**2 )*3600
        if len( np.where(dist < 0.7)[0] ) > 1:
            sel_blended.append(ii)

    # select isolated objects
    sel_isolated = []
    for ii in range(len(TAB_use)):
        dist = np.sqrt( ( TAB_use["RA_tractor.hr"]-TAB_use["RA_tractor.hr"][ii] )**2 + ( TAB_use["DEC_tractor.hr"]-TAB_use["DEC_tractor.hr"][ii] )**2 )*3600
        if len( np.where(dist < 1.5)[0] ) == 1:
            sel_isolated.append(ii)



    xx = np.asarray(TAB_use["brightness.Flux.mujy.hr"])
    #xx = ConvertToMicroJansky(flux=np.asarray(TAB_use["FLUX_AUTO.hr"]), zp=userinput["hr_zp"] )
    yy = np.asarray(TAB_use["brightness.Flux.mujy.lr"])

    xx_all = np.asarray(TABLE["brightness.Flux.mujy.hr"])
    yy_all = np.asarray(TABLE["brightness.Flux.mujy.lr"])

    xx_notuse = np.asarray(TAB_not_use["brightness.Flux.mujy.hr"])
    yy_notuse = np.asarray(TAB_not_use["brightness.Flux.mujy.lr"])

    sel_bad_fit = np.where( np.abs((yy-xx)/xx) > 0.5 )[0]
    TAB_bad = TAB_use.copy()[sel_bad_fit]


    ## Plot
    fig = plt.figure(figsize=(7,5))
    #fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.0)


    ax1 = fig.add_subplot(1,1,1)

    ax1.semilogx([1e-3,1e4],[0,0],"--",dashes=(10,2),color="black",zorder=1,linewidth=0.5)
    ax1.semilogx([1e-3,1e4],[0.2,0.2],"--",dashes=(5,5),color="black",zorder=1,linewidth=0.5)
    ax1.semilogx([1e-3,1e4],[-0.2,-0.2],"--",dashes=(5,5),color="black",zorder=1,linewidth=0.5)
    ax1.semilogx([1e-3,1e4],[0.5,0.5],"--",dashes=(5,5),color="black",zorder=1,linewidth=0.5)
    ax1.semilogx([1e-3,1e4],[-0.5,-0.5],"--",dashes=(5,5),color="black",zorder=1,linewidth=0.5)

    #ax1.semilogx(xx,np.log10(xx/yy),"o",markersize=2)
    ax1.semilogx(xx,(yy-xx)/xx,"o",markerfacecolor="white",markeredgecolor=def_cols[0],markersize=3,label="all objects",zorder=1)
    ax1.semilogx(xx[sel_compact],((yy-xx)/xx)[sel_compact],"s",markersize=3,label="compact objects",zorder=2)
    ax1.semilogx(xx[sel_blended],((yy-xx)/xx)[sel_blended],"^",markersize=3,label="blended objects",zorder=3)
    ax1.semilogx(xx[sel_isolated],((yy-xx)/xx)[sel_isolated],"o",markersize=1,label="isolated objects",zorder=4)

    ax1.semilogx(xx_notuse,(yy_notuse - xx_notuse) / xx_notuse,"x",markerfacecolor="white",markeredgecolor="black",markersize=5,label="bad fits",zorder=1)

    if plot_ids:
    	for ii in range(len(TABLE)):
    	    ax1.text(xx_all[ii],((yy_all-xx_all)/xx_all)[ii] , TABLE["NUMBER.hr"][ii],va="bottom",ha="center",fontsize=8,color="black" )
        

    #sel_show = np.where( TAB_use["NUMBER.hr"] == [91] )[0]
    #ax1.semilogx(xx[sel_show],((yy-xx)/xx)[sel_show],"o",color="red",markersize=2)


    sel_tmp = np.where( (xx > 1e-1) & (xx < 1e3) )
    txt11 = r"$\sigma = %5.3f$" % ( np.std( ((yy-xx)/xx)[sel_tmp] ) )
    txt12 = r"$\left<log(F^{hr}\,/\,F^{lr})\right> = %5.3f$" % ( np.nanmedian( np.log10(xx/yy)[sel_tmp] ) )
    txt13 = r"$F_{\rm +50%s}=%5.2f%s$" % ("\%" ,  len(np.where( ((yy-xx)/xx) > 0.5)[0])/len(xx)*100 , "\%"  )
    txt14 = r"$F_{\rm -50%s}=%5.2f%s$" % ("\%" ,  len(np.where( ((yy-xx)/xx) < -0.5)[0])/len(xx)*100 , "\%"  )
    txt15 = r"$F_{\rm in\,20%s}=%5.2f%s$" % ("\%" ,  len(np.where( np.abs((yy-xx)/xx) < 0.2)[0])/len(xx)*100 , "\%"  )
    txt16 = r"$F_{\rm in\,50%s}=%5.2f%s$" % ("\%" ,  len(np.where( np.abs((yy-xx)/xx) < 0.5)[0])/len(xx)*100 , "\%"  )
    ax1.text(0.95,0.95,txt11,transform=ax1.transAxes,ha="right",fontsize=11,color="black")
    ax1.text(0.95,0.9,txt12,transform=ax1.transAxes,ha="right",fontsize=11,color="black")
    ax1.text(0.95,0.85,txt13,transform=ax1.transAxes,ha="right",fontsize=11,color="black")
    ax1.text(0.95,0.8,txt14,transform=ax1.transAxes,ha="right",fontsize=11,color="black")
    ax1.text(0.95,0.75,txt15,transform=ax1.transAxes,ha="right",fontsize=11,color="black")
    ax1.text(0.95,0.7,txt16,transform=ax1.transAxes,ha="right",fontsize=11,color="black")
    ax1.legend(loc="lower right",fontsize=10)
    ax1.set_xlim(0.5e-2,1e3)
    ax1.set_ylim(-3,3)
    ax1.xaxis.set_tick_params(labelsize=13)
    ax1.yaxis.set_tick_params(labelsize=13)
    ax1.tick_params(axis="both",which="minor",length=2)
    ax1.tick_params(axis="both",which="major",length=4)
    ax1.minorticks_on()
    ax1.set_xlabel(r"$F_{\rm tractor}^{hr}$ [$\mu$Jy]",fontsize=15)
    #ax1.set_ylabel(r"log $F_{\rm tractor}^{lr}\,/\,F_{\rm tractor}^{hr} $",fontsize=15)
    ax1.set_ylabel(r"$(F_{\rm tractor}^{lr}-F_{\rm tractor}^{hr})\,/\,F_{\rm tractor}^{hr} $",fontsize=15)

    outfile = os.path.join(dir_this_process,"lr_fig3.pdf" )
    plt.savefig(outfile,bbox_inches='tight')
    plt.close()




    ## 3) Plot Original and residual and mask.

    ima = dict(cmap=plt.get_cmap("Greys"),vmin=-2*lr_img_pixnoise,vmax=10*lr_img_pixnoise,origin="lower",interpolation="nearest")
    ima_sqrt = dict(cmap=plt.get_cmap("Greys"),vmin=-2*np.sqrt(lr_img_pixnoise),vmax=10*np.sqrt(lr_img_pixnoise),origin="lower",interpolation="nearest")


    fig = plt.figure(figsize=(20,15))
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.1, hspace=0.1)

    ## Original image
    ax1 = fig.add_subplot(2,2,1,projection = lr_img_wcs)
    ax1.imshow(lr_img,**ima)
    for ii in range(len(TABLE)):
        #ax1.plot(TABLE["X_IMAGE_hr.lr"][ii]-1,TABLE["Y_IMAGE_hr.lr"][ii]-1,"s",markersize=5,markeredgewidth=1,fillstyle="none",alpha=1,color=def_cols[ii % len(def_cols)])
        #ax1.plot(TABLE["X_IMAGE_tractor.lr"][ii]-1,TABLE["Y_IMAGE_tractor.lr"][ii]-1,"o",markersize=5,markeredgewidth=1.0,fillstyle="none",alpha=1,color=def_cols[ii % len(def_cols)])
        colid = 0
        if TABLE["is_pointsource"][ii] == 2:
        	colid = 2
        ax1.plot(TABLE["X_IMAGE_hr.lr"][ii]-1,TABLE["Y_IMAGE_hr.lr"][ii]-1,"s",markersize=5,markeredgewidth=1,fillstyle="none",alpha=1,color=def_cols[colid])
        ax1.plot(TABLE["X_IMAGE_tractor.lr"][ii]-1,TABLE["Y_IMAGE_tractor.lr"][ii]-1,"o",markersize=5,markeredgewidth=1.0,fillstyle="none",alpha=1,color=def_cols[colid])


    # mark all the bad fits
    ax1.plot(TAB_bad["X_IMAGE_tractor.lr"]-1,TAB_bad["Y_IMAGE_tractor.lr"]-1,"o",markersize=10,markeredgewidth=1,fillstyle="none",alpha=1,color="red")

    # mark all objects that are not considered
    ax1.plot(TAB_notuse["X_IMAGE_tractor.lr"]-1,TAB_notuse["Y_IMAGE_tractor.lr"]-1,"x",markersize=7,markeredgewidth=1,fillstyle="none",alpha=1,color="red")

    ax1.set_title("Original low-resolution image",fontsize=15)
    #ax1.set_xlabel("Right-ascension (J2000)",fontsize=12,labelpad=2)
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


    ## Residual image
    ax2 = fig.add_subplot(2,2,2,projection = lr_img_wcs)
    ax2.imshow(lr_res,**ima)
    for ii in range(len(TABLE)):
        #ax2.plot(TABLE["X_IMAGE_hr.lr"][ii]-1,TABLE["Y_IMAGE_hr.lr"][ii]-1,"s",markersize=5,markeredgewidth=1,fillstyle="none",alpha=1,color=def_cols[ii % len(def_cols)])
        #ax2.plot(TABLE["X_IMAGE_tractor.lr"][ii]-1,TABLE["Y_IMAGE_tractor.lr"][ii]-1,"o",markersize=5,markeredgewidth=1.0,fillstyle="none",alpha=1,color=def_cols[ii % len(def_cols)])
        colid = 0
        if TABLE["is_pointsource"][ii] == 2:
        	colid = 2
        ax2.plot(TABLE["X_IMAGE_hr.lr"][ii]-1,TABLE["Y_IMAGE_hr.lr"][ii]-1,"s",markersize=5,markeredgewidth=1,fillstyle="none",alpha=1,color=def_cols[colid])
        ax2.plot(TABLE["X_IMAGE_tractor.lr"][ii]-1,TABLE["Y_IMAGE_tractor.lr"][ii]-1,"o",markersize=5,markeredgewidth=1.0,fillstyle="none",alpha=1,color=def_cols[colid])

    # mark all the bad fits
    ax2.plot(TAB_bad["X_IMAGE_tractor.lr"]-1,TAB_bad["Y_IMAGE_tractor.lr"]-1,"o",markersize=10,markeredgewidth=1,fillstyle="none",alpha=1,color="red")

    # mark all objects that are not considered
    ax2.plot(TAB_notuse["X_IMAGE_tractor.lr"]-1,TAB_notuse["Y_IMAGE_tractor.lr"]-1,"x",markersize=7,markeredgewidth=1,fillstyle="none",alpha=1,color="red")


    ax2.set_title("Residual low-resolution image",fontsize=15)
    #ax2.set_xlabel("Right-ascension (J2000)",fontsize=12,labelpad=2)
    #ax2.set_ylabel("Declination (J2000)",fontsize=12,labelpad=2)
    ax2.xaxis.set_tick_params(labelsize=13)
    ax2.yaxis.set_tick_params(labelsize=13)
    ax2.tick_params(axis="both",which="minor",length=2)
    ax2.tick_params(axis="both",which="major",length=4)
    ax2.coords[0].set_major_formatter('hh:mm:ss.s')
    ax2.coords[1].set_major_formatter('dd:mm:ss.s')
    ax2.coords[0].set_ticklabel(size=11)
    ax2.coords[1].set_ticklabel(size=11)
    #ax2.coords[1].set_ticklabel_visible(False)
    ax2.minorticks_on()


    ## Segmentation map
    ax3 = fig.add_subplot(2,2,3,projection = lr_img_wcs)
    ax3.imshow(lr_segmap, origin="lower",cmap=plt.get_cmap("Greys") , vmin=0,vmax=1)
    for ii in range(len(TABLE)):
        #ax3.plot(TABLE["X_IMAGE_hr.lr"][ii]-1,TABLE["Y_IMAGE_hr.lr"][ii]-1,"s",markersize=5,markeredgewidth=1,fillstyle="none",alpha=1,color=def_cols[ii % len(def_cols)])
        #ax3.plot(TABLE["X_IMAGE_tractor.lr"][ii]-1,TABLE["Y_IMAGE_tractor.lr"][ii]-1,"o",markersize=5,markeredgewidth=1.0,fillstyle="none",alpha=1,color=def_cols[ii % len(def_cols)])
        colid = 0
        if TABLE["is_pointsource"][ii] == 2:
        	colid = 2
        ax3.plot(TABLE["X_IMAGE_hr.lr"][ii]-1,TABLE["Y_IMAGE_hr.lr"][ii]-1,"s",markersize=5,markeredgewidth=1,fillstyle="none",alpha=1,color=def_cols[colid])
        ax3.plot(TABLE["X_IMAGE_tractor.lr"][ii]-1,TABLE["Y_IMAGE_tractor.lr"][ii]-1,"o",markersize=5,markeredgewidth=1.0,fillstyle="none",alpha=1,color=def_cols[colid])

    # mark all the bad fits
    ax3.plot(TAB_bad["X_IMAGE_tractor.lr"]-1,TAB_bad["Y_IMAGE_tractor.lr"]-1,"o",markersize=10,markeredgewidth=1,fillstyle="none",alpha=1,color="red")

    # mark all objects that are not considered
    ax3.plot(TAB_notuse["X_IMAGE_tractor.lr"]-1,TAB_notuse["Y_IMAGE_tractor.lr"]-1,"x",markersize=7,markeredgewidth=1,fillstyle="none",alpha=1,color="red")


    ax3.set_title("Low-resolution segmentation map",fontsize=15)
    ax3.set_xlabel("Right-ascension (J2000)",fontsize=12,labelpad=2)
    ax3.set_ylabel("Declination (J2000)",fontsize=12,labelpad=2)
    ax3.xaxis.set_tick_params(labelsize=13)
    ax3.yaxis.set_tick_params(labelsize=13)
    ax3.tick_params(axis="both",which="minor",length=2)
    ax3.tick_params(axis="both",which="major",length=4)
    ax3.coords[0].set_major_formatter('hh:mm:ss.s')
    ax3.coords[1].set_major_formatter('dd:mm:ss.s')
    ax3.coords[0].set_ticklabel(size=11)
    ax3.coords[1].set_ticklabel(size=11)
    #ax3.coords[1].set_ticklabel_visible(False)
    ax3.minorticks_on()



    ## Mask map
    ax4 = fig.add_subplot(2,2,4,projection = lr_img_wcs)
    ax4.imshow(lr_mask_bad, origin="lower",cmap=plt.get_cmap("Greys") , vmin=0,vmax=1)
    for ii in range(len(TABLE)):
        #ax4.plot(TABLE["X_IMAGE_hr.lr"][ii]-1,TABLE["Y_IMAGE_hr.lr"][ii]-1,"s",markersize=5,markeredgewidth=1,fillstyle="none",alpha=1,color=def_cols[ii % len(def_cols)])
        #ax4.plot(TABLE["X_IMAGE_tractor.lr"][ii]-1,TABLE["Y_IMAGE_tractor.lr"][ii]-1,"o",markersize=5,markeredgewidth=1.5,fillstyle="none",alpha=1,color=def_cols[ii % len(def_cols)])
        colid = 0
        if TABLE["is_pointsource"][ii] == 2:
        	colid = 2
        ax4.plot(TABLE["X_IMAGE_hr.lr"][ii]-1,TABLE["Y_IMAGE_hr.lr"][ii]-1,"s",markersize=5,markeredgewidth=1,fillstyle="none",alpha=1,color=def_cols[colid])
        ax4.plot(TABLE["X_IMAGE_tractor.lr"][ii]-1,TABLE["Y_IMAGE_tractor.lr"][ii]-1,"o",markersize=5,markeredgewidth=1.0,fillstyle="none",alpha=1,color=def_cols[colid])

    # mark all the bad fits
    ax4.plot(TAB_bad["X_IMAGE_tractor.lr"]-1,TAB_bad["Y_IMAGE_tractor.lr"]-1,"o",markersize=10,markeredgewidth=1,fillstyle="none",alpha=1,color="red")

    # mark all objects that are not considered
    ax4.plot(TAB_notuse["X_IMAGE_tractor.lr"]-1,TAB_notuse["Y_IMAGE_tractor.lr"]-1,"x",markersize=7,markeredgewidth=1,fillstyle="none",alpha=1,color="red")


    ax4.set_title("Low-resolution mask map",fontsize=15)
    ax4.set_xlabel("Right-ascension (J2000)",fontsize=12,labelpad=2)
    #ax4.set_ylabel("Declination (J2000)",fontsize=12,labelpad=2)
    ax4.xaxis.set_tick_params(labelsize=13)
    ax4.yaxis.set_tick_params(labelsize=13)
    ax4.tick_params(axis="both",which="minor",length=2)
    ax4.tick_params(axis="both",which="major",length=4)
    ax4.coords[0].set_major_formatter('hh:mm:ss.s')
    ax4.coords[1].set_major_formatter('dd:mm:ss.s')
    ax4.coords[0].set_ticklabel(size=11)
    ax4.coords[1].set_ticklabel(size=11)
    #ax4.coords[1].set_ticklabel_visible(False)
    ax4.minorticks_on()




    outfile = os.path.join(dir_this_process,"lr_fig2.pdf")
    plt.savefig(outfile,bbox_inches='tight')
    plt.close()


    print("DONE!")

    return(True)
