### CREATE HR IMAGES FOR CHECKING ####

# function making the plots
def mkplot_hr(data,residual,datatable,pixnoise,wcs,outfile,plot_ids=False):
    
    ima = dict(cmap=plt.get_cmap("Greys"),vmin=-2*pixnoise,vmax=10*pixnoise,origin="lower",interpolation="nearest")
    ima_sqrt = dict(cmap=plt.get_cmap("Greys"),vmin=-2*np.sqrt(pixnoise),vmax=10*np.sqrt(pixnoise),origin="lower",interpolation="nearest")

    v = residual.ravel()
    o = data.ravel()
    normalized_sum_of_residuals = np.nansum(v**2)/len(o)

    
    
    #sel_good = np.where( (datatable["RA_tractor.hr"] != -99) & (datatable["DEC_tractor.hr"] != -99) )
    datatable = datatable[np.where( (datatable["RA_tractor.hr"] != -99) & (datatable["DEC_tractor.hr"] != -99) )]
    
    sel_pointsource = np.where(datatable["is_pointsource"] == 2)

    ## FIGURE 1
    fig = plt.figure(figsize=(22,20))
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=0.025)

    # original
    ax1 = plt.subplot(1,2,1,projection=wcs)
    #ax1.imshow(np.sqrt(data), **ima_sqrt)
    ax1.imshow(data, **ima)
    ax1.plot(datatable["ALPHA_J2000.hr"],datatable["DELTA_J2000.hr"],"o",markersize=10,fillstyle="none",color=def_cols[0],transform=ax1.get_transform('world'),label="SExtractor")
    ax1.plot(datatable["RA_tractor.hr"],datatable["DEC_tractor.hr"],"o",markersize=14,fillstyle="none",color=def_cols[1],transform=ax1.get_transform('world'),label="Tractor")
    ax1.plot(datatable["RA_tractor.hr"][sel_pointsource],datatable["DEC_tractor.hr"][sel_pointsource],"s",markersize=14,fillstyle="none",color=def_cols[2],transform=ax1.get_transform('world'),label="point source")
    if plot_ids:
        for iii in range(len(datatable)):
            text = ax1.text(datatable["ALPHA_J2000.hr"][iii],datatable["DELTA_J2000.hr"][iii]+0.7/3600,datatable["NUMBER.hr"][iii],fontsize=7,color="white",va="bottom",ha="left",transform=ax1.get_transform('world'))
            text.set_path_effects([path_effects.Stroke(linewidth=1.0, foreground='black'),
                               path_effects.Normal()])

    ax1.set_title("Original",fontsize=15)
    ax1.set_xlabel("Right-ascension",fontsize=15)
    ax1.set_ylabel("Declination",fontsize=15)
    ax1.tick_params(axis="both",which="major",labelsize=15)
    ax1.coords[0].set_major_formatter('hh:mm:ss.s')
    ax1.coords[1].set_major_formatter('dd:mm:ss.s')
    ax1.coords[0].set_ticklabel(size=11)
    ax1.coords[1].set_ticklabel(size=11)
    ax1.legend(loc="best",fontsize=10,bbox_to_anchor=(1.015, 1.0))

    # residual
    ax2 = plt.subplot(1,2,2,projection=wcs)
    #ax2.imshow(np.sqrt(residual), **ima_sqrt)
    ax2.imshow(residual, **ima)
    ax2.plot(datatable["ALPHA_J2000.hr"],datatable["DELTA_J2000.hr"],"o",markersize=10,fillstyle="none",color=def_cols[0],transform=ax2.get_transform('world'),label="SExtractor")
    ax2.plot(datatable["RA_tractor.hr"],datatable["DEC_tractor.hr"],"o",markersize=14,fillstyle="none",color=def_cols[1],transform=ax2.get_transform('world'),label="Tractor")
    ax2.plot(datatable["RA_tractor.hr"][sel_pointsource],datatable["DEC_tractor.hr"][sel_pointsource],"s",markersize=14,fillstyle="none",color=def_cols[2],transform=ax2.get_transform('world'),label="point source")
    ax2.set_title(r"Residual ($\frac{1}{N}\,\sum\,{\rm res^{2}}$ = %5.5f)" % (normalized_sum_of_residuals) , fontsize=15)
    ax2.set_xlabel("Right-ascension",fontsize=15)
    ax2.set_ylabel("Declination",fontsize=15)
    ax2.coords[0].set_major_formatter('hh:mm:ss.s')
    ax2.coords[1].set_major_formatter('dd:mm:ss.s')
    ax2.coords[0].set_ticklabel(size=11)
    ax2.coords[1].set_ticklabel(size=11)
    #ax2.legend(loc="best",fontsize=10)

    plt.savefig(outfile,bbox_inches='tight')
    plt.close()



def plot_hr_images(userinput , tileid , plot_ids):
    
    print("PLOTTING HR IMAGES")

    ## 1. Create process ID ===========================
    this_tile_id = tileid
    #this_tile_id = 4
    #plot_ids = False


    ## Process ID
    process_id = "%04.0f_%s" % (this_tile_id , userinput["lr_large_image_name"].split("/")[-1].split(".fits")[0] ) # set the process ID to the name of the large image array + the tile number
    print("This Process is: " + str(process_id))


    ## Working directory
    dir_this_process = os.path.join(userinput["workdir"], "%s_%s" % (userinput["out_prefix"],process_id) )
    if not os.path.exists(dir_this_process):
        print("Looks like there is no data to plot")
        return(False) # UNCOMMENT



    ### things we need


    # 1. images (all in the final FITS product)
    # extensions: original, compl:[model, residual], segmentation
    with fits.open(os.path.join(dir_this_process,"hr_tractor_results.fits") ) as hdul:
        hr_img = hdul[0].data
        hr_img_h = hdul[0].header
        hr_res = hdul[2].data
        hr_segmap = hdul[3].data

    hr_pixscale = np.abs(hr_img_h["CD1_1"]*3600) # pixel scale in arcsec/px
    hr_img_wcs = wcs.WCS(hr_img_h)

    # Check some stuff
    if "CD1_2" not in hr_img_h.keys():
        hr_img_h["CD1_2"] = 0
    if "CD2_1" not in hr_img_h.keys():
        hr_img_h["CD2_1"] = 0

    # Get Pixel noise
    tmp = clip(hr_img[hr_segmap == 0], n=3, niter=10)
    hr_img_pixnoise = tmp["stdev"]
    hr_img_medbkg = tmp["med"]
    print("Pixelnoise of low-res image: %g" % hr_img_pixnoise )
    print("Median background of low-res image: %g" % hr_img_medbkg )



    # 2. Table
    TABLE = Table.read( os.path.join(dir_this_process,"hr_compl_table_final.fits") )


    ## PLOT ###
    mkplot_hr(data=hr_img,
           residual=hr_res,
           datatable=TABLE,
           pixnoise=hr_img_pixnoise,
           wcs=hr_img_wcs,
           outfile="%s/hr_fig1.pdf" % (dir_this_process),
            plot_ids=plot_ids
          )


    print("DONE!")

    return(True)


