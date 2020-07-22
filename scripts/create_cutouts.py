## Function to get grid point (center and extend + margin) of a tile from a job file
def get_tile_center_for_job(jobfile_name):
    '''
    For a given job file (including path), this function returns the center and width (+margin)
    of the tile run by the job file
    '''

    ## Read JSON file
    with open(jobfile_name,"r") as json_file:
            job_input = json.load(json_file)

    ## Get necessary information
    cutoutsize_arcsec = job_input["cutoutsize_arcsec"] # in arcseconds
    margin_arcsec = job_input["cutoutoverlap_arcsec"] # in arcseconds
    tile_nbr = job_input["tile_id"]
    img_large_lr_path = job_input["lr_large_image_name"]
    
    ## Read LR image (only header)
    with fits.open(img_large_lr_path) as hdul:
        img_large_lr_header = hdul[1].header
    lr_pixscale = np.abs(img_large_lr_header["CD1_1"]*3600) # pixel scale in arcsec/px

    ## Create grid (centers)

    # get width
    width_RA = img_large_lr_header["NAXIS2"]*lr_pixscale/60 # in arcminutes
    width_DEC = img_large_lr_header["NAXIS1"]*lr_pixscale/60 # in arcminutes

    # get number of tiles that fit the low-res image (high-res image has to cover same area or larger)
    ntiles_RA = int(np.ceil(width_RA / (cutoutsize_arcsec[0]/60))) # width is in arcmin and cutout size in arcsec -> convert
    ntiles_RA_rest = np.abs(ntiles_RA*cutoutsize_arcsec[0] - width_RA*60) # in arcsec
    ntiles_DEC = int(np.ceil(width_DEC / (cutoutsize_arcsec[1]/60))) # width is in arcmin and cutout size in arcsec -> convert
    ntiles_DEC_rest = np.abs(ntiles_DEC*cutoutsize_arcsec[1] - width_DEC*60) # in arcsec
    
    # get starting point
    ref_point_lr = wcs.WCS(img_large_lr_header).all_pix2world(0,0,0)

    # get grids for lr (used in the end)
    grid_lr_ra = np.asarray([ ref_point_lr[0] + (cutoutsize_arcsec[0]/3600)/2 - (cutoutsize_arcsec[0]/3600)*(ii+1) for ii in range(ntiles_RA) ]) + (ntiles_RA_rest/3600)/2
    grid_lr_dec = np.asarray([ ref_point_lr[1] - (cutoutsize_arcsec[1]/3600)/2 + (cutoutsize_arcsec[1]/3600)*(ii+1) for ii in range(ntiles_DEC) ]) - (ntiles_DEC_rest/3600)/2
    grid_lr_all = np.meshgrid(grid_lr_ra,grid_lr_dec)
    grid_lr_all[0] = grid_lr_all[0].flatten()
    grid_lr_all[1] = grid_lr_all[1].flatten()

    grid_lr = [[],[]]
    grid_lr[0] = np.asarray([grid_lr_all[0][tile_nbr]])
    grid_lr[1] = np.asarray([grid_lr_all[1][tile_nbr]])
        
    out = {"ra_center":grid_lr[0][0],
          "dec_center":grid_lr[1][0],
           "ra_width_arcsec":cutoutsize_arcsec[0],
           "dec_width_arcsec":cutoutsize_arcsec[1],
            "ra_width_margin_arcsec":cutoutsize_arcsec[0]+margin_arcsec[0],
           "dec_width_margin_arcsec":cutoutsize_arcsec[1]+margin_arcsec[1],
          }

    #plt.plot(grid_lr_all[0] , grid_lr_all[1] , "o")
    #plt.plot(out["ra_center"] , out["dec_center"] , "o" , color="red" , markersize=10)

    return(out)




## Cut out an image slice.
# Center is [RA,DEC] in degrees
# size is [dRA,dDEC] in arcseconds
# outfile is the path and name where the file is saved to
def cutout_to_file(img,img_header,center,size,outfile,frame='icrs',overwrite=True,mode="partial", fill_value=np.nan):
    
    ## Cut out
    position = SkyCoord(center[0],center[1], unit="deg", frame=frame)
    size = u.Quantity((size[0], size[1]), u.arcsec)
    cutout = Cutout2D(data=img, position=position, size=size, wcs=wcs.WCS(img_header) , mode=mode, fill_value=fill_value)
    
    ## Put the cutout image in the FITS HDU
    hdu_new = fits.PrimaryHDU(cutout.data)
    hdu_new.header = img_header.copy()
    
    ## Add new header (meaning: update it)
    hdu_new.header.update(cutout.wcs.to_header())
    
    ## Add to hdu list and save
    hdul_new = fits.HDUList([hdu_new])    
    hdul_new.writeto(outfile,overwrite=overwrite)
    
    return(True)


## Cut out an image slice. Same as "cutout_to_file" but also makes cutout of mask image and puts them together with hdu0=image, hdu1=mask
# Center is [RA,DEC] in degrees
# size is [dRA,dDEC] in arcseconds
# outfile is the path and name where the file is saved to
def cutout_to_file_with_mask(img,img_header,mask,mask_header,center,size,outfile,frame='icrs',overwrite=True,mode="partial", fill_value=np.nan):
    
    ## Cut out
    position = SkyCoord(center[0],center[1], unit="deg", frame=frame)
    size = u.Quantity((size[0], size[1]), u.arcsec)
    cutout = Cutout2D(data=img, position=position, size=size, wcs=wcs.WCS(img_header) , mode=mode, fill_value=fill_value)
    cutout_mask = Cutout2D(data=mask, position=position, size=size, wcs=wcs.WCS(img_header) , mode=mode, fill_value=fill_value)
    
    ## Put the cutout image in the 0th FITS HDU
    hdu_new = fits.PrimaryHDU(cutout.data)
    hdu_new.header = img_header.copy()
    
    ## Put the cutout mask in the 1st FITS HDU
    hdu1_new = fits.ImageHDU(cutout_mask.data)
    hdu1_new.header = mask_header.copy()
    hdu1_new.header["EXTNAME"] = "MASK"
    
    
    ## Add new header (meaning: update it)
    hdu_new.header.update(cutout.wcs.to_header())
    hdu1_new.header.update(cutout_mask.wcs.to_header())
    
    ## Add to hdu list and save
    hdu_new.verify('silentfix')
    hdu1_new.verify('silentfix')
    hdul_new = fits.HDUList([hdu_new , hdu1_new])
    hdul_new.writeto(outfile,overwrite=overwrite)
    
    return(True)



## This function cuts a big tile in small tiles given a cutout size in arcseconds.
# number_of_cutouts = Number of cutouts. Set to -99 if all should be generated.
# MODIFIED: does not use HR weight image!
def create_cutouts(userinput,overwrite_cutouts,tile_nbr):
    
    ## Get necessary information
    cutoutsize_arcsec = userinput["cutoutsize_arcsec"] # in arcseconds
    margin_arcsec = userinput["cutoutoverlap_arcsec"] # in arcseconds
    img_large_hr_path = userinput["hr_large_image_name"]
    img_large_lr_path = userinput["lr_large_image_name"]
    #print(cutoutsize_arcsec)
    
    ## Create directory
    cutout_output_dir = os.path.join( userinput["cutout_output_dir"] , img_large_lr_path.split("/")[-1].split(".fits")[0] )
    sh.mkdir('-p', cutout_output_dir)
    
    ## Load large images
    
    # 1) high-resolution image
    with fits.open(img_large_hr_path) as hdul:
        img_large_hr = hdul[0].data
        img_large_hr_header = hdul[0].header
    hr_pixscale = np.abs(img_large_hr_header["CD1_1"]*3600) # pixel scale in arcsec/px
    
    
    # 2) low-resolution image (need to fix header because SIMPLE keyword is in primary hdu but is needed by cutout2D. Therefore add it to 1 hdu)
    with fits.open(img_large_lr_path) as hdul:
        img_large_lr = hdul[1].data
        img_large_lr_header = hdul[1].header
        img_large_lr_header0 = hdul[0].header
        mask_large_lr = hdul[2].data
        mask_large_lr_header = hdul[2].header
    keys_add = ["SIMPLE"]
    for key in keys_add:
        img_large_lr_header.append(  (key, img_large_lr_header0[key], img_large_lr_header0.comments[key] )  )
    hdu_tmp = fits.PrimaryHDU() # create placeholder hdu to run "silentfix"
    hdu_tmp.header = img_large_lr_header
    hdu_tmp.verify('silentfix')
    img_large_lr_header = hdu_tmp.header.copy()    
    lr_pixscale = np.abs(img_large_lr_header["CD1_1"]*3600) # pixel scale in arcsec/px
    
    
    ## Create grid (centers)
    
    # get width
    width_RA = img_large_lr.shape[0]*lr_pixscale/60 # in arcminutes
    width_DEC = img_large_lr.shape[1]*lr_pixscale/60 # in arcminutes
    
    # get number of tiles that fit the low-res image (high-res image has to cover same area or larger)
    ntiles_RA = int(np.ceil(width_RA / (cutoutsize_arcsec[0]/60))) # width is in arcmin and cutout size in arcsec -> convert
    ntiles_RA_rest = np.abs(ntiles_RA*cutoutsize_arcsec[0] - width_RA*60) # in arcsec
    ntiles_DEC = int(np.ceil(width_DEC / (cutoutsize_arcsec[1]/60))) # width is in arcmin and cutout size in arcsec -> convert
    ntiles_DEC_rest = np.abs(ntiles_DEC*cutoutsize_arcsec[1] - width_DEC*60) # in arcsec
    
    # get starting point
    #ref_point_hr = wcs.WCS(img_large_hr_header).all_pix2world(0,0,0)
    ref_point_lr = wcs.WCS(img_large_lr_header).all_pix2world(0,0,0)
    
    # get grids for hr (not needed in the end)
    #grid_hr_ra = np.asarray([ ref_point_hr[0] + (cutoutsize_arcsec[0]/3600)/2 - (cutoutsize_arcsec[0]/3600)*(ii+1) for ii in range(ntiles_RA) ]) + (ntiles_RA_rest/3600)/2
    #grid_hr_dec = np.asarray([ ref_point_hr[1] - (cutoutsize_arcsec[1]/3600)/2 + (cutoutsize_arcsec[1]/3600)*(ii+1) for ii in range(ntiles_DEC) ]) - (ntiles_DEC_rest/3600)/2
    #grid_hr = np.meshgrid(grid_hr_ra,grid_hr_dec)
    #grid_hr[0] = grid_hr[0].flatten()
    #grid_hr[1] = grid_hr[1].flatten()
    
    # get grids for lr (used in the end)
    grid_lr_ra = np.asarray([ ref_point_lr[0] + (cutoutsize_arcsec[0]/3600)/2 - (cutoutsize_arcsec[0]/3600)*(ii+1) for ii in range(ntiles_RA) ]) + (ntiles_RA_rest/3600)/2
    grid_lr_dec = np.asarray([ ref_point_lr[1] - (cutoutsize_arcsec[1]/3600)/2 + (cutoutsize_arcsec[1]/3600)*(ii+1) for ii in range(ntiles_DEC) ]) - (ntiles_DEC_rest/3600)/2
    grid_lr_all = np.meshgrid(grid_lr_ra,grid_lr_dec)
    grid_lr_all[0] = grid_lr_all[0].flatten()
    grid_lr_all[1] = grid_lr_all[1].flatten()
    

    # test if the tile id is valid
    if tile_nbr >= len(grid_lr_all[0]):
    	print("Tile ID is invalid")
    	out = {"ntiles_RA":ntiles_RA,
    	"ntiles_DEC":ntiles_DEC,
    	"cutout_output_dir":cutout_output_dir,
    	"tile_nbr":tile_nbr,
    	"success":False
    	}
    	return(out)


    # extract part of grid that should be used to create cutouts
    grid_lr = [[],[]]
    grid_lr[0] = np.asarray([grid_lr_all[0][tile_nbr]])
    grid_lr[1] = np.asarray([grid_lr_all[1][tile_nbr]])

    
    ## Create Cutouts (use lr grid!)

    # 1) for high-res (including weight image)
    print("Creating tile %g (out of in total %g cutouts) for high-res image (including weight images) . . . " % (tile_id,ntiles_RA*ntiles_DEC), end="\n")
    for ra_this,dec_this,tile_nbr_this in zip(grid_lr[0],grid_lr[1],[tile_nbr]):
        
        if (overwrite_cutouts == False) & (os.path.exists( os.path.join(cutout_output_dir , "%04.0f_%s.fits" % ( tile_nbr_this, userinput["hr_large_image_name"].split("/")[-1].split(".fits")[0]) ) )):
            print(" image file exists, skip. ",end="")
            pass # file exists, therefore do not create cutout
        else:
            cutout_to_file(img=img_large_hr.copy(),
                           img_header=img_large_hr_header.copy(),
                          center=[ra_this,dec_this],
                          size=[cutoutsize_arcsec[0] + margin_arcsec[0] , cutoutsize_arcsec[1] + margin_arcsec[1]],
                          outfile= os.path.join(cutout_output_dir , "%04.0f_%s.fits" % ( tile_nbr_this, userinput["hr_large_image_name"].split("/")[-1].split(".fits")[0]) ),
                           mode="trim",
                           fill_value = np.nan
                          )

    print("done!")
    
    # 2) for low-res (here now also cut the mask)
    print("Creating tile %g (out of in total %g cutouts) for low-res image . . . " % (tile_id,ntiles_RA*ntiles_DEC), end="\n")

    for ra_this,dec_this,tile_nbr_this in zip(grid_lr[0],grid_lr[1],[tile_nbr]):
        
        if (overwrite_cutouts == False) &  (os.path.exists( os.path.join(cutout_output_dir , "%04.0f_%s.fits" % ( tile_nbr_this, userinput["lr_large_image_name"].split("/")[-1].split(".fits")[0]) ) )):
            print(" image file exists, skip.",end="")
            pass # file exists, therefore do not create cutout
        else:
            cutout_to_file_with_mask(img=img_large_lr.copy(),
                                     img_header=img_large_lr_header.copy(),
                                     mask=mask_large_lr.copy(),
                                     mask_header=mask_large_lr_header.copy(),
                          center=[ra_this,dec_this],
                          size=[cutoutsize_arcsec[0] + margin_arcsec[0] , cutoutsize_arcsec[1] + margin_arcsec[1]],
                           outfile= os.path.join(cutout_output_dir , "%04.0f_%s.fits" % ( tile_nbr_this, userinput["lr_large_image_name"].split("/")[-1].split(".fits")[0]) ),
                           mode="trim",
                           fill_value = np.nan
                          )
        
    print("done!")
    
    
    ## propagate parameters
    out = {"ntiles_RA":ntiles_RA,
           "ntiles_DEC":ntiles_DEC,
           "cutout_output_dir":cutout_output_dir,
           "tile_nbr":tile_nbr,
           "success":True
          }
    
    return(out)
