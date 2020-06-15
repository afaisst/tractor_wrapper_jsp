## This function runs the code
def runTractor(userinput , tileids_array, doplots):
    
    
    for tileid in list(tileids_array):
        
        print("WORKING ON TILE ID %g" % tileid)
        
        get_HR_model(userinput=userinput,tileids=[tileid])
        
        if doplots:
        	plot_hr_images(userinput=userinput , tileid=tileid , plot_ids=False)
        
        create_LR_segmap_from_SExtractor(userinput=userinput,tileid=tileid , doplot=doplots)
        
        #create_mock_image_from_HR_image(userinput=userinput , tileid=tileid , astrometric_offset_sigma_mas=[50,50])
        
        fit_LR_image(userinput=userinput , tileid=tileid , seg_map_type="sextractor")
        
        if doplots:
        	plot_lr_images(userinput=userinput , tileid=tileid)
        
        run_sextractor_residual(userinput=userinput,tileid=tileid , doplot=doplots)
        
        print("DONE (TILE ID %g)" % tileid)
    
    return(True)




## This is the Tractor helper function that starts everyting
def runTractor_helper(userinput , usemultiproc , tile_id_list , doplots):

	print("+++++ RUNNING TRACTOR (single core) . . . +++++")
	start_time = time.time() # get time
	results = runTractor(userinput,[tile_id_list],doplots)
	elapsed_time = time.time() - start_time
	print("\n\n Finished. TOTAL ELAPSED TIME: " + str(round(elapsed_time/60,2)) + " minutes")
	print("+++++ DONE! +++++")

	return(True)


    