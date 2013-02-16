# -*- coding: utf-8 -*-
import AllocationModule as AM
import sys
import numpy as np
sys.path.append('../../Wrapper/')
import archive

def main():

    #load all the inputs
    data_bank = '../../../data/data_bank.h5'
    with archive.archive(data_bank, 'r') as ar:
        # numbers coming from param.ini
        fiber_pitch 	= np.float_(ar['/fiber/parameters/fiber_pitch'])
        fiber_size 	= np.float_(ar['/fiber/parameters/fiber_size'])
        Ndiameter 	= np.int_(ar['/fiber/parameters/num_fibers_on_diameter'])
        n_pass_per_tile = np.int_(ar['/allocation/parameters/n_pass_per_tile'])
        patrol_radius   = np.float_(ar['/allocation/parameters/patrol_radius'])
        allocation_method = ar['/allocation/parameters/allocation_method']

        # numbers coming from the previous step in the pipeline
        gals_x_in	= ar['/gal/ra_true']		
        gals_y_in	= ar['/gal/dec_true']		
        gals_id_in	= ar['/gal/galaxy_index']
        number_of_tiles = ar['/tiling/number_of_tiles'] 
        tiles_centers_ra = ar['/tiling/tile_centers_ra'] 
        tiles_centers_dec = ar['/tiling/tile_centers_dec'] 
        tile_observation_time =ar['/tiling/tile_observation_time'] 
        
    #initialize the galaxy structures
    all_gals = AM.MockGalaxyCatalog(x_in=gals_x_in, y_in=gals_y_in, id_in=gals_id_in)

    #loop over the tiles
    tile_ID = 0
    for ra_center, dec_center in zip(tiles_centers_ra, tiles_centers_dec):
        for pass_ID in range(n_pass_per_tile):
        # create the unperturbed set of fibers for this pass in this tile
            fibers = AM.FiberSet(Ndiameter=Ndiameter, fiber_pitch=fiber_pitch,\
                             fiber_size=fiber_size, center_x=ra_center, center_y=dec_center)
            AM.make_fiber_allocation(fibers, all_gals, tile_visit_ID=tile_ID, rank_criterion=allocation_method,\
                                     patrol_radius=fibers.fiber_pitch, exclusion_radius=fibers.fiber_size * 0.5)

        #update the tile_ID
        tile_ID = tile_ID + 1

# ===================================================== 
if __name__ == '__main__':
    main()


