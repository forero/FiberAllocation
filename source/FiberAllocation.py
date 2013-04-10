# -*- coding: utf-8 -*-
import sys
sys.path.append('../../Wrapper/')
sys.path.append('../../Plotters/')
sys.path.append('../../Utilities/')

import AllocationModule as AM
import numpy as np
import math

import diagnostic_plot as diagplot

import Plotters
import archive
import Utilities




def main():

	print "Started fiber allocation. Importing from data bank..."
	#load all the inputs
	data_bank = '../../../data/data_bank.h5'
	with archive.archive(data_bank, 'r') as ar:
		# numbers coming from param.ini
		fiber_pitch 	  = ar['/Fiber_Allocation/fiber_pitch']
		fiber_size 		  = ar['/Fiber_Allocation/fiber_size']
		Ndiameter 		  = ar['/Fiber_Allocation/num_fibers_on_diameter']
		n_pass_per_tile   = ar['/Fiber_Allocation/n_pass_per_tile']
		patrol_radius	 = ar['/Fiber_Allocation/patrol_radius']
		allocation_method = ar['/Fiber_Allocation/allocation_method']

		# numbers coming from the previous step in the pipeline
		gals_x_in				= ar['/gal/ra_true']
		gals_y_in				= ar['/gal/dec_true']
		gals_id_in				= ar['/gal/galaxy_index']
		number_of_tiles				= ar['/tiling/number_of_tiles']
		tiles_centers_ra		= ar['/tiling/tile_centers_ra']
		tiles_centers_dec			= ar['/tiling/tile_centers_dec']
		galaxy_tile_id_list			= ar['/gal/tile_ID']
		tile_id_list			= ar['/tiling/tile_ID']
		survey_selection_flag = ar['/gal/survey_selection_flag']

		survey_selection_flag = np.array(survey_selection_flag)
		gals_x_in = gals_x_in[survey_selection_flag]
		gals_y_in = gals_y_in[survey_selection_flag]
		gals_id_in = gals_id_in[survey_selection_flag]
		galaxy_tile_id_list = galaxy_tile_id_list[survey_selection_flag]



	print "...done."


	selected = np.where(galaxy_tile_id_list>0)
	n_in_tiles = np.size(selected)
	print 'There are %d galaxies in tiles'%(n_in_tiles)
	print 'Min-Max RA Tiles: %f %f'%(np.amin(tiles_centers_ra), np.amax(tiles_centers_ra))
	print 'Min-Max DEC Tiles: %f %f'%(np.amin(tiles_centers_dec), np.amax(tiles_centers_dec))
	print 'Min-Max RA Galaxies: %f %f'%(np.amin(gals_x_in), np.amax(gals_x_in))
	print 'Min-Max DEC Galaxies: %f %f'%(np.amin(gals_y_in), np.amax(gals_y_in))

#	plot_name = "all_galaxies"
#	diagplot.plot_positions(gals_x_in[selected], gals_y_in[selected], plot_name)

	# redefine the positions of the galaxies and tile centers to a 'planar' frame
	gals_x_in = gals_x_in * np.cos(gals_y_in * math.pi / 180.0)
	tiles_centers_ra  = tiles_centers_ra * np.cos(tiles_centers_dec * math.pi / 180.0)

	#initialize full galaxy structure
	all_gals = AM.MockGalaxyCatalog(x_in=gals_x_in, y_in=gals_y_in, id_in=gals_id_in)

	fraction_allocated = np.empty((0))

	#loop over the tiles
	for ra_center, dec_center, tile_id in zip(tiles_centers_ra, tiles_centers_dec, tile_id_list):
		# create the unperturbed set of fibers for this pass in this tile
		fibers = AM.FiberSet(Ndiameter=Ndiameter, fiber_pitch=fiber_pitch,\
						 fiber_size=fiber_size, center_x=ra_center, center_y=dec_center)

		# set to -1 the fiber_ID of the galaxies that should fall inside a tile
		index_inside = np.where((galaxy_tile_id_list==tile_id) & (all_gals.tile_ID < 0))
		all_gals.tile_ID[index_inside] = -1

		if np.size(index_inside)>0:
			print "Galaxies to allocate", np.size(index_inside)

			# Finally. Make the allocation!
			AM.make_fiber_allocation(fibers, all_gals, tile_ID=tile_id, visit_ID=1, rank_criterion=allocation_method,\
						 patrol_radius=fibers.fiber_pitch, exclusion_radius=fibers.fiber_size * 0.5)

			n_alloc = np.where(all_gals.fiber_ID>0)

		#update the number of allocated galaxies for a control plot
		gal_alloc = np.where(all_gals.fiber_ID>0)
		n_total = np.size(all_gals.fiber_ID)
		n_alloc = np.size(gal_alloc)
		fraction_allocated = np.append(fraction_allocated, (1.0 * n_alloc / (1.0*n_total)))



	# final stats
	gal_alloc = np.where(all_gals.fiber_ID>0)
	n_total = np.size(all_gals.fiber_ID)
	n_alloc = np.size(gal_alloc)

	print "Allocation efficiency", (1.0*n_alloc)/(1.0*n_total)
	print "Initial number of galaxies (after target selection)", n_total
	print "Initial number of galaxies (in tiles)", n_in_tiles
	print "Galaxies that were allocated:", n_alloc



	inner_fiber_selection_flag = np.array((all_gals.fiber_ID > 0))
	print "survey selection:", survey_selection_flag.shape, np.sum(survey_selection_flag)
	print "inner fiber selection:", inner_fiber_selection_flag.shape, np.sum(inner_fiber_selection_flag)

	fiber_selection_flag = survey_selection_flag
	fiber_selection_flag[survey_selection_flag] = inner_fiber_selection_flag
	print "fiber selection:", fiber_selection_flag.shape, np.sum(fiber_selection_flag)

	# Write result to data bank
	with archive.archive(data_bank, 'a') as ar:
		ar['/gal/tile_id'] = all_gals.tile_ID
		ar['/gal/visit_id'] = all_gals.visit_ID
		ar['/gal/fiber_id'] = all_gals.fiber_ID
		ar['/gal/fiber_selection_flag'] = fiber_selection_flag



	# =====================================================
	# Diagnostics
	# =====================================================
	print '... diagnostic plots'

	# fiberid distribution
	try:
		check_fiber_sel = numpy.where( fiber_selection_flag  == True)[0]
		fiber_id_sub = all_gals.fiber_ID[check_fiber_sel]
		width, center, bin_edges, hist = Utilities.make_histogram(fiber_id_sub, int(float(Ngal)/10.))
		Plotters.plot_fa_hist_fiberid(center, width, hist,
									   fig_number=0, plot_number=0,  title='n(fiberid)', base_directory='.', show_plot=True)
	except (RuntimeError, TypeError, NameError):
		print 'bad plot'
		pass


	# galaxy positions for those that are allocated (ra/dec position plot)
	try:
		Plotters.plot_fa_gal_position(all_gals.x[gal_alloc],all_gals.y[gal_alloc],
									 fig_number=1, plot_number=0,  title='galaxy position', base_directory='.', show_plot=True)
	except (RuntimeError, TypeError, NameError):
		print 'bad plot'
		pass

	# fraction allocated vs. tile number
	try:
		tile_list = np.arange(np.size(fraction_allocated))
		Plotters.plot_fa_completeness(tile_list, fraction_allocated,
									   fig_number=2, plot_number=0,  title='fiber completeness', base_directory='.', show_plot=False)
	except (RuntimeError, TypeError, NameError):
		print 'bad plot'
		pass

	#plot_name = "galaxies_with_fiber"
	#diagplot.plot_positions(all_gals.x[gal_alloc],all_gals.y[gal_alloc],plot_name)

	#plot_name = "completeness fraction"
	#tile_list = np.arange(np.size(fraction_allocated))
	#diagplot.plot_fraction(tile_list,fraction_allocated, plot_name)




# =====================================================
if __name__ == '__main__':
	main()


