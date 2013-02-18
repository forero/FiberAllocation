import numpy as np
import matplotlib.pyplot as plt
import sys
import AllocationModule as AM


plate_scale = 0.27/60. / (15. / 1e3) / 60.0# convert arcsec to arcmin; micron to mm; arcmin to deg --> arcmin/mm --> deg/mm
fiber_pitch_mm = 6.0
fiber_size_mm = 0.7 * 2.0

#Generate the fibers
fibers = AM.FiberSet(Ndiameter=73, fiber_pitch=fiber_pitch_mm * plate_scale, fiber_size=fiber_size_mm * plate_scale)

#Generates a mock catalog
fiber_number_density = fibers.Nfibers/((3.0 * np.sqrt(3.0)/2.0)*(fibers.radius **2))
print fiber_number_density, fibers.radius
gals = AM.MockGalaxyCatalog(number_density=2.0*fiber_number_density, radius_fov=fibers.radius, random=True)



#runs the allocation
method="galaxy_density"
for i in range(2):
    fibers = AM.FiberSet(Ndiameter=73, fiber_pitch=fiber_pitch_mm * plate_scale, fiber_size=fiber_size_mm * plate_scale)
    AM.make_fiber_allocation(fibers, gals, tile_visit_ID=i, rank_criterion=method, patrol_radius=fibers.fiber_pitch,\
                    exclusion_radius=fibers.fiber_size * 0.5)

#prints some statistics
allocated_gals = np.where(gals.tile_visit_ID!=-1)
n_allocated = np.size(allocated_gals)
n_total_gals = 1.0*gals.Ngalaxies* (3.0 * np.sqrt(3.0)/2.0/np.pi)
print 'fraction of allocated galaxies', n_allocated/n_total_gals

#makes a nice plot of allocated galaxies
fig = plt.figure(1, figsize=(9.5,9.0))
ax = plt.axes()
ax.set_xlabel("$\mathrm{ra\ [deg]}$",fontsize=25)
ax.set_ylabel("$\mathrm{dec\ [deg]}$",fontsize=25)
ax.set_title("$\mathrm{allocated\ galaxies}$", fontsize=25)
ticklabels_x = ax.get_xticklabels()
ticklabels_y = ax.get_yticklabels()

for label_x in ticklabels_x:
    label_x.set_fontsize(18)
    label_x.set_family('serif')
for label_y in ticklabels_y:
    label_y.set_fontsize(18)
    label_y.set_family('serif')

#this is the key part where the allocated galaxies are selected
allocated_gals = np.where(gals.tile_visit_ID!=-1)
plt.scatter(gals.x[allocated_gals], gals.y[allocated_gals], alpha=1.0, s=1.0)

ax.set_xlim([-1.1,1.1])
ax.set_ylim([-1.1,1.1])

ax.set_aspect('equal', 'datalim')

filename='Allocated_galaxies'

plt.savefig(filename + '.eps',format = 'eps', transparent=True)
plt.savefig(filename + '.pdf',format = 'pdf', transparent=True)
