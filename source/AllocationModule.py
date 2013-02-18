# -*- coding: utf-8 -*-

import numpy as np
import sys

class FiberSet(object):
    """
    Initialize the fibers in a hexagonal tile.
    Flat sky and cartesian positions are assumed

    Inputs:
    Ndiameter: number of fiabers along the tile's diameter. 
    Fiber_pitch: fiber_pitch in degrees
    Fiber_size: fiber_size in degrees
    center_x: fiber center in degrees
    center_y: fiber center in degrees

    Outputs:
    An object self with:
    - self.fiber_size in degrees
    - self.fiber_pitch in degrees
    - self.radius: tile radius in degrees
    - self.x, self.y: fiber positions in degrees
    - self.ID: unique fiber ID in the tile
    - self.Nfibers: total number of fibers in the tile
    - self.allocated_galaxy_ID: ID of the galaxy to be observed.
    """
    def __init__(self,
                 Ndiameter=73,
                 fiber_pitch = 6.0,
                 fiber_size = 0.05,
                 center_x = 0.0,
                 center_y = 0.0,
                ):
        self.center_x = center_x
        self.center_y = center_y
        self.fiber_size = fiber_size
        self.fiber_pitch = fiber_pitch # This quantity has to be modif
        self.Ndiameter = Ndiameter
        self.make_hexagon_tile()
            
    def make_hexagon_tile(self):
        self.Ndiameter = self.Ndiameter - (self.Ndiameter%2 -1) #makes the number of spines in the main axis an odd number
        self.radius = (self.Ndiameter - 1.0) * self.fiber_pitch / 2.0
                
        x_baseline = np.arange(self.Ndiameter) * self.fiber_pitch
        x_baseline = x_baseline - np.median(x_baseline)
        y_baseline = np.zeros(self.Ndiameter, dtype=float)

        self.x = np.empty((0))
        self.y = np.empty((0))

        self.x = np.append(self.x, x_baseline)
        self.y = np.append(self.y, y_baseline)

        # first pass (dec positive direction)
        n_vertical =  ((self.Ndiameter-1) - ((self.Ndiameter-1)/2))
        
        for i in range(n_vertical):
            new_n_points = self.Ndiameter - i -1
            new_x_line = np.zeros((new_n_points), dtype=float)
            new_y_line = np.zeros((new_n_points), dtype=float)
            new_x_line = x_baseline[0:self.Ndiameter-i-1] + (i+1) * self.fiber_pitch * 0.5
            new_y_line = y_baseline[0:self.Ndiameter-i-1] + (i+1) * (np.sqrt(3.0)/2.0) * self.fiber_pitch
            self.x = np.append(self.x, new_x_line)
            self.y = np.append(self.y, new_y_line)

        # second pass (dec negative direction)
        for i in range(n_vertical):
            new_n_points = self.Ndiameter - i -1
            new_x_line = np.zeros((new_n_points), dtype=float)
            new_y_line = np.zeros((new_n_points), dtype=float)
            new_x_line = x_baseline[0:self.Ndiameter-i-1] + (i+1) * self.fiber_pitch * 0.5
            new_y_line = y_baseline[0:self.Ndiameter-i-1] - (i+1) * (np.sqrt(3.0)/2.0) * self.fiber_pitch
            self.x = np.append(self.x, new_x_line)
            self.y = np.append(self.y, new_y_line)
        
        self.x[:] = self.x[:] + self.center_x
        self.y[:] = self.y[:] + self.center_y

        self.original_x = self.x.copy()
        self.original_y = self.y.copy()
        
        
        self.Nfibers = self.x.size    
        self.ID = np.arange(self.Nfibers)
        self.allocated_galaxy_ID = np.zeros(self.Nfibers)
        self.allocated_galaxy_ID[:] = -1 

# <codecell>

class MockGalaxyCatalog(object):
    """
    Generates a mock galaxy catalog over a circular footprint.
    """
    def __init__(self,
                 radius_fov = 1.0,
                 number_density = 1000,
                 galaxy_type_abundance = np.array([1.0]), 
                 center_x = 0.0,
                 center_y = 0.0, 
                 filename = "data.dat",
                 x_in = np.array([0.0]),
                 y_in = np.array([0.0]),
                 id_in = np.array([0]),
                 random = False
                ):        
        #if there is a filename to load do it
        
        if(filename!="data.dat"):
            self.filename = filename
            self.load_galaxy_catalog()
        elif(random==True):
            #basic intrisic properties of the galaxy distribution
            self.radius_fov = radius_fov
            self.Ngalaxies = np.int_(3.14159 * (self.radius_fov**2) * number_density)
        
            self.galaxy_type_abundance = galaxy_type_abundance.copy()
            self.n_types = np.size(self.galaxy_type_abundance)
            
            self.center_x = center_x
            self.center_y = center_y
            self.make_random_mock_catalog()
        else:
            self.Ngalaxies = np.size(x_in)
            self.x = x_in
            self.y = y_in
            self.ID = id_in
            self.center_x = np.mean(x_in)
            self.center_y = np.mean(y_in)

        #properties related the fiber allocation
        self.tile_visit_ID = np.zeros(self.Ngalaxies)
        self.tile_visit_ID[:] = -1
        self.fibers_per_galaxy = np.zeros(self.Ngalaxies)
        self.galaxies_around_fiber_size = np.zeros(self.Ngalaxies)
        self.galaxies_around_fiber_pitch = np.zeros(self.Ngalaxies)
            
        
    def load_galaxy_catalog(self):
        print "Loading", self.filename
        raw_data = np.loadtxt(self.filename)
        selected = raw_data[:,7]

        
        self.x = raw_data[:,8]
        self.y = raw_data[:,9]
        self.ID = raw_data[:,1]
        self.Ngalaxies = np.size(self.ID)
        
        self.gal_type = raw_data[:,0]
        
        self.center_x = np.average(self.x)
        self.center_y = np.average(self.y)
        self.radius_fov = -1.0
        print np.shape(self.x)
        
    def make_random_mock_catalog(self):
        radius_fov_rand = self.radius_fov * np.sqrt(np.random.rand(self.Ngalaxies)) # Radius in degrees
        theta_fov_rand = np.random.rand(self.Ngalaxies) * 2.0 * np. pi
         
        self.ID = np.arange(self.Ngalaxies)
        self.x = radius_fov_rand * np.cos(theta_fov_rand)
        self.y = radius_fov_rand * np.sin(theta_fov_rand)
        
        self.x[:] = self.x[:] + self.center_x
        self.y[:] = self.y[:] + self.center_y

def count_fibers_per_galaxy(galaxies_x, galaxies_y,  fibers_x, fibers_y, fiber_pitch, fibers_per_galaxy):
    n_items = galaxies_x.size
    for i in range(n_items):
        x = galaxies_x[i]
        y = galaxies_y[i]
        delta_fiber_radius = np.sqrt((fibers_x - x)**2 + (fibers_y - y)**2)
        index = np.where(delta_fiber_radius < fiber_pitch)
        fibers_per_galaxy[i] = np.size(index)
    return fibers_per_galaxy

# <codecell>

def count_galaxies_per_galaxy(galaxies_x, galaxies_y, min_radius, max_radius, galaxies_per_galaxy):
    n_items = galaxies_x.size
    for i in range(n_items):
        x = galaxies_x[i]
        y = galaxies_y[i]
        delta_radius = np.sqrt((galaxies_x - x)**2 + (galaxies_y - y)**2)
        index = np.where((delta_radius > min_radius) & (delta_radius<max_radius))
        galaxies_per_galaxy[i] = np.size(index)
    return galaxies_per_galaxy

# <codecell>

def close_match_xy(tx, ty, fx, fy, t_id, f_id, f_per_target, t_per_target, t_close, exclusion_epsilon=0.01,epsilon=0.1, rank_criterion="distance"):
    """
    Finds close matches between two sets of points in euclidian space: (tx,ty) and (fx,fy) (targets and fibers).
    Each target and fiber has a well defined id (t_id, f_id)
    A close match is defined below a level of epsilon in delta_x and delta_y.
    """
    n_f = fx.size
    id_match_f = np.arange(n_f)
    id_match_t = np.arange(n_f)
    index_match_t = np.arange(n_f)
    distance_t_to_f = np.zeros(n_f)
    
    #find the numer of fibers per target
    f_per_target = count_fibers_per_galaxy(tx, ty, fx, fy, epsilon, f_per_target)
    t_per_target = count_galaxies_per_galaxy(tx, ty, exclusion_epsilon, epsilon, t_per_target)
    t_close = count_galaxies_per_galaxy(tx, ty, 0.0, epsilon * 0.1, t_close)
    
    #loop over the number of fibers to find the number of targets within radius epsilon
    for i in range(n_f):
        f_x = fx[i]
        f_y = fy[i]
        n_target = tx.size
        initial_id = np.copy(t_id)
        
        
        #find the index of the targets around epsilon of the fiber
        delta_radius = np.sqrt((tx - f_x)**2 + (ty - f_y)**2)
        index = np.where(delta_radius<epsilon)
        
        #cut the targets to those close to the fibers
        id_close = initial_id[index]
        close_tx = tx[index]
        close_ty = ty[index]
        close_delta_radius = delta_radius[index]
        fiber_per_target = f_per_target[index]
        galaxy_density = t_per_target[index]
        galaxy_in_fiber = t_close[index]
        index_close = index[0]
    
        # sort the IDs by some criterion
        if(rank_criterion=="distance"):
            sorted_delta_radius = np.argsort(close_delta_radius)
        elif(rank_criterion=="fibers_per_target"):
            sorted_delta_radius = np.argsort(fiber_per_target)
        elif(rank_criterion=="galaxy_density"):
            # hight prob indicates a high number density of targets, meaning that it is hard to get it with out collisions.
            prob = np.maximum(galaxy_density/7.0, galaxy_in_fiber) 
            priority = np.exp(-prob) 
            sorted_delta_radius = np.argsort(priority)
        else:
            print "The priority criterion does not exist:", rank_criterion
            sys.exit()
            
            
        sorted_id = id_close[sorted_delta_radius]
        sorted_index = index_close[sorted_delta_radius]
        sorted_tx = close_tx[sorted_delta_radius]
        sorted_ty = close_ty[sorted_delta_radius]
        if(sorted_id.size):
            id_match_f[i] = f_id[i] # this is the fiber looking for its target
            id_match_t[i] = sorted_id[0] # this is the ID of the galaxy it can observe
            index_match_t[i] = sorted_index[0] # this is the index of the target in the original input array
            distance_t_to_f[i] = close_delta_radius[sorted_delta_radius][0]
            fx[i] = sorted_tx[0]
            fy[i] = sorted_ty[0]
#            print close_delta_radius[sorted_delta_radius][0]
        else:
            id_match_f[i] = f_id[i]
            id_match_t[i] = -1
            index_match_t[i] = -1
            distance_t_to_f[i] = -1.0
        #target_in_f[i] = sorted_id[0]
    return id_match_f, id_match_t, index_match_t, distance_t_to_f

# <codecell>

def reset_fiber_collisions(id_t, f_x, f_y, f_x0, f_y0, epsilon=0.0):
    """
    Fiber collisions are found if the distance between two fibers is less than epsilon
    f_x0, f_y0 are the unperturbed fiber positions.
    """
    
    #find the cases of colliding fibers
    u = np.empty((0))
    test = np.copy(id_t)
    n_reset = 0
    for i in range(test.size):
        fx = f_x[i]
        fy = f_y[i]
        delta_radius = np.sqrt((fx - f_x)**2 + (fy - f_y)**2)
        index = np.where(delta_radius < epsilon)
        if np.size(index)>1: #the fiber 'i' always collides with itself
            id_delete = np.argmin(delta_radius[index])
            index_collide = np.delete(index,id_delete)

            #resets the ID of the targets and moves the fibers back to its original position
            id_t[index_collide] = -1.0
            f_x[i] = f_x0[i]
            f_y[i] = f_y0[i]
            n_reset = n_reset+1
    print n_reset, 'fibers were reset'
    return n_reset
            

# <codecell>

def make_fiber_allocation(fibers, gals, tile_visit_ID=1, patrol_radius=0.0, rank_criterion="distance", exclusion_radius=0.0):
    max_n_iterations = 10

    #check that there are fibers to be allocated
    n_to_allocate = np.size(np.where(fibers.allocated_galaxy_ID==-1))
    n_to_find = np.size(np.where(gals.tile_visit_ID==-1))
    iteration = 0
    gal_assigned = 0
    old_n_to_allocate = 0
    if(n_to_allocate != fibers.Nfibers):
        print 'All the fibers must be free before starting its allocation'
        sys.exit()
    while (n_to_allocate > 0 and n_to_find > 0 and iteration < max_n_iterations and old_n_to_allocate!=n_to_allocate):
       # print 'There are ', n_to_allocate, 'fibers to allocate'
       # print 'There are ', n_to_find, 'galaxies to find'
        old_n_to_allocate = n_to_allocate

        # Initialize the arrays with the fibers that have to be allocated
        to_allocate = np.where(fibers.allocated_galaxy_ID==-1)
        to_allocate = to_allocate[0]
        fibers_x = fibers.x[to_allocate]
        fibers_y = fibers.y[to_allocate]
        fibers_x0 = fibers.original_x[to_allocate]
        fibers_y0 = fibers.original_y[to_allocate]
        fibers_ID = fibers.ID[to_allocate]
        
        #initialize the arrays with galaxies that have not been observed
        free_gals = np.where(gals.tile_visit_ID==-1)
        free_gals = free_gals[0]
        gals_x = gals.x[free_gals]
        gals_y = gals.y[free_gals]
        gals_ID = gals.ID[free_gals]
        fibers_per_galaxy = gals.fibers_per_galaxy[free_gals]
        galaxies_per_galaxy = gals.galaxies_around_fiber_pitch[free_gals]
        galaxies_per_fiber_size = gals.galaxies_around_fiber_size[free_gals]
        
        
        # Make a pass of allocation 
        alloc_id_fibers, id_matching_gals, index_matching_gals, distance_fiber_to_gal = \
        close_match_xy(gals_x, gals_y, fibers_x, fibers_y, gals_ID, fibers_ID, \
        fibers_per_galaxy, galaxies_per_galaxy, galaxies_per_fiber_size, \
        exclusion_epsilon = exclusion_radius, epsilon=patrol_radius, rank_criterion=rank_criterion)
        
        # Reset the colliding fibers
        n_reset = reset_fiber_collisions(id_matching_gals, fibers_x, fibers_y, fibers_x0, fibers_y0, epsilon=exclusion_radius)

        #update the ID of the galaxy to observe
        j=0
        n_sin_match = np.size(np.where(id_matching_gals==-1))

        #loop over the fibers to assign them the galaxy and find the galaxy to assign it its corresponding fiber and tile ID
        for i in range(np.size(to_allocate)):
            fibers.allocated_galaxy_ID[to_allocate[i]] = id_matching_gals[i]
            fibers.x[to_allocate[i]] = fibers_x[i]
            fibers.y[to_allocate[i]] = fibers_y[i]
            if(id_matching_gals[i]!=-1):
                gal_assigned = gal_assigned +1
                gals.tile_visit_ID[free_gals[index_matching_gals[i]]] = tile_visit_ID
             
        
        #recalculate the number of fibers that need to be allocated.
        n_to_allocate = np.size(np.where(fibers.allocated_galaxy_ID==-1))
        n_to_find = np.size(np.where(gals.tile_visit_ID==-1))
        iteration = iteration + 1

        if(n_reset==0):
            return
    return
