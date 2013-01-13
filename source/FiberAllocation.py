# <codecell>
#!/usr/env/python

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

        self.Nfibers = self.x.size    
        self.ID = np.arange(self.Nfibers)
        self.allocated_galaxy_ID = np.zeros(self.Nfibers)
        self.allocated_galaxy_ID[:] = -1 


class MockGalaxyCatalog(object):
    """
    Generates a mock galaxy catalog over a circular footprint.
    """
    def __init__(self,
                 Ngalaxies=100,
                 radius_fov = 1.0,
                 priority_levels = 3, 
                 center_x = 0.0,
                 center_y = 0.0
                ):
        self.radius_fov = radius_fov
        self.Ngalaxies = Ngalaxies
        self.priority_levels = priority_levels
        self.center_x = center_x
        self.center_y = center_y
        self.tile_visit_ID = np.zeros(self.Ngalaxies)
        self.tile_visit_ID[:] = -1
        
        self.make_random_mock_catalog()
        
        
    def make_random_mock_catalog(self):
        radius_fov_rand = self.radius_fov * np.sqrt(np.random.rand(self.Ngalaxies)) # Radius in degrees
        theta_fov_rand = np.random.rand(self.Ngalaxies) * 2.0 * np.pi
         
        self.ID = np.arange(self.Ngalaxies)
        self.x = radius_fov_rand * np.cos(theta_fov_rand)
        self.y = radius_fov_rand * np.sin(theta_fov_rand)
        
        self.x[:] = self.x[:] + self.center_x
        self.y[:] = self.y[:] + self.center_y
        
        # random asignment with priority levels for each galaxy, higher values correspond to higher priority
        # the number of priority levels is set by the variable priority_levels
        self.priority = np.random.randn(self.Ngalaxies) 
        self.priority = 0.99*(self.priority - np.amin(self.priority))/(np.amax(self.priority) - np.amin(self.priority))
        self.priority = np.int_(self.priority * self.priority_levels)



def close_match_xy(tx, ty, fx, fy, t_id, f_id, epsilon=0.1):
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
    
    #loop over the number of fibers to find the number of targets within radius epsilon
    for i in range(n_f):
        f_x = fx[i]
        f_y = fy[i]
        n_target = tx.size
        initial_id = np.copy(t_id)
        
        
        #select the index range for the good targets
        delta_radius = np.sqrt((tx - f_x)**2 + (ty - f_y)**2)
        index = np.where(delta_radius<epsilon)
        
        #rank the targets in distance to the fiber and return the sorted indices
        id_close = initial_id[index]
        close_delta_radius = delta_radius[index]
        index_close = index[0]
    #    print 'index_close', index_close
            
        sorted_delta_radius = np.argsort(close_delta_radius)
        sorted_id = id_close[sorted_delta_radius]
        sorted_index = index_close[sorted_delta_radius]
        if(sorted_id.size):
            id_match_f[i] = f_id[i] # this is the fiber looking for its target
            id_match_t[i] = sorted_id[0] # this is the ID of the galaxy it can observe
            index_match_t[i] = sorted_index[0] # this is the index of the target in the original input array
            distance_t_to_f[i] = close_delta_radius[sorted_delta_radius][0]
#            print close_delta_radius[sorted_delta_radius][0]
        else:
            id_match_f[i] = f_id[i]
            id_match_t[i] = -1
            index_match_t[i] = -1
            distance_t_to_f[i] = -1.0
        #target_in_f[i] = sorted_id[0]
    return id_match_f, id_match_t, index_match_t, distance_t_to_f

# <codecell>

def reset_fiber_collisions(id_t, distance_f_to_t):
    """
    Fibers collisions are found if there is an ID in id_t (different from -1) that is duplicated.
    It means that the two different fibers are allocated to the same target.
    
    This function takes as an input the array of targets IDs assigned to fibers, the corresponding distance fiber_to_target.
    The function resets the id_t to -1 for the fibers with collisions that are the farthest from the target
    """
    
    #find the cases of colliding fibers
    u = np.empty((0))
    test = np.copy(id_t)
    for i in range(test.size):
        less_test = np.delete(test,i)
        x = test[i]
        if (x in less_test) and (x!=-1) and (x not in u):
            u = np.append(u,x)
    #do something if there is at least one case of colliding fibers
    if(u.size):
        n_colliding  = u.size
        print 'There are', n_colliding, 'colliding fibers to reset'
        for targets_id in u:
            #selects the fibers that are colliding
            collision_id = np.where(id_t==targets_id) 
            collision_id = collision_id[0]
            #find the distance from the fiber to the target
            collision_distance = distance_f_to_t[collision_id]
            #resets all the fibers except the closest one (there might be more than two fibers colliding)
            id_to_reset = np.delete(collision_id, np.array(np.argmin(collision_distance)))
            id_t[id_to_reset] = -1
    else:
        print 'The are zero colliding fibers'
            

# <codecell>

def make_fiber_allocation(fibers, gals, tile_visit_ID=1):
    max_n_iterations = 5
    
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
        old_n_to_allocate = n_to_allocate
       # print 'There are ', n_to_allocate, 'fibers to allocate'
       # print 'There are ', n_to_find, 'galaxies to find'

        # Initialize the arrays with the fibers that have to be allocated
        to_allocate = np.where(fibers.allocated_galaxy_ID==-1)
        to_allocate = to_allocate[0]
        fibers_x = fibers.x[to_allocate]
        fibers_y = fibers.y[to_allocate]
        fibers_ID = fibers.ID[to_allocate]
        
        #initialize the arrays with galaxies that have not been observed
        free_gals = np.where(gals.tile_visit_ID==-1)
        free_gals = free_gals[0]
        gals_x = gals.x[free_gals]
        gals_y = gals.y[free_gals]
        gals_ID = gals.ID[free_gals]
        
        # Make a pass of allocation 
        alloc_id_fibers, id_matching_gals, index_matching_gals, distance_fiber_to_gal = \
        close_match_xy(gals_x, gals_y, fibers_x, fibers_y, gals_ID, fibers_ID, epsilon=fibers.fiber_pitch)
        # Reset the colliding fibers
        reset_fiber_collisions(id_matching_gals, distance_fiber_to_gal)

        #update the ID of the galaxy to observe
        j=0
        n_sin_match = np.size(np.where(id_matching_gals==-1))

        for i in range(np.size(to_allocate)):
            fibers.allocated_galaxy_ID[to_allocate[i]] = id_matching_gals[i]
            if(id_matching_gals[i]!=-1):
                gal_assigned = gal_assigned +1
                gals.tile_visit_ID[free_gals[index_matching_gals[i]]] = tile_visit_ID
             
        
        #recalculate the numbers that still have to be allocated
        n_to_allocate = np.size(np.where(fibers.allocated_galaxy_ID==-1))
        n_to_find = np.size(np.where(gals.tile_visit_ID==-1))
        iteration = iteration + 1
    #simple test: the number of galaxies in fibers must be equal to the number of assigned fibers
    galaxies_with_fiber = np.size(np.where(gals.tile_visit_ID==tile_visit_ID))
    fibers_with_galaxy = np.size(np.where(fibers.allocated_galaxy_ID!=-1))
    if(galaxies_with_fiber!=fibers_with_galaxy):
        print 'The number of galaxies observed in this tile is not the same as the number of allocated fibers'
        sys.exit()




# Example

# Creates a set of fibers
fibers = FiberSet(Ndiameter=73, fiber_pitch=0.1, fiber_size=0.01)

# Creates a set of galaxies
gals = MockGalaxyCatalog(Ngalaxies=50000, radius_fov=fibers.radius, priority_levels=3)

# Runs the fiber allocation
make_fiber_allocation(fibers, gals, tile_visit_ID=340)

print np.size(np.where(fibers.allocated_galaxy_ID!=-1)), 'fibers will be used'
