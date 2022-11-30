"""
Honeycomb Lattice Generator

Author: Gokhan Oztarhan
Created date: 09/06/2019
Last modified: 21/07/2022
"""

import time
import logging

import numpy as np


logger = logging.getLogger(__name__)


def honeycomb(
    flk_type, a, n_side, width=2, bc='', com_to_origin=False, **kwargs
):
    """
    Main function for honeycomb lattice generator.
    
    Parameters
    ----------
    flk_type : str
        flake type
    a : float, 
        dot-to-dot distance,
        lattice constant
    n_side : int
        number of points on one side of the flk_type
    width : int, optional
        vertical width of nanoribbon flk_type,
        vertical number of points
    bc : str, optional
        boundary condition for nanoribbon flk_type
    com_to_origin : bool, optional
        shift center of mass to origin

    Returns
    -------
    pos[:,:] : float
        positions
    sub[:] : int
        sub lattice type,
        0 or 1 depending on A or B sublattice
    m[:], n[:] : int 
        m and n types of honeycomb lattice
    ind_A[:], ind_B[:] : int
        indices of A and B sublattice
    distances[:,:] : int
        distances between all pairs of lattice sites
    ind_NN[:,:] : int
        indices of nearest neighbors,
        only for upper triangle indices
    ind_NN_2nd[:,:] : int
        indices of 2nd nearest neighbors,
        only for upper triangle indices
    n_side: int
        number of points on one side of the flk_type
    """   
    _FLAKE_FUNC = {
        'hexagonal_zigzag': hexagonal_zigzag,
        'hexagonal_armchair': hexagonal_armchair,
        'triangular_zigzag': triangular_zigzag,
        'triangular_armchair': triangular_armchair,
        'nanoribbon': nanoribbon,
    }  

    ind, pos, sub, m, n, \
    n_total, n_NN, n_side, cpu_time = _FLAKE_FUNC[flk_type](
        a, n_side, width=width, **kwargs
    )
                              
    # Select interior elements only
    pos, sub, m, n = pos[ind], sub[ind], m[ind], n[ind]
    
    # Sublattice indices
    ind_A, ind_B = np.where(sub == 1)[0], np.where(sub == 0)[0]
    
    logger.info( 
        'lat_type = honeycomb (%.3f s)\n' %(cpu_time) \
        + 'a = %.5e\n' %(a) \
        + 'flk_type = %s\n' %(flk_type) \
        + 'n_side = %i\n' %(n_side) \
        + 'n_total calc/found = %i/%i\n' %(n_total, pos.shape[0])   
    )  
    
    # Distances between all pair of sites
    distances, cpu_time = _distances(pos)
    logger.info('Distances calculated. (%.3f s)\n' %cpu_time)
    
    # Get nearest neighbor indices
    ind_NN, cpu_time = _nearest_neighbor(a, distances) 
    
    logger.info(
        '1st NN calc/found = %i/%i (%.3f s)\n' \
        %(n_NN, 2 * ind_NN.shape[0], cpu_time)
    ) 
    
    # Get second nearest neighbor indices
    ind_NN_2nd, cpu_time = _nearest_neighbor_2nd(a, distances)
    
    logger.info(
        '2nd NN found = %i (%.3f s)\n' \
        %(2 * ind_NN_2nd.shape[0], cpu_time)
    ) 
    
    # Apply boundary conditions for nanoribbon
    if flk_type == 'nanoribbon':
        ind_drop, ind_NN, ind_NN_2nd, cpu_time = \
            _boundary_condition(a, bc, pos, ind_NN, ind_NN_2nd)
        
        if ind_drop is not None:
            # Drop excluded elements
            ind = np.arange(pos.shape[0])
            ind = np.delete(ind, ind_drop)
            pos, sub, m, n = pos[ind], sub[ind], m[ind], n[ind] 
            
            # Correct sublattice indices
            ind_A, ind_B = np.where(sub == 1)[0], np.where(sub == 0)[0]
            
            # Shift elements by keeping track of dropped indices
            for i, ind in enumerate(ind_drop):
                ind_NN[ind_NN[:,0] > ind - i,0] -= 1
                ind_NN[ind_NN[:,1] > ind - i,1] -= 1
                ind_NN_2nd[ind_NN_2nd[:,0] > ind - i,0] -= 1
                ind_NN_2nd[ind_NN_2nd[:,1] > ind - i,1] -= 1
            
            # Drop duplicate indices
            ind_NN = np.unique(ind_NN, axis=0)
            ind_NN_2nd = np.unique(ind_NN_2nd, axis=0)

        logger.info(
            'width = %i\n' %(width) \
            + 'bc = %s\n' %(bc) \
            + 'n_total after boundary conditions = %i\n' \
                %(pos.shape[0])
            + '1st NN after boundary conditions = %i\n' \
                %(2 * ind_NN.shape[0])
            + '2nd NN after boundary conditions = %i\n' \
                %(2 * ind_NN_2nd.shape[0])
        ) 
        
    # Shift lattice center to origin
    pos = _shift_to_origin(pos)
    
    # Shift lattice center of mass to origin
    if com_to_origin:
        pos = _shift_center_of_mass_to_origin(pos)
    
    return pos, sub, m, n, ind_A, ind_B, distances, ind_NN, ind_NN_2nd, n_side

#------------------------------------------------------------------------------

def _generate(a, mn_range, vertical):
    """
    Main generator of honeycomb lattice.
    
    Lattice vectors with nearest neighbor distance a:
    r1 and r2
    
    Generation of lattice:
    posA = m*r1 + n*r2
    posB = m*r1 + n*r2 + r3
    
    where m and n are integers.
    """
    r1x, r1y = a*np.sqrt(3)/2, a*3/2
    r2x, r2y = -a*np.sqrt(3)/2, a*3/2
    r3x, r3y = 0, a
    shiftx, shifty = -a*np.sqrt(3)/2, -a/2
    
    if (vertical != 0):
       r1x, r1y = r1y, r1x
       r2x, r2y = r2y, r2x
       r3x, r3y = r3y, r3x
       shiftx, shifty = shifty, shiftx
        
    mm = np.arange(-mn_range,mn_range+1,dtype=int)
    
    m, n = np.meshgrid(mm,mm)
    m = m.flatten()[:,None]
    n = n.flatten()[:,None]
        
    posAx = m*r1x + n*r2x + shiftx
    posAy = m*r1y + n*r2y + shifty
    posBx = posAx + r3x
    posBy = posAy + r3y
    
    ones = np.full(m.shape[0],1)
    zeros = np.zeros(m.shape[0])
    
    pos = np.vstack([np.hstack([posAx,posAy]),
                     np.hstack([posBx,posBy])]) 
    sub = np.hstack([ones,zeros])
    m = np.vstack([m,m])   
    n = np.vstack([n,n])    
    
    return pos, sub, m[:,0], n[:,0]
    
    
def _shift_to_origin(pos):
    """Shift flake center to origin."""
    ymax = np.max(pos[:,1])
    ymin = np.min(pos[:,1])
    yshift = (ymax + ymin) / 2.0
    pos[:,1] = pos[:,1] - yshift
    
    xmax = np.max(pos[:,0])
    xmin = np.min(pos[:,0])
    xshift = (xmax + xmin) / 2.0
    pos[:,0] = pos[:,0] - xshift
    
    return pos
    
    
def _shift_center_of_mass_to_origin(pos):
    """Shift flake center of mass to origin."""
    # Center of mass of the lattice geometry
    com = pos.mean(axis=0)
    
    # Shift center of mass of the lattice to origin
    pos[:,0] -= com[0]
    pos[:,1] -= com[1]
            
    return pos
    
  
def _flake_vertex(r, theta, xshift=0, yshift=0):
    """Generate flake vertices."""
    vertex = np.zeros([np.size(theta),2])
    vertex[:,0] = r*np.cos(theta) + xshift
    vertex[:,1] = r*np.sin(theta) + yshift
   
    return vertex
    

def _interior_elements_convex(pos, vertex):
    """
    Indices of interior elements for a convex polygon.
    Vertex positions should be in counter-clockwise order.
    """
    x = vertex[:,0]
    y = vertex[:,1] 
    
    # elimination on the index array is faster than the whole array
    ind = np.arange(pos.shape[0])
    
    # all points are left to the line segment if checker > 0
    for i in range(vertex.shape[0]):
        checker = (pos[ind,1] - y[i-1]) * (x[i] - x[i-1]) \
                  - (pos[ind,0] - x[i-1]) * (y[i] - y[i-1])               
        ind = ind[checker > 0]
        
    return ind
    
#------------------------------------------------------------------------------

def _distances(pos):
    """Distance calculator between all pairs of lattice sites."""
    tic = time.time()
    
    distances = np.zeros((pos.shape[0], pos.shape[0]))
    
    # Looping over indices is both memory friendly and faster.
    for i in range(pos.shape[0]):
        distances[i,:] = np.sqrt(
            (pos[i,0] - pos[:,0])**2 \
            + (pos[i,1] - pos[:,1])**2
        )

    toc = time.time()
    return distances, toc - tic
    
    
def _nearest_neighbor(a, distances):
    """1st Nearest neighbor finder."""
    tic = time.time()
    
    eps = a * 0.001
    ind = np.arange(distances.shape[0])
    ind_NN = np.zeros((3 * distances.shape[0], 2), dtype=int)
    
    # Looping over indices is suprisingly faster than the vectorized
    # calculation over all matrix. Here, the loop is over only on
    # the upper triangular part of the matrix since it is symmetrical.
    j = 0
    for i in ind[:-1]:
        ind_nn = (distances[i,i+1:] < a + eps).nonzero()[0]
        if ind_nn.shape[0] > 0:
            ind_nn += i + 1
            
            ind_NN[j:j+ind_nn.shape[0],0] = i
            ind_NN[j:j+ind_nn.shape[0],1] = ind_nn

            j += ind_nn.shape[0]
    
    ind_NN = ind_NN[:j,:]
    
    toc = time.time()
    return ind_NN, toc - tic


def _nearest_neighbor_2nd(a, distances):
    """2nd Nearest neighbor finder."""
    tic = time.time()
    
    eps = a * 0.001
    a_sqrt_3 = a * np.sqrt(3)
    ind = np.arange(distances.shape[0])
    ind_NN_2nd = np.zeros((6 * distances.shape[0], 2), dtype=int)
    
    # Looping over indices is suprisingly faster than the vectorized
    # calculation over all matrix. Here, the loop is over only on
    # the upper triangular part of the matrix since it is symmetrical.
    j = 0
    for i in ind[:-1]:
        ind_nn = (distances[i,i+1:] < a_sqrt_3 + eps).nonzero()[0]
        ind_ind_nn = (distances[i,i+1+ind_nn] > a_sqrt_3 - eps).nonzero()[0]
        ind_nn = ind_nn[ind_ind_nn]
        if ind_nn.shape[0] > 0:
            ind_nn += i + 1
            
            ind_NN_2nd[j:j+ind_nn.shape[0],0] = i
            ind_NN_2nd[j:j+ind_nn.shape[0],1] = ind_nn

            j += ind_nn.shape[0]
    
    ind_NN_2nd = ind_NN_2nd[:j,:]
    
    toc = time.time()
    return ind_NN_2nd, toc - tic


def _boundary_condition(a, bc, pos, ind_NN, ind_NN_2nd):
    """Boundary condition generator."""
    tic = time.time()
    
    ind_drop = None
    
    if 'x' in bc: 
        xr = np.max(pos[:,0]) - a * np.sqrt(3) / 2 / 2
        xl = np.min(pos[:,0]) + a * np.sqrt(3) / 2 / 2
        indr = np.where(pos[:,0] > xr)[0] # right
        indl = np.where(pos[:,0] < xl)[0] # left
        
        # Sort with respect to y position
        xr_sorted = np.argsort(pos[indr,1])
        xl_sorted = np.argsort(pos[indl,1])
        
        # Indices on the left will be replaced by indices on the right
        # since indices on the left will be dropped from lattice.
        for i, ind_L in enumerate(indl[xl_sorted]):
            ind_R = indr[xr_sorted[i]]
            
            # 1st nearest neighbor
            ind_xl = np.where(ind_NN[:,0] == ind_L)[0]
            ind_NN[ind_xl,0] = ind_R
            ind_xl = np.where(ind_NN[:,1] == ind_L)[0]
            ind_NN[ind_xl,1] = ind_R
            
            # 2nd nearest neighbor
            ind_xl = np.where(ind_NN_2nd[:,0] == ind_L)[0]
            ind_NN_2nd[ind_xl,0] = ind_R
            ind_xl = np.where(ind_NN_2nd[:,1] == ind_L)[0]
            ind_NN_2nd[ind_xl,1] = ind_R
            
        # Sort with respect to index number
        # These are the indices being dropped
        ind_drop = np.sort(indl)
        
        # Remaining 2nd nearest neighbor
        xr2 = np.max(pos[:,0]) - 3 * a * np.sqrt(3) / 2 / 2
        xl2 = np.min(pos[:,0]) + 3 * a * np.sqrt(3) / 2 / 2
        indr = np.where((pos[:,0] > xr2) & (pos[:,0] < xr))[0]
        indl = np.where((pos[:,0] < xl2) & (pos[:,0] > xl))[0]
        
        # Sort with respect to y position
        xr_sorted = np.argsort(pos[indr,1])
        xl_sorted = np.argsort(pos[indl,1])
        
        # Append to 2nd nearest neighbors
        ind_bcx = np.hstack(
            [indr[xr_sorted][:,None], indl[xl_sorted][:,None]]
        )
        ind_NN_2nd = np.append(ind_NN_2nd, ind_bcx, 0)
    
    if 'y' in bc:
        yt = np.max(pos[:,1]) - a / 2 / 2
        yb = np.min(pos[:,1]) + a / 2 / 2
        indt = np.where(pos[:,1] > yt)[0] # top
        indb = np.where(pos[:,1] < yb)[0] # bottom
        
        # Sort with respect to x position
        yt_sorted = np.argsort(pos[indt,0])
        yb_sorted = np.argsort(pos[indb,0])
        
        # Append to 1st nearest neighbors
        ind_bcy = np.hstack(
            [indt[yt_sorted][:,None], indb[yb_sorted][:,None]]
        )
        ind_NN = np.append(ind_NN,ind_bcy,0)
        
        # 2nd nearest neighbors are the 1st nearest neighbors 
        # of the opposite side lattice sites. 
        ind_bcy = np.zeros((6 * indt.shape[0],2), dtype=int)
        j = 0
        for i, ind_T in enumerate(indt[yt_sorted]):
            ind_B = indb[yb_sorted[i]]
            
            # Find 1st nearest neighbor of top
            ind_yt = np.where(ind_NN[:,0] == ind_T)[0]
            col_T, col_T_NN = 0, 1
            if ind_yt.shape[0] < 1:
                ind_yt = np.where(ind_NN[:,1] == ind_T)[0]
                col_T, col_T_NN = 1, 0
            
            # Set 2nd nearest neighbor of bottom using 1st NN of top
            ind_bcy[j:j+ind_yt.shape[0],0] = np.full(ind_yt.shape[0], ind_B)
            ind_bcy[j:j+ind_yt.shape[0],1] = ind_NN[ind_yt,col_T_NN]
            j += ind_yt.shape[0]
            
            # Find 1st nearest neighbor of bottom
            ind_yb = np.where(ind_NN[:,0] == ind_B)[0]
            col_B, col_B_NN = 0, 1
            if ind_yb.shape[0] < 1:
                ind_yb = np.where(ind_NN[:,1] == ind_B)[0]
                col_B, col_B_NN = 1, 0
            
            # Set 2nd nearest neighbor of top using 1st NN of bottom   
            ind_bcy[j:j+ind_yb.shape[0],0] = np.full(ind_yb.shape[0], ind_T)
            ind_bcy[j:j+ind_yb.shape[0],1] = ind_NN[ind_yb,col_B_NN]
            j += ind_yb.shape[0]
        
        ind_bcy = ind_bcy[:j,:]  
        
        # Drop the same indices e.g. (3,3) or (42,42) if exist.
        # They are diagonal elements. 
        ind_bcy = np.delete(
            ind_bcy, np.where(ind_bcy[:,0] == ind_bcy[:,1])[0], axis=0
        ) 
        
        # Append to 2nd nearest neighbors
        ind_NN_2nd = np.append(ind_NN_2nd, ind_bcy, 0)
    
    toc = time.time()        
    return ind_drop, ind_NN, ind_NN_2nd, toc-tic

#------------------------------------------------------------------------------

def hexagonal_zigzag(a, n_side, **kwargs):
    """Generator for hexagonal flake with zigzag edges."""
    tic = time.time()  
    
    # Calculate the total number of dots
    n_total = (n_side**2)*6
      
    # Generate honeycomb lattice
    pos, sub, m, n = _generate(a,n_side+2,0)
    
    # Hexagon for zigzag edge
    r = (2*n_side - 1)*a*np.sqrt(3)/2 + a*np.sqrt(3)/2
    theta = np.radians(np.arange(0,360,60))
    vertex = _flake_vertex(r,theta)

    # Get only interior elements of lattice positions
    ind = _interior_elements_convex(pos,vertex)

    # 1st nearest neighbor count
    n_NN = 2*(n_side*6) + 3*(n_total-(n_side*6))   
    
    toc = time.time()            
    return ind, pos, sub, m, n, n_total, n_NN, n_side, toc-tic 
    
    
def hexagonal_armchair(a, n_side, **kwargs):
    """Generator for hexagonal flake with armchair edges."""
    tic = time.time()  
    
    if n_side < 4:
        n_side = 4
    elif n_side % 2 != 0:
        n_side = n_side + 1        
    
    # Calculate the total number of dots
    n_total = (np.sum(np.arange(n_side//2,0,-1) - 1)*6 + 1)*6  
    
    # Generate honeycomb lattice
    pos, sub, m, n = _generate(a,n_side+2,0)
    
    # Hexagon for armchair edge
    r = ( 3*(n_side//2 - 1) + 1 )*a + a/2    
    theta = np.radians(np.arange(30,361,60))
    vertex = _flake_vertex(r,theta)

    # Get only interior elements of lattice positions
    ind = _interior_elements_convex(pos,vertex)

    # 1st nearest neighbor count
    n_NN = 2*((n_side-1)*6) + 3*(n_total-((n_side-1)*6))
    
    toc = time.time()          
    return ind, pos, sub, m, n, n_total, n_NN, n_side, toc-tic 
    
    
def triangular_zigzag(a, n_side, **kwargs):
    """Generator for triangular flake with zigzag edges."""
    tic = time.time()       
    
    # Calculate the total number of dots
    n_total = np.sum(np.arange(3,3+n_side*2,2)) + n_side*2 + 1 

    # Generate honeycomb lattice
    pos, sub, m, n = _generate(a,n_side+2,0)
    
    # Triangle for zigzag edge
    r = (n_side + 1)*a + a
    theta = np.radians(np.arange(90,360,120))
    yshift = (1 - np.mod(n_side,3))*a
    vertex = _flake_vertex(r,theta,yshift=yshift)

    # Get only interior elements of lattice positions
    ind = _interior_elements_convex(pos,vertex)
    
    # Removing the elements which are out of the triangle boundaries
    ind = np.delete(ind,np.argmax(pos[ind,0]),0)
    ind = np.delete(ind,np.argmin(pos[ind,0]),0)
    ind = np.delete(ind,np.argmax(pos[ind,1]),0)

    # 1st nearest neighbor count
    n_NN = 2*(3 + n_side*3) + 3*(n_total-(3 + n_side*3))
    
    toc = time.time()  
    return ind, pos, sub, m, n, n_total, n_NN, n_side, toc-tic 
    
    
def triangular_armchair(a, n_side, **kwargs):
    """Generator triangular flake with armchair edges"""
    tic = time.time()    
    
    if n_side < 4:
        n_side = 4
    elif n_side % 2 != 0:
        n_side = n_side + 1       
        
    # Calculate the total number of dots    
    n_total = np.sum(np.arange(1,n_side//2 + 1))*6
        
    # Generate honeycomb lattice
    pos, sub, m, n = _generate(a,n_side+2,1)
    
    # Triangle for armchair edge
    r = n_side*a*np.sqrt(3)/2 + a*np.sqrt(3)/2
    theta = np.radians(np.arange(90,360,120))
    vertex = _flake_vertex(r,theta)

    # Get only interior elements of lattice positions
    ind = _interior_elements_convex(pos,vertex)
    
    # 1st nearest neighbor count
    n_NN = 2*(n_side*3) + 3*(n_total-n_side*3)
    
    toc = time.time()
    return ind, pos, sub, m, n, n_total, n_NN, n_side, toc-tic  
    
    
def nanoribbon(a, n_side, width, **kwargs):
    """Generator for nanoribbon"""
    tic = time.time()
   
    if width % 2 != 0: 
        width = width + 1 
  
    # Calculate the total number of dots  
    n_total = ((2*n_side) + 1)*width

    # Generate honeycomb lattice
    pos, sub, m, n = _generate(a,np.max([n_side,width])+2,0)
    
    # y and x sides (yt=ytop, yb=ybottom)
    yt = (width + width/2 - 1)*a/2 \
         + (1 - np.mod(width/2,2))*1.5*a + a/2
    yb = -(width + width/2 - 1)*a/2 \
         + (1 - np.mod(width/2,2))*1.5*a - a/2
    xr = n_side*a*np.sqrt(3)/2 \
         - (1 + np.mod(n_side,2))*a*np.sqrt(3)/2 + a*np.sqrt(3)/2/2
    xl = -n_side*a*np.sqrt(3)/2 \
         - (1 + np.mod(n_side,2))*a*np.sqrt(3)/2 - a*np.sqrt(3)/2/2
    
    # Removing the elements which are out of the nanoribbon boundaries
    ind = np.arange(pos.shape[0])
    ind = ind[pos[ind,0] < xr]
    ind = ind[pos[ind,0] > xl]
    ind = ind[pos[ind,1] < yt]
    ind = ind[pos[ind,1] > yb]
    
    # 1st nearest neighbor count
    n_NN = 2*(2*n_side+2*width) + 3*(n_total-(2*n_side+2*width))
        
    toc = time.time()
    return ind, pos, sub, m, n, n_total, n_NN, n_side, toc-tic 

#------------------------------------------------------------------------------

########################
# DEPRECATED FUNCTIONS #
########################

def _nearest_neighbor_deprecated(m, n, ind_A, ind_B):
    """1st Nearest neighbor finder for honeycomb lattice."""
    tic = time.time()
    
    ind_NN = np.zeros([m.shape[0]*3,2],dtype=np.int64) 
    
    j = 0
    for i in ind_A:
        mA, nA = m[i], n[i]
    
        ind_B_m = ind_B[m[ind_B] == mA] 
        ind_B_n = ind_B_m[n[ind_B_m] == nA]
        if ind_B_n:
            ind_NN[j,0], ind_NN[j,1] = i, ind_B_n[0]
            j += 1
            
        ind_B_n = ind_B_m[n[ind_B_m] == nA-1]
        if ind_B_n:
            ind_NN[j,0], ind_NN[j,1] = i, ind_B_n[0]
            j += 1
    
        ind_B_n = ind_B[n[ind_B] == nA] 
        ind_B_m = ind_B_n[m[ind_B_n] == mA-1]
        if ind_B_m:
            ind_NN[j,0], ind_NN[j,1] = i, ind_B_m[0]
            j += 1

    ind_NN = ind_NN[:j]
    
    toc = time.time()
    return ind_NN, toc-tic
    
    
def _nearest_neighbor_2nd_deprecated(m, n, ind_A, ind_B):
    """2nd Nearest neighbor finder for honeycomb lattice."""
    tic = time.time()
    
    ind_NN_2nd = np.zeros([m.shape[0]*6,2],dtype=np.int64) 
    
    j = 0
    for ind_sub in [ind_A, ind_B]:
        for i in ind_sub:
            ind_m = ind_sub[m[ind_sub] == m[i]]
            ind_n = ind_m[n[ind_m] == n[i] + 1]
            if ind_n:
                ind_NN_2nd[j,0], ind_NN_2nd[j,1] = i, ind_n[0]
                j += 1
            
            ind_n = ind_m[n[ind_m] == n[i] - 1]
            if ind_n:
                ind_NN_2nd[j,0], ind_NN_2nd[j,1] = i, ind_n[0]
                j += 1
                
            ind_n = ind_sub[n[ind_sub] == n[i]]
            ind_m = ind_n[m[ind_n] == m[i] + 1]
            if ind_m:
                ind_NN_2nd[j,0], ind_NN_2nd[j,1] = i, ind_m[0]
                j += 1
                
            ind_m = ind_n[m[ind_n] == m[i] - 1]
            if ind_m:
                ind_NN_2nd[j,0], ind_NN_2nd[j,1] = i, ind_m[0]
                j += 1
        
            ind_m = ind_sub[m[ind_sub] == m[i] + 1]
            ind_n = ind_m[n[ind_m] == n[i] - 1]
            if ind_n:
                ind_NN_2nd[j,0], ind_NN_2nd[j,1] = i, ind_n[0]
                j += 1
                
            ind_m = ind_sub[m[ind_sub] == m[i] - 1]
            ind_n = ind_m[n[ind_m] == n[i] + 1]
            if ind_n:
                ind_NN_2nd[j,0], ind_NN_2nd[j,1] = i, ind_n[0]
                j += 1
    
    ind_NN_2nd = ind_NN_2nd[:j]
    
    # Drop the lower triangular part
    ind_NN_2nd = ind_NN_2nd[ind_NN_2nd[:,0] < ind_NN_2nd[:,1],:]
    
    toc = time.time()
    return ind_NN_2nd, toc-tic
    

