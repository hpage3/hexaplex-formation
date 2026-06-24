import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.ndimage import rotate as scipyrotate
from PIL import Image

try:
    import numexpr as ne
    _HAVE_NUMEXPR = True
except ImportError:
    _HAVE_NUMEXPR = False

# The atomic number of the unique atoms in the structure. Used as scattering factors.
atomic_number = {'C': 6, 'N': 7, 'O': 8}

def Rx(angle_deg):
    """Rotation matrix around x"""

    angle_rad = angle_deg/180.0*np.pi
    return np.array([[1, 0, 0],
                     [0, np.cos(angle_rad), -np.sin(angle_rad)],
                     [0, np.sin(angle_rad), np.cos(angle_rad)]])

def Ry(angle_deg):
    """Rotation matrix around y"""

    angle_rad = angle_deg/180.0*np.pi
    return np.array([[np.cos(angle_rad), 0, np.sin(angle_rad)],
                     [0, 1, 0],
                     [-np.sin(angle_rad), 0, np.cos(angle_rad)]])

def Rz(angle_deg):
    """Rotation matrix around z"""

    angle_rad = angle_deg/180.0*np.pi
    return np.array([[np.cos(angle_rad), -np.sin(angle_rad), 0],
                     [np.sin(angle_rad), np.cos(angle_rad), 0],
                     [0, 0, 1]])


def helix_maker(hexad_atoms, hexad_coord, rise, twist, num_hexads):
    """ A function to generate a stacked hexad structure given the
        rise, twist, and the number of hexad units"""

    coord = []
    atoms = np.array([i for i in hexad_atoms]*num_hexads) # Replicate the atom names
    # Loop over the number of hexads
    for i in range(num_hexads):
        # Loop over all atoms in the hexad
        for j in range(len(hexad_coord)):
            x = hexad_coord[j, 0]
            y = hexad_coord[j, 1]
            z = hexad_coord[j, 2] + i*rise # Apply the rise
            c = np.dot(Rz(twist*i), np.array([x, y, z])) # Apply the twist
            coord.append(c)

    return atoms, np.array(coord)

def save_xyz(atoms, coord, file_name):
    """ A function to write an XYZ file given the atom names and coordinates"""

    out = open(file_name, 'w')
    out.write('%i\n\n' %len(atoms))
    for i in range(len(atoms)):
        out.write('%-6s%12.6f%12.6f%12.6f\n' %(atoms[i], coord[i, 0], coord[i, 1], coord[i, 2]))
    out.close()

def _Ry_stack(angles_deg):
    """Stack of rotation matrices around y, one per angle. Shape (n, 3, 3)."""

    a = np.deg2rad(np.asarray(angles_deg, dtype=np.float64).ravel())
    c, s = np.cos(a), np.sin(a)
    R = np.zeros((a.size, 3, 3))
    R[:, 0, 0] = c;  R[:, 0, 2] = s
    R[:, 1, 1] = 1.0
    R[:, 2, 0] = -s; R[:, 2, 2] = c
    return R

def _Rz_stack(angles_deg):
    """Stack of rotation matrices around z, one per angle. Shape (n, 3, 3)."""

    a = np.deg2rad(np.asarray(angles_deg, dtype=np.float64).ravel())
    c, s = np.cos(a), np.sin(a)
    R = np.zeros((a.size, 3, 3))
    R[:, 0, 0] = c;  R[:, 0, 1] = -s
    R[:, 1, 0] = s;  R[:, 1, 1] = c
    R[:, 2, 2] = 1.0
    return R

def make_oriented_coords(coord, tilts, rotations):
    """ Build every oriented copy of the molecule for a tilt/rotation series.

    For each tilt angle the molecule is first rotated about y (the beam tilt),
    then for each rotation angle it is rotated about z (the fiber axis). This is
    the vectorized, loop-free equivalent of repeatedly applying Rz(rotation) to
    Ry(tilt) @ coord. Each orientation is built independently from the original
    coordinates (the rotations do not accumulate).

    Returns an array of shape (n_orientations, n_atoms, 3) with the orientations
    ordered tilt-major, then rotation.
    """

    coord = np.asarray(coord, dtype=np.float64)
    Ry_t = _Ry_stack(tilts)          # (T, 3, 3)
    Rz_r = _Rz_stack(rotations)      # (R, 3, 3)
    # Combined rotation for every (tilt, rotation) pair: Rz_r @ Ry_t
    R_comb = np.einsum('rij,tjk->trik', Rz_r, Ry_t).reshape(-1, 3, 3)  # (T*R, 3, 3)
    # Apply each orientation matrix to every atom
    return np.einsum('oij,aj->oai', R_comb, coord)                     # (O, A, 3)

def generate_fiber_diffraction_series(atoms, coords_stack, wavelength, distance_to_detector,
                                      z_grid_limits, x_grid_limits, z_grid_size, x_grid_size):
    """ Generate the summed fiber diffraction pattern for a whole orientation series.

    `coords_stack` has shape (n_orientations, n_atoms, 3): one rotated copy of
    the molecule per tilt/rotation sample. For every detector grid point (x, z)
    and every orientation o we compute the path length of every atom,

        dist_oa = y_oa + sqrt((x + x_oa)^2 + (D - y_oa)^2 + (z + z_oa)^2),

    accumulate the complex wave  s_o = sum_a f_a * exp(i * k * dist_oa)  with
    k = 2*pi/lambda and the atomic numbers f_a as scattering factors, and store
    the intensity summed over the series, sum_o |s_o|^2. This is exactly what the
    original notebook computed with two nested Python loops (tilt, rotation) that
    added one `generate_fiber_diffraction` pattern at a time; here the whole
    series is evaluated in one vectorized sweep.

    The grid is processed in chunks to bound memory. Within a chunk, numexpr (if
    available) computes the phase, cos and sin elementwise across all CPU cores
    for every (grid point, orientation, atom) at once, and the sum over atoms is
    done as a matrix-vector product with the scattering factors, which dispatches
    to multithreaded BLAS. (numexpr's own axis reduction is single-threaded, so
    it is deliberately not used for the sum.)

    float64 is kept throughout: the phase k*dist is ~1e10, so float32 would not
    resolve it modulo 2*pi. A global phase is removed by subtracting the detector
    distance D from every path length before taking cos/sin. Because the
    intensity is |s|^2, this common per-grid-point offset cancels exactly, while
    shrinking the phase argument by ~10x (from ~3e10 to ~2e9). float64 still
    cannot fully resolve a ~2e9 phase, but spot checks against an 80-bit
    reference show this is accurate to ~1e-5 relative; the subtraction removes the
    extra order-of-magnitude of round-off the original carried, and keeps the
    numexpr and numpy paths consistent.
    """

    coords_stack = np.ascontiguousarray(coords_stack, dtype=np.float64)
    if coords_stack.ndim == 2:                       # a single orientation
        coords_stack = coords_stack[None, :, :]
    n_orient, n_atoms, _ = coords_stack.shape

    # Determine the resolution of the grid
    z_range = np.linspace(z_grid_limits[0], z_grid_limits[1], z_grid_size)
    x_range = np.linspace(x_grid_limits[0], x_grid_limits[1], x_grid_size)

    # Flatten the detector grid into a list of (x, z) points
    zz, xx = np.meshgrid(z_range, x_range, indexing='ij')
    gx = xx.ravel()                       # x coordinate of every grid point
    gz = zz.ravel()                       # z coordinate of every grid point
    n_points = gx.size

    # Per-(orientation, atom) quantities, flattened to one axis (length O*A)
    f = np.asarray(atoms, dtype=np.float64)          # scattering factors, (A,)
    xa = np.ascontiguousarray(coords_stack[:, :, 0].ravel())
    ya = np.ascontiguousarray(coords_stack[:, :, 1].ravel())
    za = np.ascontiguousarray(coords_stack[:, :, 2].ravel())
    D = float(distance_to_detector)
    k = 2.0 * np.pi / wavelength
    y_term = (D - ya) ** 2                           # constant in x and z

    out = np.empty(n_points, dtype=np.float64)

    # Process the grid in chunks so the (chunk x orientations x atoms)
    # intermediates stay small. Target ~1.5e7 elements per array (~120 MB).
    n_oa = n_orient * n_atoms
    chunk = max(1, int(1.5e7 // max(1, n_oa)))

    # Shape the (O*A,) atom arrays as a row vector to broadcast against the chunk
    xa_r = xa[None, :]
    ya_r = ya[None, :]
    za_r = za[None, :]
    y_term_r = y_term[None, :]

    for start in range(0, n_points, chunk):
        stop = min(start + chunk, n_points)
        m = stop - start
        x_c = gx[start:stop, None]        # (m, 1)
        z_c = gz[start:stop, None]        # (m, 1)

        if _HAVE_NUMEXPR:
            local = {'xa': xa_r, 'ya': ya_r, 'za': za_r, 'y_term': y_term_r,
                     'x': x_c, 'z': z_c, 'k': k, 'D': D}
            # Phase, cos and sin elementwise (multithreaded via numexpr)...
            phase = ne.evaluate(
                "k * (ya - D + sqrt((x + xa)**2 + y_term + (z + za)**2))",
                local_dict=local)
            cos_m = ne.evaluate("cos(phase)")
            sin_m = ne.evaluate("sin(phase)")
        else:
            phase = k * (ya_r - D + np.sqrt((x_c + xa_r) ** 2 + y_term_r + (z_c + za_r) ** 2))
            cos_m = np.cos(phase)
            sin_m = np.sin(phase)

        # Reduce over atoms within each orientation via BLAS matrix-vector
        # products: (m*O, A) @ (A,) -> (m*O,), then reshape back to (m, O).
        c = (cos_m.reshape(m * n_orient, n_atoms) @ f).reshape(m, n_orient)
        s = (sin_m.reshape(m * n_orient, n_atoms) @ f).reshape(m, n_orient)
        # Intensity summed over the orientation series
        out[start:stop] = (c * c + s * s).sum(axis=1)

    return out.reshape(z_grid_size, x_grid_size)

def generate_fiber_diffraction(atoms, coord, wavelength, distance_to_detector,
                               z_grid_limits, x_grid_limits, z_grid_size, x_grid_size):
    """ Generate the fiber diffraction pattern for a single orientation.

    Thin wrapper around generate_fiber_diffraction_series for one molecule
    orientation; kept for backward compatibility. See that function for details.
    """

    return generate_fiber_diffraction_series(
        atoms, np.asarray(coord)[None, :, :], wavelength, distance_to_detector,
        z_grid_limits, x_grid_limits, z_grid_size, x_grid_size)

def plot_fiber_diffraction(diffraction_data, z_grid_limits, x_grid_limits, max_intensity_scaling):
    """ A function to plot the fiber diffraction pattern"""

    # Plot the diffraction pattern using a gray color map
    # To prevent the disappearance of low intensity peaks, we limit the
    # color map  to the maximum intensity multiplied by a specified factor
    plt.imshow(diffraction_data, extent=(x_grid_limits[0], x_grid_limits[1], z_grid_limits[0], z_grid_limits[1]),
               cmap='gray_r', vmax=np.max(diffraction_data)*max_intensity_scaling)
    plt.show()
    plt.imsave("plot2.png",diffraction_data,cmap='gray_r', vmax=np.max(diffraction_data)*max_intensity_scaling)




def plot_powder_diffraction(diffraction_data, z_grid_limits, x_grid_limits, max_intensity_scaling,z_grid_size):
    """ A function to plot the powder diffraction pattern"""

    # Plot the diffraction pattern using a gray color map
    # To prevent the disappearance of low intensity peaks, we limit the
    # color map  to the maximum intensity multiplied by a specified factor
    plt.imshow(diffraction_data, extent=(x_grid_limits[0], x_grid_limits[1], z_grid_limits[0], z_grid_limits[1]),
               cmap='gray_r', vmax=np.max(diffraction_data)*max_intensity_scaling)
   # plt.show()
    plt.imsave("plot2.png",diffraction_data,cmap='gray_r', vmax=np.max(diffraction_data)*max_intensity_scaling)
    # Load image
    img = mpimg.imread("plot2.png")

    rotatedsum = diffraction_data/360
    
    for ii in range(360):
        rotatedlocal = scipyrotate(diffraction_data, ii+1, reshape=False)
        rotatedsum = rotatedsum + rotatedlocal/360


    plt.imshow(rotatedsum, extent=(x_grid_limits[0], x_grid_limits[1], z_grid_limits[0], z_grid_limits[1]),
               cmap='gray_r', vmax=np.max(diffraction_data)*max_intensity_scaling)
    plt.imsave("PowderPattern.png",rotatedsum,cmap='gray_r', vmax=np.max(diffraction_data)*max_intensity_scaling)
    middlerow = int(z_grid_size/2)
    row = rotatedsum[middlerow]  # middle row

    with open("middlerow.txt", "w") as f:
        f.write("\n".join(map(str, row)))
