import numpy as np
import matplotlib.pyplot as plt

try:
    import numexpr as ne

    _HAVE_NUMEXPR = True
except ImportError:
    _HAVE_NUMEXPR = False

# The atomic number of the unique atoms in the structure. Used as scattering factors.
atomic_number = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'P': 15}

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
    """Stack rotation matrices around y, one per angle. Shape (n, 3, 3)."""

    a = np.deg2rad(np.asarray(angles_deg, dtype=np.float64).ravel())
    c, s = np.cos(a), np.sin(a)
    R = np.zeros((a.size, 3, 3))
    R[:, 0, 0] = c
    R[:, 0, 2] = s
    R[:, 1, 1] = 1.0
    R[:, 2, 0] = -s
    R[:, 2, 2] = c
    return R

def _Rz_stack(angles_deg):
    """Stack rotation matrices around z, one per angle. Shape (n, 3, 3)."""

    a = np.deg2rad(np.asarray(angles_deg, dtype=np.float64).ravel())
    c, s = np.cos(a), np.sin(a)
    R = np.zeros((a.size, 3, 3))
    R[:, 0, 0] = c
    R[:, 0, 1] = -s
    R[:, 1, 0] = s
    R[:, 1, 1] = c
    R[:, 2, 2] = 1.0
    return R

def make_oriented_coords(coord, tilts, rotations):
    """Build independent oriented copies for every tilt/azimuth pair.

    Each azimuthal rotation is applied to the same tilted coordinate frame for a
    given tilt; rotations do not accumulate across the series.
    """

    coord = np.asarray(coord, dtype=np.float64)
    Ry_t = _Ry_stack(tilts)
    Rz_r = _Rz_stack(rotations)
    combined = np.einsum("rij,tjk->trik", Rz_r, Ry_t).reshape(-1, 3, 3)
    return np.einsum("oij,aj->oai", combined, coord)

def generate_fiber_diffraction_series(atoms, coords_stack, wavelength, distance_to_detector,
                                      z_grid_limits, x_grid_limits, z_grid_size, x_grid_size):
    """Generate summed diffraction for one or more independent orientations.

    `coords_stack` may be a single `(n_atoms, 3)` coordinate array or an
    `(n_orientations, n_atoms, 3)` stack. The return value is the summed
    detector image over orientations, matching repeated calls to
    `generate_fiber_diffraction` with each oriented coordinate frame.
    """

    coords_stack = np.ascontiguousarray(coords_stack, dtype=np.float64)
    if coords_stack.ndim == 2:
        coords_stack = coords_stack[None, :, :]
    n_orient, n_atoms, _ = coords_stack.shape

    z_range = np.linspace(z_grid_limits[0], z_grid_limits[1], z_grid_size)
    x_range = np.linspace(x_grid_limits[0], x_grid_limits[1], x_grid_size)
    zz, xx = np.meshgrid(z_range, x_range, indexing="ij")
    gx = xx.ravel()
    gz = zz.ravel()
    n_points = gx.size

    f = np.asarray(atoms, dtype=np.float64)
    xa = np.ascontiguousarray(coords_stack[:, :, 0].ravel())
    ya = np.ascontiguousarray(coords_stack[:, :, 1].ravel())
    za = np.ascontiguousarray(coords_stack[:, :, 2].ravel())
    D = float(distance_to_detector)
    k = 2.0 * np.pi / wavelength
    y_term = (D - ya) ** 2

    out = np.empty(n_points, dtype=np.float64)
    n_oa = n_orient * n_atoms
    chunk = max(1, int(1.5e7 // max(1, n_oa)))

    xa_r = xa[None, :]
    ya_r = ya[None, :]
    za_r = za[None, :]
    y_term_r = y_term[None, :]

    for start in range(0, n_points, chunk):
        stop = min(start + chunk, n_points)
        m = stop - start
        x_c = gx[start:stop, None]
        z_c = gz[start:stop, None]

        if _HAVE_NUMEXPR:
            local = {
                "xa": xa_r,
                "ya": ya_r,
                "za": za_r,
                "y_term": y_term_r,
                "x": x_c,
                "z": z_c,
                "k": k,
                "D": D,
            }
            phase = ne.evaluate(
                "k * (ya - D + sqrt((x + xa)**2 + y_term + (z + za)**2))",
                local_dict=local,
            )
            cos_m = ne.evaluate("cos(phase)")
            sin_m = ne.evaluate("sin(phase)")
        else:
            phase = k * (
                ya_r
                - D
                + np.sqrt((x_c + xa_r) ** 2 + y_term_r + (z_c + za_r) ** 2)
            )
            cos_m = np.cos(phase)
            sin_m = np.sin(phase)

        c = (cos_m.reshape(m * n_orient, n_atoms) @ f).reshape(m, n_orient)
        s = (sin_m.reshape(m * n_orient, n_atoms) @ f).reshape(m, n_orient)
        out[start:stop] = (c * c + s * s).sum(axis=1)

    return out.reshape(z_grid_size, x_grid_size)

def generate_fiber_diffraction(atoms, coord, wavelength, distance_to_detector,
                               z_grid_limits, x_grid_limits, z_grid_size, x_grid_size):
    """A function to generate the fiber diffraction pattern."""

    return generate_fiber_diffraction_series(
        atoms,
        np.asarray(coord)[None, :, :],
        wavelength,
        distance_to_detector,
        z_grid_limits,
        x_grid_limits,
        z_grid_size,
        x_grid_size,
    )

def plot_fiber_diffraction(diffraction_data, z_grid_limits, x_grid_limits, max_intensity_scaling):
    """ A function to plot the fiber diffraction pattern"""

    # Plot the diffraction pattern using a gray color map
    # To prevent the disappearance of low intensity peaks, we limit the
    # color map  to the maximum intensity multiplied by a specified factor
    plt.imshow(diffraction_data, extent=(x_grid_limits[0], x_grid_limits[1], z_grid_limits[0], z_grid_limits[1]),
               cmap='gray_r', vmax=np.max(diffraction_data)*max_intensity_scaling)
    plt.show()
