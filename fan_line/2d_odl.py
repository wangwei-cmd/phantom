"""
Example using a filtered back-projection (FBP) in fan beam using `fbp_op`.

Note that the FBP is only approximate in this geometry, but still gives a
decent reconstruction that can be used as an initial guess in more complicated
methods.

Here we look at a partial scan, where the angular interval is not 2 * pi.
This caues issues for the regular FBP reconstruction, but can be improved
via a Parker weighting.
"""

import numpy as np
import odl
import h5py
import matplotlib.pyplot as plt

# --- Set up geometry of the problem --- #


# Reconstruction space: discretized functions on the cube
# [-20, 20]^2 with 300 samples per dimension.
reco_space = odl.uniform_discr(
    min_pt=[-255.5, -255.5], max_pt=[255.5, 255.5], shape=[512, 512],
    dtype='float64')

# Make a circular cone beam geometry with flat detector
# Angles: uniformly spaced, n = 360, min = 0, max = pi + fan angle
angle_partition = odl.uniform_partition(0, 220*np.pi/180, 221)
# Detector: uniformly sampled, n = 512, min = -40, max = 40
detector_partition = odl.uniform_partition(-494, 494, 989)
# Geometry with large fan angle
geometry = odl.tomo.FanBeamGeometry(
    angle_partition, detector_partition, src_radius=1000, det_radius=363)


# --- Create Filtered Back-projection (FBP) operator --- #


# Ray transform (= forward projection). We use the ASTRA CUDA backend.
ray_trafo = odl.tomo.RayTransform(reco_space, geometry, impl='astra_cuda')

# Create FBP operator using utility function
# We select a Hann filter, and only use the lowest 80% of frequencies to avoid
# high frequency noise.
fbp = odl.tomo.fbp_op(ray_trafo, filter_type='Ram-Lak', frequency_scaling=1.0)

# Apply parker weighting in order to improve reconstruction
parker_weighting = odl.tomo.parker_weighting(ray_trafo)
parker_weighting.show()
parker_weighted_fbp = fbp * parker_weighting


# --- Show some examples --- #


# Create a discrete Shepp-Logan phantom (modified version)
# phantom = odl.phantom.shepp_logan(reco_space, modified=True)

data = h5py.File('f.mat', 'r')
phantom = data['f'][:]
del data
phantom= np.transpose(phantom, [1, 0])
phantom = phantom.astype('float64')



# Create projection data by calling the ray transform on the phantom
proj_data = ray_trafo(phantom)

# Calculate filtered back-projection of data
fbp_reconstruction = fbp(proj_data)
pw_fbp_reconstruction = parker_weighted_fbp(proj_data)

# Shows a slice of the phantom, projections, and reconstruction
# phantom.show(title='Phantom')
# proj_data.show(title='Projection Data (Sinogram)')
# fbp_reconstruction.show(title='Filtered Back-projection')
# pw_fbp_reconstruction.show(title='Parker-weighted Filtered Back-projection',
#                            force_show=True)

plt.imshow(phantom,'gray')
plt.show()
plt.imshow(proj_data,'gray')
plt.show()
plt.imshow(pw_fbp_reconstruction,'gray')
plt.show()

np.save('parker_weighting.npy',parker_weighting)
np.save('sin.npy',proj_data)
np.save('ct.npy',pw_fbp_reconstruction)