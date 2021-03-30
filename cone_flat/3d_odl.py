"""
Example using a filtered back-projection (FBP) in cone-beam 3d using `fbp_op`.

Note that the FBP is only approximate in this geometry, but still gives a
decent reconstruction that can be used as an initial guess in more complicated
methods.

Here we look at a partial scan, where the angular interval is not 2 * pi.
This caues issues for the regular FBP reconstruction, but can be improved
via a Parker weighting.

Note that since this is a fully 3d example, it may take some time to run,
about ~20s.
"""

import numpy as np
import odl
import h5py
import matplotlib.pyplot as plt


# --- Set up geometry of the problem --- #


# Reconstruction space: discretized functions on the cube
# [-20, 20]^3 with 300 samples per dimension.
reco_space = odl.uniform_discr(
    min_pt=[-255.5, -255.5, -12], max_pt=[255.5, 255.5, 12], shape=[512, 512, 25],
    dtype='float64')

# Make a circular cone beam geometry with flat detector and with a short scan
# Angles: uniformly spaced, n = 360, min = 0, max = 1.3 * pi
angle_partition = odl.uniform_partition(0, 274*np.pi/180, 275)
# angle_partition = odl.uniform_partition(0, 359*np.pi/180, 360)
# Detector: uniformly sampled, n = (512, 512), min = (-60, -60), max = (60, 60)
#detector_partition = odl.uniform_partition([-25.5, -493.50], [25.5, 493.50], [52, 988])
detector_partition = odl.uniform_partition([-906, -76], [906, 76], [1813, 153])
# Geometry with large cone and fan angle and tilted axis.
geometry = odl.tomo.ConeBeamGeometry(
    angle_partition, detector_partition, src_radius=500, det_radius=363)


# --- Create Filtered Back-projection (FBP) operator --- #


# Ray transform (= forward projection). We use the ASTRA CUDA backend.
ray_trafo = odl.tomo.RayTransform(reco_space, geometry, impl='astra_cuda')

# Create FBP operator using utility function
# We select a Shepp-Logan filter, and only use the lowest 80% of frequencies to
# avoid high frequency noise.
fbp = odl.tomo.fbp_op(ray_trafo,filter_type='Shepp-Logan', frequency_scaling=0.8)

# Apply parker weighting in order to improve reconstruction
parker_weighting = odl.tomo.parker_weighting(ray_trafo)
parker_weighted_fbp = fbp * parker_weighting


# --- Show some examples --- #


# Create a discrete Shepp-Logan phantom (modified version)

# phantom = odl.phantom.shepp_logan(reco_space, modified=True)

data = h5py.File('img_512x512x25.mat', 'r')
phantom = data['img'][:]
del data
phantom= np.transpose(phantom, [2, 1, 0])
phantom = phantom.astype('float64')
phantom = phantom[:,:,0:25]


# Create projection data by calling the ray transform on the phantom
proj_data = ray_trafo(phantom)

# Calculate filtered back-projection of data
fbp_reconstruction = fbp(proj_data)
pw_fbp_reconstruction = parker_weighted_fbp(proj_data)

# Shows a slice of the phantom, projections, and reconstruction
#phantom.show(title='Phantom')
#proj_data.show(title='Simulated Data (Sinogram)')
#fbp_reconstruction.show(title='Filtered Back-projection')
#pw_fbp_reconstruction.show(title='Parker-weighted Filtered Back-projection',
#                          force_show=True)
plt.imshow(phantom[:,:,4],'gray')
plt.show()
plt.imshow(proj_data[:,:,4],'gray')
plt.show()
plt.imshow(pw_fbp_reconstruction[:,:,4],'gray')
plt.show()
print('debug')
np.save('phantom.npy',phantom)
np.save('sin.npy',proj_data)
np.save('ct.npy',pw_fbp_reconstruction)
np.save('parker_weighting',parker_weighting)