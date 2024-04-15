import nibabel as nib
import numpy as np

# Load the N images
n_img_lh_path = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/t_map/N_image_LH_resampled.nii"
n_img_rh_path = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/t_map/N_image_RH_resampled.nii"
n_img_lh = nib.load(n_img_lh_path).get_fdata()
n_img_rh = nib.load(n_img_rh_path).get_fdata()

# Check if there are any voxels with N_i >= 6
print("Max N LH:", np.max(n_img_lh))
print("Max N RH:", np.max(n_img_rh))
print("Count of N_i >= 6 LH:", np.sum(n_img_lh >= 6))
print("Count of N_i >= 6 RH:", np.sum(n_img_rh >= 6))
