from scipy.stats import t
import numpy as np
import nibabel as nib
import os

def compute_p_values(t_statistic_img_path, n_img_path,p_value_img_path):
    # Load the t-statistic and N images
    t_statistic_img = nib.load(t_statistic_img_path)
    t_statistic_data = t_statistic_img.get_fdata()
    
    n_img = nib.load(n_img_path)
    n_data = n_img.get_fdata()
    
    # Initialize the p-value array with the same shape as the t-statistic image
    p_values = np.zeros_like(t_statistic_data)
    
    # Calculate p-values for each voxel
    for i in np.ndindex(t_statistic_data.shape):
        if n_data[i] > 1:  # Need at least 2 observations to have a degree of freedom
            # Degrees of freedom
            df = n_data[i] - 1
            # Convert t-statistic to p-value using a two-tailed test
            p_values[i] = 2 * t.sf(np.abs(t_statistic_data[i]), df)
        else:
            p_values[i] = 1  # Assign a p-value of 1 where we don't have enough data
    
    # Save the p-value image
    p_value_img = nib.Nifti1Image(p_values, t_statistic_img.affine)
    nib.save(p_value_img, p_value_img_path)
    
    return p_value_img_path

# Calculate and save p-value images for LH and RH
p_value_img_path_lh = compute_p_values("/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/t_map/t_statistic2_LH.nii", "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/t_map/N_image_LH_resampled.nii","/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/P_map/p_values2_LH.nii")
p_value_img_path_rh = compute_p_values("/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/t_map/t_statistic2_RH.nii", "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/t_map/N_image_RH_resampled.nii","/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/P_map/p_values2_RH.nii")

print(f"P-value image saved for LH: {p_value_img_path_lh}")
print(f"P-value image saved for RH: {p_value_img_path_rh}")
