import nibabel as nib
from nilearn.image import resample_to_img
import os

# Define your images and paths
image_paths = [
    "BINc_LH_VTA_BINc_gaussian.nii",
    "BINc_RH_VTA_BINc_gaussian.nii",
    "BOId_LH_VTA_BOId_gaussian.nii",
    "BOId_RH_VTA_BOId_gaussian.nii",
    "BOUl_LH_VTA_BOUl_gaussian.nii",
    "BOUl_RH_VTA_BOUl_gaussian.nii",
    "COCj_LH_VTA_COCj_gaussian.nii",
    "COCj_RH_VTA_COCj_gaussian.nii",
    "DAUd_LH_VTA_DAUd_gaussian.nii",
    "DAUd_RH_VTA_DAUd_gaussian.nii",
    "DERc_LH_VTA_DERc_gaussian.nii",
    "DERc_RH_VTA_DERc_gaussian.nii",
    "ERIi_LH_VTA_ERIi_gaussian.nii",
    "ERIi_RH_VTA_ERIi_gaussian.nii",
    "FORa_LH_VTA_FORa_gaussian.nii",
    "FORa_RH_VTA_FORa_gaussian.nii",
    "FRAr_LH_VTA_FRAr_gaussian.nii",
    "FRAr_RH_VTA_FRAr_gaussian.nii",
    "GATd_LH_VTA_GATd_gaussian.nii",
    "GATd_RH_VTA_GATd_gaussian.nii",
    "HEMg_LH_VTA_HEMg_gaussian.nii",
    "HEMg_RH_VTA_HEMg_gaussian.nii",
    "ILAf_LH_VTA_ILAf_gaussian.nii",
    "ILAf_RH_VTA_ILAf_gaussian.nii",
    "JAKv_LH_VTA_JAKv_gaussian.nii",
    "JAKv_RH_VTA_JAKv_gaussian.nii",
    "LAMp_LH_VTA_LAMp_gaussian.nii",
    "LAMp_RH_VTA_LAMp_gaussian.nii",
    "LEBm_LH_VTA_LEBm_gaussian.nii",
    "LEBm_RH_VTA_LEBm_gaussian.nii",
    "PASj_LH_VTA_PASj_gaussian.nii",
    "PASj_RH_VTA_PASj_gaussian.nii",
    "RAFc_LH_VTA_RAFc_gaussian.nii",
    "RAFc_RH_VTA_RAFc_gaussian.nii",
    "SAMm_LH_VTA_SAMm_gaussian.nii",
    "SAMm_RH_VTA_SAMm_gaussian.nii",
    "SILf_LH_VTA_SILf_gaussian.nii",
    "SILf_RH_VTA_SILf_gaussian.nii",
    "VALa_LH_VTA_VALa_gaussian.nii",
    "VALa_RH_VTA_VALa_gaussian.nii"
]

reference_image_path = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/74_t1mri_blank_resampled.nii"
output_dir = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/create_gaussian_from_FATs/gaussianized_resampled"

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Load the reference image
reference_img = nib.load(reference_image_path)

# Process each image
# for img_name in image_paths:
#     # Construct the full path for the image
#     img_path = os.path.join("/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/create_gaussian_from_FATs/gaussianized_vats", img_name)
    
#     # Load the image
#     img = nib.load(img_path)
    
#     # Resample the image to match the reference
#     resampled_img = resample_to_img(img, reference_img)
    
#     # Construct the output path
#     output_path = os.path.join(output_dir, img_name)
    
#     # Save the resampled image
#     nib.save(resampled_img, output_path)

#     print(f"Resampled and saved: {output_path}")


# resample /Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/Mean_effect_image/mean_effect_image_LH.nii
# /Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/Mean_effect_image/mean_effect_image_RH.nii
# /Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/Threshold_N_image_0_05/N_image_LH_thresholded.nii
# /Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/Threshold_N_image_0_05/N_image_RH_thresholded.nii

mean_effect_img_path_lh = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/Mean_effect_image/mean_effect_image_LH.nii"
mean_effect_img_path_rh = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/Mean_effect_image/mean_effect_image_RH.nii"
n_img_path_lh = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/Threshold_N_image_0_05/N_image_LH_thresholded.nii"
n_img_path_rh = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/Threshold_N_image_0_05/N_image_RH_thresholded.nii"

mean_effect_img_lh= nib.load(mean_effect_img_path_lh)
n_img_lh= nib.load(n_img_path_lh)
mean_effect_img_rh= nib.load(mean_effect_img_path_rh)
n_img_rh= nib.load(n_img_path_rh)

# Resample the images
resampled_mean_effect_img_lh = resample_to_img(mean_effect_img_lh, reference_img, interpolation='nearest')
resampled_n_img_lh = resample_to_img(n_img_lh, reference_img, interpolation='nearest')
resampled_mean_effect_img_rh = resample_to_img(mean_effect_img_rh, reference_img, interpolation='nearest')
resampled_n_img_rh = resample_to_img(n_img_rh, reference_img, interpolation='nearest')

output_dir2 = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/t_map"

# nib.save(resampled_mean_effect_img_lh, os.path.join(output_dir2, "mean_effect_image_LH_resampled.nii"))
nib.save(resampled_n_img_lh, os.path.join(output_dir2, "N_image_LH_resampled.nii"))
# nib.save(resampled_mean_effect_img_rh, os.path.join(output_dir2, "mean_effect_image_RH_resampled.nii"))
nib.save(resampled_n_img_rh, os.path.join(output_dir2, "N_image_RH_resampled.nii"))
