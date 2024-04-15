import numpy as np
import nibabel as nib
import os
from scipy.stats import t

def compute_weighted_t_statistics(mean_effect_path, n_image_path, gaussian_folder, side, n_permutations=100):
    mean_effect_img = nib.load(mean_effect_path)
    mean_effect_data = mean_effect_img.get_fdata()

    n_img = nib.load(n_image_path)
    n_data = n_img.get_fdata()

    SSE_weighted = np.zeros_like(mean_effect_data)
    variance_weighted = np.zeros_like(mean_effect_data)
    t_statistic = np.zeros_like(mean_effect_data)
    total_weights = np.zeros_like(mean_effect_data)
    p_values = np.zeros_like(mean_effect_data)

    # Collect all patient data for permutation
    patient_data = []
    for filename in os.listdir(gaussian_folder):
        patient_id = filename.split('_')[0]
        if patient_id in clinical_improvement:
            gaussian_path = os.path.join(gaussian_folder, filename)
            gaussian_img = nib.load(gaussian_path)
            gaussian_data = gaussian_img.get_fdata()
            patient_data.append((gaussian_data, clinical_improvement[patient_id]))

    # Compute real t-statistic
    for gaussian_data, gk in patient_data:
        deviation = gk - mean_effect_data
        SSE_weighted += gaussian_data * (deviation ** 2)
        total_weights += gaussian_data

    # Compute t-statistics and prepare for permutations
    real_t_statistics = np.zeros_like(mean_effect_data)
    for i in np.ndindex(mean_effect_data.shape):
        if n_data[i] >= 6 and total_weights[i] > 0:
            variance_weighted[i] = SSE_weighted[i] / total_weights[i]
            SD_weighted = np.sqrt(variance_weighted[i])
            real_t_statistics[i] = mean_effect_data[i] / (SD_weighted / np.sqrt(total_weights[i]))

    # Permutation testing
    max_t_stats = []
    for _ in range(n_permutations):
        np.random.shuffle(patient_data)  # Shuffle patient data for permutation
        SSE_weighted_perm = np.zeros_like(mean_effect_data)

        for gaussian_data, gk in patient_data:
            deviation_perm = gk - mean_effect_data  # Use permuted labels
            SSE_weighted_perm += gaussian_data * (deviation_perm ** 2)

        max_t_stat = 0
        for i in np.ndindex(mean_effect_data.shape):
            if n_data[i] >= 6 and total_weights[i] > 0:
                variance_weighted_perm = SSE_weighted_perm[i] / total_weights[i]
                SD_weighted_perm = np.sqrt(variance_weighted_perm)
                t_stat_perm = mean_effect_data[i] / (SD_weighted_perm / np.sqrt(total_weights[i]))
                max_t_stat = max(max_t_stat, abs(t_stat_perm))
        max_t_stats.append(max_t_stat)

    # Calculate p-values based on permutation distribution
    for i in np.ndindex(mean_effect_data.shape):
        p_values[i] = np.sum(np.array(max_t_stats) >= abs(real_t_statistics[i])) / n_permutations

    # Save the t-statistic and p-value images
    t_statistic_img = nib.Nifti1Image(real_t_statistics, mean_effect_img.affine)
    t_statistic_path = os.path.join("/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/Multiple_comparisons", f"t_statistic2_{side}.nii")
    nib.save(t_statistic_img, t_statistic_path)

    p_value_img = nib.Nifti1Image(p_values, mean_effect_img.affine)
    p_value_path = os.path.join("/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/Multiple_comparisons", f"p_values_{side}.nii")
    nib.save(p_value_img, p_value_path)

    return t_statistic_path, p_value_path

# Run for both LH and RH

clinical_improvement = {
    'DERc': 60, 'SILf': 80, 'COCj': 94, 'HEMg': 14, 'SAMm': 95,
    'PASj': 97, 'FRAr': 96, 'BOUl': 80, 'LEBm': 11, 'RAFc': 31,
    'GATd': 95, 'JAKv': 0, 'VALa': 88, 'DAUd': 81, 'LAMp': 26,
    'FORa': 76, 'ILAf': 0, 'ERIi': 50, 'BOId': 33
}
mean_effect_img_path_lh = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/t_map/mean_effect_image_LH_resampled.nii"
mean_effect_img_path_rh = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/t_map/mean_effect_image_RH_resampled.nii"
n_img_path_lh = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/t_map/N_image_LH_resampled.nii"
n_img_path_rh = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/t_map/N_image_RH_resampled.nii"
gaussian_folders = {
    'LH': "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/create_gaussian_from_FATs/gaussianized_resampled_LH",
    'RH': "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/create_gaussian_from_FATs/gaussianized_resampled_RH"
}
t_statistic_path_lh, p_value_path_lh = compute_weighted_t_statistics(mean_effect_img_path_lh, n_img_path_lh, gaussian_folders['LH'], "LH")
t_statistic_path_rh, p_value_path_rh = compute_weighted_t_statistics(mean_effect_img_path_rh, n_img_path_rh, gaussian_folders['RH'], "RH")

print(f"Weighted t-statistic and p-value images saved for LH: {t_statistic_path_lh}, {p_value_path_lh}")
print(f"Weighted t-statistic and p-value images saved for RH: {t_statistic_path_rh}, {p_value_path_rh}")
