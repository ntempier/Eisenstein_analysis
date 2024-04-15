import numpy as np
import nibabel as nib
import os
from scipy.stats import t
from statsmodels.stats.multitest import multipletests
import random

def compute_weighted_t_statistics(mean_effect_path, n_image_path, gaussian_folder, side, n_permutations=100):
    mean_effect_img = nib.load(mean_effect_path)
    mean_effect_data = mean_effect_img.get_fdata()

    n_img = nib.load(n_image_path)
    n_data = n_img.get_fdata()

    SSE_weighted = np.zeros_like(mean_effect_data)
    variance_weighted = np.zeros_like(mean_effect_data)
    total_weights = np.zeros_like(mean_effect_data)

    # Collect all patient data for permutation
    patient_data = []
    for filename in os.listdir(gaussian_folder):
        patient_id = filename.split('_')[0]
        if patient_id in clinical_improvement:
            gaussian_path = os.path.join(gaussian_folder, filename)
            gaussian_img = nib.load(gaussian_path)
            gaussian_data = gaussian_img.get_fdata()
            patient_data.append((gaussian_data, clinical_improvement[patient_id]))

    # Compute real t-statistics and p-values
    real_t_statistics = np.zeros_like(mean_effect_data)
    for gaussian_data, gk in patient_data:
        deviation = gk - mean_effect_data
        SSE_weighted += gaussian_data * (deviation ** 2)
        total_weights += gaussian_data

    for i in np.ndindex(mean_effect_data.shape):
        if n_data[i] >= 6 and total_weights[i] > 0:
            variance_weighted[i] = SSE_weighted[i] / total_weights[i]
            SD_weighted = np.sqrt(variance_weighted[i])
            real_t_statistics[i] = mean_effect_data[i] / (SD_weighted / np.sqrt(total_weights[i]))

    real_p_values = 2 * t.sf(np.abs(real_t_statistics), n_data - 1)  # two-tailed p-value

    # Permutation testing for Q
    Q_observed = np.sum(-np.log10(real_p_values[real_p_values <= 0.05]))
    Q_permutations = []

    for _ in range(n_permutations):
        random.shuffle(patient_data)  # Permute clinical improvements
        SSE_weighted_perm = np.zeros_like(mean_effect_data)

        for gaussian_data, gk in patient_data:
            deviation_perm = gk - mean_effect_data
            SSE_weighted_perm += gaussian_data * (deviation_perm ** 2)

        t_statistics_perm = np.zeros_like(mean_effect_data)
        for i in np.ndindex(mean_effect_data.shape):
            if n_data[i] >= 6 and total_weights[i] > 0:
                variance_weighted_perm = SSE_weighted_perm[i] / total_weights[i]
                SD_weighted_perm = np.sqrt(variance_weighted_perm)
                t_statistics_perm[i] = mean_effect_data[i] / (SD_weighted_perm / np.sqrt(total_weights[i]))
        p_values_perm = 2 * t.sf(np.abs(t_statistics_perm), n_data - 1)
        Q_perm = np.sum(-np.log10(p_values_perm[p_values_perm <= 0.05]))
        Q_permutations.append(Q_perm)

    p_value = np.sum(np.array(Q_permutations) >= Q_observed) / n_permutations

    return p_value



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
