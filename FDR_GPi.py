import numpy as np
import nibabel as nib
import os
from scipy.stats import t
from statsmodels.stats.multitest import multipletests

# Définition des chemins et autres éléments comme précédemment

def compute_weighted_t_statistics(mean_effect_path, n_image_path, gaussian_folder, side):
    mean_effect_img = nib.load(mean_effect_path)
    mean_effect_data = mean_effect_img.get_fdata()

    n_img = nib.load(n_image_path)
    n_data = n_img.get_fdata()

    SSE_weighted = np.zeros_like(mean_effect_data)
    variance_weighted = np.zeros_like(mean_effect_data)
    t_statistic = np.zeros_like(mean_effect_data)
    p_values = np.zeros_like(mean_effect_data)

    total_weights = np.zeros_like(mean_effect_data)

    for filename in os.listdir(gaussian_folder):
        patient_id = filename.split('_')[0]
        if patient_id in clinical_improvement:
            gk = clinical_improvement[patient_id]
            gaussian_path = os.path.join(gaussian_folder, filename)
            gaussian_img = nib.load(gaussian_path)
            gaussian_data = gaussian_img.get_fdata()

            deviation = gk - mean_effect_data
            SSE_weighted += gaussian_data * (deviation ** 2)
            total_weights += gaussian_data

    for i in np.ndindex(mean_effect_data.shape):
        if n_data[i] >= 6:
            if total_weights[i] > 0:
                variance_weighted[i] = SSE_weighted[i] / total_weights[i]
                if variance_weighted[i] > 0:  # Vérifiez que la variance est positive avant de calculer SD
                    SD_weighted = np.sqrt(variance_weighted[i])
                    t_statistic[i] = mean_effect_data[i] / (SD_weighted / np.sqrt(total_weights[i]))
                    df = n_data[i] - 1  # Degrés de liberté
                    p_values[i] = 2 * t.sf(np.abs(t_statistic[i]), df)  # p-value bilatérale
                else:
                    t_statistic[i] = 0
                    p_values[i] = 1
            else:
                t_statistic[i] = 0
                p_values[i] = 1
        else:
            t_statistic[i] = 0
            p_values[i] = 1

    # Correction FDR
    alpha_value = 0.1  # Modifiable selon vos besoins
    rejected, p_values_corrected, _, _ = multipletests(p_values.flatten(), alpha=alpha_value, method='fdr_bh')
    p_values_corrected = p_values_corrected.reshape(p_values.shape)

    # Sauvegarde des images
    t_statistic_img = nib.Nifti1Image(t_statistic, mean_effect_img.affine)
    p_value_corrected_img = nib.Nifti1Image(p_values_corrected, mean_effect_img.affine)
    # Chemins pour sauvegarder les résultats
    t_statistic_path = os.path.join("/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/FDR", f"t_statistic2_{side}_FDR.nii")
    p_value_path = os.path.join("/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/FDR", f"p_values_corrected_{side}_FDR.nii")

    nib.save(t_statistic_img, t_statistic_path)
    nib.save(p_value_corrected_img, p_value_path)

    return t_statistic_path, p_value_path

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
# Compute and save weighted t-statistics and corrected p-values for LH and RH
t_statistic_path_lh, p_value_path_lh = compute_weighted_t_statistics(mean_effect_img_path_lh, n_img_path_lh, gaussian_folders['LH'], "LH")
t_statistic_path_rh, p_value_path_rh = compute_weighted_t_statistics(mean_effect_img_path_rh, n_img_path_rh, gaussian_folders['RH'], "RH")

print(f"Weighted t-statistic and corrected p-value images saved for LH: {t_statistic_path_lh}, {p_value_path_lh}")
print(f"Weighted t-statistic and corrected p-value images saved for RH: {t_statistic_path_rh}, {p_value_path_rh}")
