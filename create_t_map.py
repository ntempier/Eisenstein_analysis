import numpy as np
import nibabel as nib
import os

# Paths
mean_effect_img_path_lh = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/t_map/mean_effect_image_LH_resampled.nii"
mean_effect_img_path_rh = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/t_map/mean_effect_image_RH_resampled.nii"
n_img_path_lh = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/t_map/N_image_LH_resampled.nii"
n_img_path_rh = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/t_map/N_image_RH_resampled.nii"

gaussian_folders = {
    'LH': "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/create_gaussian_from_FATs/gaussianized_resampled_LH",
    'RH': "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/create_gaussian_from_FATs/gaussianized_resampled_RH"
}

clinical_improvement = {
    'DERc': 60, 'SILf': 80, 'COCj': 94, 'HEMg': 14, 'SAMm': 95,
    'PASj': 97, 'FRAr': 96, 'BOUl': 80, 'LEBm': 11, 'RAFc': 31,
    'GATd': 95, 'JAKv': 0, 'VALa': 88, 'DAUd': 81, 'LAMp': 26,
    'FORa': 76, 'ILAf': 0, 'ERIi': 50, 'BOId': 33
}

def compute_weighted_t_statistics(mean_effect_path, n_image_path, gaussian_folder, side):
    mean_effect_img = nib.load(mean_effect_path)
    mean_effect_data = mean_effect_img.get_fdata()
    
    n_img = nib.load(n_image_path)
    n_data = n_img.get_fdata()

    SSE_weighted = np.zeros_like(mean_effect_data)
    variance_weighted = np.zeros_like(mean_effect_data)
    t_statistic = np.zeros_like(mean_effect_data)

    total_weights = np.zeros_like(mean_effect_data)

    for filename in os.listdir(gaussian_folder):
        patient_id = filename.split('_')[0]
        if patient_id in clinical_improvement:
            gk = clinical_improvement[patient_id]
            gaussian_path = os.path.join(gaussian_folder, filename)
            gaussian_img = nib.load(gaussian_path)
            gaussian_data = gaussian_img.get_fdata()
            
            deviation = gk - mean_effect_data  # Calcul de l'écart à la moyenne
            SSE_weighted += gaussian_data * (deviation ** 2)  # Mise à jour de la somme des carrés des écarts pondérés
            total_weights += gaussian_data  # Mise à jour des poids totaux pour chaque voxel

    for i in np.ndindex(mean_effect_data.shape):
        if n_data[i] >= 6:  # Appliquer un seuil sur le nombre de participants
            Ni = n_data[i]
            if total_weights[i] > 0:  # Éviter la division par zéro
                variance_weighted[i] = SSE_weighted[i] / total_weights[i]
                SD_weighted = np.sqrt(variance_weighted[i])
                t_statistic[i] = mean_effect_data[i] / (SD_weighted / np.sqrt(total_weights[i]))
            else:
                t_statistic[i] = 0
        else:
            t_statistic[i] = 0

    t_statistic_img = nib.Nifti1Image(t_statistic, mean_effect_img.affine)
    t_statistic_path = os.path.join("/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/t_map", f"t_statistic2_{side}.nii")
    nib.save(t_statistic_img, t_statistic_path)

    return t_statistic_path

# Compute and save weighted t-statistics for LH and RH
t_statistic_path_lh = compute_weighted_t_statistics(mean_effect_img_path_lh, n_img_path_lh, gaussian_folders['LH'],"LH")
t_statistic_path_rh = compute_weighted_t_statistics(mean_effect_img_path_rh, n_img_path_rh, gaussian_folders['RH'],"RH")

print(f"Weighted t-statistic image saved for LH: {t_statistic_path_lh}")
print(f"Weighted t-statistic image saved for RH: {t_statistic_path_rh}")
