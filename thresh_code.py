import nibabel as nib
import numpy as np
import os 

# Chemins des images N pour chaque côté
path_lh = '/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/N_image/N_image_LH.nii'
path_rh = '/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/N_image/N_image_RH.nii'

def apply_threshold_to_n_image(path, threshold=0.05):
    # Charger l'image N
    n_img = nib.load(path)
    n_data = n_img.get_fdata()

    # Déterminer le seuil basé sur le poids maximal * 0.05
    max_weight = np.max(n_data) * threshold

    # Appliquer le seuil
    n_data[n_data <= max_weight] = 0

    # Sauvegarder l'image modifiée
    output_dir = '/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/Threshold_N_image_0_05'
    output_path = path.replace('.nii', '_thresholded.nii')
    output_path = os.path.join(output_dir, os.path.basename(output_path))
    thresholded_img = nib.Nifti1Image(n_data, n_img.affine, n_img.header)
    nib.save(thresholded_img, output_path)
    return output_path

# Appliquer le seuil pour chaque image
output_lh = apply_threshold_to_n_image(path_lh)
output_rh = apply_threshold_to_n_image(path_rh)

print("Les images traitées ont été sauvegardées à :")
print("Côté gauche :", output_lh)
print("Côté droit :", output_rh)
