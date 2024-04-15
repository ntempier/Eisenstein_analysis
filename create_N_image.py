import os
import nibabel as nib
import numpy as np
from scipy.ndimage import affine_transform

# Fonction pour rééchantillonner les images des participants à la résolution de l'image de référence
def resample_to_ref(participant_img, ref_img):
    participant_data = participant_img.get_fdata()
    participant_affine = participant_img.affine
    ref_affine = ref_img.affine
    ref_shape = ref_img.shape

    # Calculer la transformation affine relative entre les espaces de l'image de référence et du participant
    relative_affine = np.linalg.inv(ref_affine) @ participant_affine

    # Rééchantillonner l'image du participant à l'espace de l'image de référence
    resampled_data = affine_transform(
        participant_data,
        matrix=relative_affine[:3, :3],
        offset=np.dot(relative_affine, np.array([0, 0, 0, 1]))[:3],
        output_shape=ref_shape,
        order=1,  # ordre d'interpolation linéaire
        mode='nearest',  # gestion des valeurs hors limites
        cval=0  # valeur pour les points hors de l'image d'origine
    )
    return resampled_data

# Chemin vers l'image de référence vierge
reference_path = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/74_t1mri_blank.nii.gz"

# Charger l'image de référence
ref_img = nib.load(reference_path)
ref_data = np.zeros(ref_img.shape, dtype=int)  # Initialiser avec des zéros

# Liste des chemins des volumes à traiter
participant_images = [
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/COCj_LH_VTA-segmentation-COCj_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/BOId_LH_VTA-segmentation-BOId_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/LAMp_LH_VTA-segmentation-LAMp_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/FRAr_LH_VTA-segmentation-FRAr_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/GATd_LH_VTA-segmentation-GATd_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/FORa_LH_VTA-segmentation-FORa_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/HEMg_LH_VTA-segmentation-HEMg_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/SAMm_LH_VTA-segmentation-SAMm_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/DERc_LH_VTA-segmentation-DERc_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/VALa_LH_VTA-segmentation-VALa_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/BINc_LH_VTA-segmentation-BINc_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/DAUd_LH_VTA-segmentation-DAUd_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/SILf_LH_VTA-segmentation-SILf_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/BOUl_LH_VTA-segmentation-BOUl_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/ERIi_LH_VTA-segmentation-ERIi_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/PASj_LH_VTA-segmentation-PASj_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/RAFc_LH_VTA-segmentation-RAFc_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/ILAf_LH_VTA-segmentation-ILAf_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/JAKv_LH_VTA-segmentation-JAKv_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/LEBm_LH_VTA-segmentation-LEBm_LH_VTA-label.nii",
]

# Parcourir chaque image de participant et comptabiliser les contributions
for img_path in participant_images:
    participant_img = nib.load(img_path)
    resampled_data = resample_to_ref(participant_img, ref_img)
    ref_data += (resampled_data > 0)  # Incrémenter les compteurs pour les voxels actifs

# Créer une nouvelle image Nifti avec les données comptabilisées
n_image = nib.Nifti1Image(ref_data, ref_img.affine, ref_img.header)

# Définir le répertoire de sortie
output_dir = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/N_image"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Sauvegarder l'image N dans le dossier spécifié
output_path = os.path.join(output_dir, "N_image.nii")
nib.save(n_image, output_path)

print("L'image N a été créée et sauvegardée à :", output_path)

