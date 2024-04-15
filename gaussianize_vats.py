# le but de ce code est de partir des volumes, de tissu activé binaire et de les pondéré par la distance au centre du volume sous la forme d'une Gaussienne 3D avec une largeur à mi-hauteur de 3,0 mm
import nibabel as nib
import numpy as np
from scipy.ndimage import gaussian_filter, center_of_mass



def gaussian_weighting(volume_path, fwhm=3.0):
    # Charger le volume binaire
    img = nib.load(volume_path)
    data = img.get_fdata()
    
    # Obtenir la résolution des voxels à partir de l'entête de l'image
    voxel_dims = img.header.get_zooms()  # Taille des voxels en mm (dx, dy, dz)
    
    # Convertir FWHM en sigma pour la gaussienne
    sigma_mm = fwhm / (2 * np.sqrt(2 * np.log(2)))  # Sigma en mm
    sigma_voxels = sigma_mm / np.mean(voxel_dims)  # Convertir sigma en unités de voxel
    
    # Calculer le centre de masse du volume activé
    center = center_of_mass(data)
    
    # Créer une grille de coordonnées 3D
    x = np.arange(data.shape[0])
    y = np.arange(data.shape[1])
    z = np.arange(data.shape[2])
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    # Calculer la distance de chaque voxel au centre en voxels
    distances_voxels = np.sqrt(((X - center[0]) ** 2 + (Y - center[1]) ** 2 + (Z - center[2]) ** 2))
    
    # Appliquer la pondération gaussienne en unités de voxel
    weighted_data = np.exp(-(distances_voxels ** 2) / (2 * sigma_voxels ** 2))
    
    # Appliquer la pondération uniquement aux voxels activés
    weighted_data *= data
    
    # Sauvegarder le résultat pondéré
    output_dir = "/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/create_gaussian_from_FATs/gaussianized_vats"
    output_filename = volume_path.split("/")[-1].replace('.nii', '_weighted.nii')
    output_path = f"{output_dir}/{output_filename}"
    weighted_img = nib.Nifti1Image(weighted_data, img.affine, img.header)
    nib.save(weighted_img, output_path)


# Liste des chemins des volumes à traiter
volume_paths = [
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/BINc_LH_VTA-segmentation-BINc_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/BINc_RH_VTA-segmentation-BINc_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/BOId_LH_VTA-segmentation-BOId_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/BOId_RH_VTA-segmentation-BOId_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/BOUl_LH_VTA-segmentation-BOUl_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/BOUl_RH_VTA-segmentation-BOUl_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/BRUl_LH_VTA-segmentation-BRUl_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/BRUl_RH_VTA-segmentation-BRUl_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/COCj_LH_VTA-segmentation-COCj_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/COCj_RH_VTA-segmentation-COCj_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/DAUd_LH_VTA-segmentation-DAUd_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/DAUd_RH_VTA-segmentation-DAUd_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/DERc_LH_VTA-segmentation-DERc_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/DERc_RH_VTA-segmentation-DERc_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/DUCm_LH_VTA-segmentation-DUCm_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/DUCm_RH_VTA-segmentation-DUCm_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/ERIi_LH_VTA-segmentation-ERIi_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/ERIi_RH_VTA-segmentation-ERIi_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/FORa_LH_VTA-segmentation-FORa_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/FORa_RH_VTA-segmentation-FORa_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/FRAr_LH_VTA-segmentation-FRAr_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/FRAr_RH_VTA-segmentation-FRAr_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/GATd_LH_VTA-segmentation-GATd_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/GATd_RH_VTA-segmentation-GATd_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/GONd_LH_VTA-segmentation-GONd_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/GONd_RH_VTA-segmentation-GONd_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/HEMg_LH_VTA-segmentation-HEMg_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/HEMg_RH_VTA-segmentation-HEMg_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/ILAf_LH_VTA-segmentation-ILAf_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/ILAf_RH_VTA-segmentation-ILAf_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/JAKv_LH_VTA-segmentation-JAKv_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/JAKv_RH_VTA-segmentation-JAKv_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/LAMp_LH_VTA-segmentation-LAMp_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/LAMp_RH_VTA-segmentation-LAMp_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/LEBm_LH_VTA-segmentation-LEBm_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/LEBm_RH_VTA-segmentation-LEBm_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/PASj_LH_VTA-segmentation-PASj_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/PASj_RH_VTA-segmentation-PASj_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/RAFc_LH_VTA-segmentation-RAFc_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/RAFc_RH_VTA-segmentation-RAFc_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/SAMm_LH_VTA-segmentation-SAMm_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/SAMm_RH_VTA-segmentation-SAMm_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/SILf_LH_VTA-segmentation-SILf_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/SILf_RH_VTA-segmentation-SILf_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/TAFc_LH_VTA-segmentation-TAFc_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/TAFc_RH_VTA-segmentation-TAFc_RH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/VALa_LH_VTA-segmentation-VALa_LH_VTA-label.nii",
    "/Users/nicolas.tempier/Desktop/These/Gaetan/VATs_in_YEB/VALa_RH_VTA-segmentation-VALa_RH_VTA-label.nii"
]

# Appliquer la pondération gaussienne à chaque volume
for path in volume_paths:
    gaussian_weighting(path)
