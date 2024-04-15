import nibabel as nib
import matplotlib.pyplot as plt

# Load the NIfTI file
nii_file = '/Users/nicolas.tempier/Desktop/These/Gaetan/Methode_eisenstein/FDR/p_values_corrected_RH.nii'
image = nib.load(nii_file)
data = image.get_fdata()

# Flatten the image data to 1D for histogram plotting
data_flattened = data.flatten()

# Plot the histogram
plt.hist(data_flattened, bins=50, color='blue', alpha=0.7)
plt.title('Histogram of Image Values')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.show()
