# lowfield_maxgirf

This repository contains the code and datasets for
**"MaxGIRF: Image Reconstruction Incorporating Concomitant Field and Gradient Impulse Response Function Effects"**, by Nam G. Lee, Rajiv Ramasawmy, Yongwan Lim, Adrienne E. Campbell-Washburn, and Krishna S. Nayak.

Nam Gyun Lee, University of Southern California, Dec 2021.

## ISMRMRD datasets
Example human and phantom datasets can be found on Zenodo:

- Human: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5830910.svg)](https://doi.org/10.5281/zenodo.5830910)
- Phantom: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5874604.svg)](https://doi.org/10.5281/zenodo.5874604)

* **Axial Cartesian spin-echo dataset**
  - meas_MID00273_FID03656_se_15b130_tra.h5
  - noise_meas_MID00273_FID03656_se_15b130_tra.h5

* **Axial spiral spin-echo dataset** 
  - meas_MID00275_FID03658_se_spiral_1102_ax_s24.h5
  - noise_meas_MID00275_FID03658_se_spiral_1102_ax_s24.h5

* **Axial Cartesian GRE datasets for static off-resonance map estimation**

  - meas_MID00276_FID03659_gre_TE1.h5
  - meas_MID00277_FID03660_gre_TE2.h5
  - meas_MID00278_FID03661_gre_TE3.h5
  - meas_MID00279_FID03662_gre_TE4.h5
  - meas_MID00280_FID03663_gre_TE5.h5
  - noise_meas_MID00276_FID03659_gre_TE1.h5
  - noise_meas_MID00276_FID03660_gre_TE2.h5
  - noise_meas_MID00276_FID03661_gre_TE3.h5
  - noise_meas_MID00276_FID03662_gre_TE4.h5
  - noise_meas_MID00276_FID03663_gre_TE5.h5

* **Sagittal Cartesian spin-echo dataset**
  - meas_MID00258_FID03641_se_15b130.h5
  - noise_meas_MID00258_FID03641_se_15b130.h5
  
* **Sagittal spiral spin-echo dataset**
  - meas_MID00260_FID03643_se_spiral_1102_sag_s24.h5
  - noise_meas_MID00260_FID03643_se_spiral_1102_sag_s24.h5
  
* **Sagittal Cartesian GRE datasets for static off-resonance map estimation**
  - meas_MID00261_FID03644_gre_TE1.h5
  - meas_MID00262_FID03645_gre_TE2.h5
  - meas_MID00263_FID03646_gre_TE3.h5
  - meas_MID00264_FID03647_gre_TE4.h5
  - meas_MID00265_FID03648_gre_TE5.h5
  - noise_meas_MID00261_FID03644_gre_TE1.h5
  - noise_meas_MID00261_FID03645_gre_TE2.h5
  - noise_meas_MID00261_FID03646_gre_TE3.h5
  - noise_meas_MID00261_FID03647_gre_TE4.h5
  - noise_meas_MID00261_FID03648_gre_TE5.h5

* **Mat files containing a static off-resonance map obtained with [NLINV estimation](https://github.com/usc-mrel/nlinv_estimation)**
  - B0map_nlinv_min1.0e-06_axial.mat
  - B0map_nlinv_min1.0e-06_sagittal_ro0.mat (ro := remove oversampling: 1=yes, 0=no)
  - B0map_nlinv_min1.0e-06_sagittal_ro1.mat

Once datasets are downloaded, organize datasets in this structure:
 
     |---path-to-dataset (e.g., D:\lowfield\NHLBI\data\20201102_NV_brain)
         |---meas_MID00260_FID03643_se_spiral_1102_sag_s24.h5
         |---meas_MID00275_FID03658_se_spiral_1102_ax_s24.h5
         |---noise_meas_MID00260_FID03643_se_spiral_1102_sag_s24.h5
         |---noise_meas_MID00275_FID03658_se_spiral_1102_ax_s24.h5
         |---B0map_nlinv_min1.0e-06_axial.mat
         |---B0map_nlinv_min1.0e-06_sagittal_ro0.mat
         |---B0map_nlinv_min1.0e-06_sagittal_ro1.mat

## Figure 3 (numerical simulation)

Update variables `src_directory` and `mida_directory` in `figure3\demo_2d_spiral_simulation.m`.

    src_directory = 'path-to-this-package';
    mida_directory = 'path-to-this-package\thirdparty\MIDA_v1.0\MIDA_v1_voxels';
 
Run `figure3\demo_batch_2d_spiral_simulation.m`

Run `figure3\demo_generate_figure3.m` to generate Figure 3.
 
## Figures 5 and 6 (human axial spiral spin-echo imaging)

Update variables in `demo_non_cartesian_recon_human_axial.m`.

    src_directory          = 'path-to-this-package';
    ismrmrd_directory      = 'path-to-ISMRMRD-package';
    ismrmrd_noise_fullpath = 'path-to-dataset\noise_meas_MID00275_FID03658_se_spiral_1102_ax_s24.h5';
    ismrmrd_data_fullpath  = 'path-to-dataset\meas_MID00275_FID03658_se_spiral_1102_ax_s24.h5';
    siemens_dat_fullpath   = 'path-to-dataset\meas_MID00275_FID03658_se_spiral_1102_ax_s24.dat';
    B0map_fullpath         = 'path-to-dataset\B0map_nlinv_min1.0e-06_axial.mat';

Run `demo_non_cartesian_recon_human_axial.m`.

Run `demo_cartesian_recon_SE_axial.m`.

Run `figure5\demo_generate_figure5.m` and `figure6\demo_generate_figure6.m` to generate Figures 5 and 6, respectively.

## Figures 7 and 8 (human sagittal spiral spin-echo imaging)

Update variables in `demo_non_cartesian_recon_human_sagittal.m`.

    src_directory          = 'path-to-this-package';
    ismrmrd_directory      = 'path-to-ISMRMRD-package';
    ismrmrd_noise_fullpath = 'path-to-dataset\noise_meas_MID00260_FID03643_se_spiral_1102_sag_s24.h5';
    ismrmrd_data_fullpath  = 'path-to-dataset\meas_MID00260_FID03643_se_spiral_1102_sag_s24.h5';
    siemens_dat_fullpath   = 'path-to-dataset\meas_MID00260_FID03643_se_spiral_1102_sag_s24.dat';
    B0map_fullpath         = 'path-to-dataset\B0map_nlinv_min1.0e-06_sagittal_ro0.mat';

Run `demo_non_cartesian_recon_human_sagittal.m`.

Run `demo_cartesian_recon_SE_sagittal.m`.

Run `figure7\demo_generate_figure7.m` and `figure8\demo_generate_figure8.m` to generate Figures 7 and 8, respectively.
