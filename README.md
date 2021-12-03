# lowfield_maxgirf

This repository contains code and datasets for

"MaxGIRF: Image Reconstruction Incorporating Concomitant
Field and Gradient Impulse Response Function Effects".

This code is distributed under the BSD license.

Nam Gyun Lee, University of Southern California, Dec 2021.

Download ISMRMD datasets:

Directory structure for the datasets:

e.g.) |---20201102_NV_brain
           |---d20201102_NV_brain.mat
           |---h5
		        |---meas_MID00260_FID03643_se_spiral_1102_sag_s24.h5
				|---meas_MID00275_FID03658_se_spiral_1102_ax_s24.h5
		   |---noise
		        |---noise_meas_MID00260_FID03643_se_spiral_1102_sag_s24.h5
		        |---noise_meas_MID00275_FID03658_se_spiral_1102_ax_s24.h5
 

## Figure 3 (numerical simulation)

Update variables `src_directory` and `mida_directory` in `figure3\demo_2d_spiral_simulation.m`.

    src_directory = 'path-to-this-package';
    mida_directory = 'path-to-this-package\thirdparty\MIDA_v1.0\MIDA_v1_voxels';
 
Run `figure3\demo_batch_2d_spiral_simulation.m`

Run `figure3\demo_generate_figure3.m`
 
## Figures 5 and 6 (human axial spiral spin-echo imaging)

Update variables `src_directory` and `ismrmrd_directory` in `demo_maxgirf_spiral_se_human.m`.

    src_directory = 'path-to-this-package';
    ismrmrd_directory = 'path-to-ISMRMRD-package';

Update variables `B0map_fullpath` and `data_directory` in `demo_batch_maxgirf_spiral_se_human_axial.m`.

    B0map_fullpath = 'path-to-this-package\B0map_nlinv_min1.0e-06_axial.mat';
    data_directory = 'path-to-human-dataset\20201102_NV_brain';

Run `demo_batch_maxgirf_spiral_se_human_axial.m`.

Run `figures5&6\demo_generate_figure5.m` and `figures5&6\demo_generate_figure6.m` to generate Figures 5 and 6, respectively.

## Figures 7 and 8 (human sagittal spiral spin-echo imaging)

Update variables `src_directory` and `ismrmrd_directory` in `demo_maxgirf_spiral_se_human.m`.

    src_directory = 'E:\lowfield_maxgirf';
    ismrmrd_directory = 'D:\ismrmrd\ismrmrd';

Update variables `B0map_fullpath` and `data_directory` in `demo_batch_maxgirf_spiral_se_human_sagittal.m`.

    B0map_fullpath = 'E:\lowfield_maxgirf\B0map_nlinv_min1.0e-06_sagittal.mat';
    data_directory = 'D:\lowfield\NHLBI\data\20201102_NV_brain';

Run `demo_batch_maxgirf_spiral_se_human_sagittal.m`.

Run `figures7&8\demo_generate_figure7.m` and `figures7&8\demo_generate_figure8.m` to generate Figures 7 and 8, respectively.
