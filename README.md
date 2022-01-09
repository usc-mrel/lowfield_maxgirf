# lowfield_maxgirf

This repository contains the code and datasets for
**"MaxGIRF: Image Reconstruction Incorporating Concomitant
Field and Gradient Impulse Response Function Effects"**, by Nam G. Lee, Rajiv Ramasawmy, Yongwan Lim, Adrienne E. Campbell-Washburn, and Krishna S. Nayak.

This code is distributed under the BSD license.

Nam Gyun Lee, University of Southern California, Dec 2021.

## ISMRMRD datasets
Example human and phantom datasets can be found on Zenodo:

- Human: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5830910.svg)](https://doi.org/10.5281/zenodo.5830910)
- Phantom:

* **Axial Cartesian spin-echo dataset**
  - [meas_MID00273_FID03656_se_15b130_tra.h5]

* **Axial spiral spin-echo dataset** 
  - [meas_MID00275_FID03658_se_spiral_1102_ax_s24.h5](https://drive.google.com/file/d/1M5bMNL2bWOsEqaKKpHBPQLHVSesq-Lx2/view?usp=sharing)

* **Axial Cartesian GRE datasets for static off-resonance map estimation**

  - [meas_MID00276_FID03659_gre_TE1.h5](https://drive.google.com/file/d/1oNQJP_fau6dZoS5EXGtUk8yQ75Bkb7GV/view?usp=sharing)  
  - [meas_MID00277_FID03660_gre_TE2.h5](https://drive.google.com/file/d/1T5czOtlOezgHWib8Teoz1JBZE_bMtg4F/view?usp=sharing)
  - [meas_MID00278_FID03661_gre_TE3.h5](https://drive.google.com/file/d/1zhbBxD9RJ0v4JIcsmbH6_Rv0NfkSDggb/view?usp=sharing) 
  - [meas_MID00279_FID03662_gre_TE4.h5](https://drive.google.com/file/d/1hE9suU9RN0c8LdPuwzgOt6N4M_0y9daV/view?usp=sharing)
  - [meas_MID00280_FID03663_gre_TE5.h5](https://drive.google.com/file/d/1Stp4XnRI91sbJ2eBa6VLALIoKrKnKsfY/view?usp=sharing)

* **Axial noise only datasets**
  - [noise_meas_MID00273_FID03656_se_15b130_tra.h5](https://drive.google.com/file/d/1WXC0YfERL8yTGisEU1suVL-EZel8DgQy/view?usp=sharing)
  - [noise_meas_MID00275_FID03658_se_spiral_1102_ax_s24.h5](https://drive.google.com/file/d/1LpWtlNOvPFWWV_6lCJ2GMwjP0AyM1qt8/view?usp=sharing)
  - [noise_meas_MID00276_FID03659_gre_TE1.h5](https://drive.google.com/file/d/1EGSKo5qSrLfusbYnECRTG3Jqib-DcPMp/view?usp=sharing)
  - [noise_meas_MID00276_FID03660_gre_TE2.h5](https://drive.google.com/file/d/1hTlesOJvDV6R4aQduOP_1T9QdT-sKGcQ/view?usp=sharing)
  - [noise_meas_MID00276_FID03661_gre_TE3.h5](https://drive.google.com/file/d/1arMDTIoAMz3-_otCQlQPsFhLKchEdG5k/view?usp=sharing)
  - [noise_meas_MID00276_FID03662_gre_TE4.h5](https://drive.google.com/file/d/1WL_zA2CmvMHc8DjriCUU9D-A61ZicfkK/view?usp=sharing)
  - [noise_meas_MID00276_FID03663_gre_TE5.h5](https://drive.google.com/file/d/1SQaDExFkGq9fUoPZnJoorW0jWbNSpzc9/view?usp=sharing)

* **Sagittal Cartesian spin-echo dataset**
  - [meas_MID00258_FID03641_se_15b130.h5](https://drive.google.com/file/d/1PK79_QYW82A33nqXM2sTkXpb11HsKxrD/view?usp=sharing)

* **Sagittal spiral spin-echo dataset**
  - [meas_MID00260_FID03643_se_spiral_1102_sag_s24.h5](https://drive.google.com/file/d/1NOh64QtBmbyImTiHzzqwuV_GXiyfvvAR/view?usp=sharing)

* **Sagittal Cartesian GRE datasets for static off-resonance map estimation**
  - [meas_MID00261_FID03644_gre_TE1.h5](https://drive.google.com/file/d/11wJrXucCl9j7Q1LXP3SClGrycmmaatWu/view?usp=sharing)
  - [meas_MID00262_FID03645_gre_TE2.h5](https://drive.google.com/file/d/1aep6XWc8Ijjjw2nFGCCagTu2wLzLFP6W/view?usp=sharing)
  - [meas_MID00263_FID03646_gre_TE3.h5](https://drive.google.com/file/d/1d3-0V6vniSekO160inIsLdCujF20UjjS/view?usp=sharing)
  - [meas_MID00264_FID03647_gre_TE4.h5](https://drive.google.com/file/d/1xMWlx6UWcBdMAcqjkT49vLj-3JYHYGJO/view?usp=sharing)
  - [meas_MID00265_FID03648_gre_TE5.h5](https://drive.google.com/file/d/1TySk-_X00Sx7cXTfWDRc7YeH1U2P-SNH/view?usp=sharing)

* **Sagittal noise only datasets**
  - [noise_meas_MID00258_FID03641_se_15b130.h5](https://drive.google.com/file/d/1flg3lO7K3dxnZOigjV7xIkifH1UE--kO/view?usp=sharing)
  - [noise_meas_MID00260_FID03643_se_spiral_1102_sag_s24.h5](https://drive.google.com/file/d/1ZCwd5p3zht53_2bskfFNuAZnXUXfOIHl/view?usp=sharing)
  - [noise_meas_MID00261_FID03644_gre_TE1.h5](https://drive.google.com/file/d/1rO0JGqbDSCRjkSE6xv0bInTgfWxkdYxl/view?usp=sharing)
  - [noise_meas_MID00261_FID03645_gre_TE2.h5](https://drive.google.com/file/d/1ZkfxvxBhUhZVFkKcS0FdMBPkmI8--xTi/view?usp=sharing)
  - [noise_meas_MID00261_FID03646_gre_TE3.h5](https://drive.google.com/file/d/1w-ZTvYwrdeMe5KtCgPzVOtrAnRJEVdAE/view?usp=sharing)
  - [noise_meas_MID00261_FID03647_gre_TE4.h5](https://drive.google.com/file/d/1U3joLbogiJZEySx-08UfXzKVz1Hwkq_9/view?usp=sharing)
  - [noise_meas_MID00261_FID03648_gre_TE5.h5](https://drive.google.com/file/d/17LRHCkYv5OL-7agkSe336m90mgp49UOR/view?usp=sharing)

* **Mat file containing reconstructed Cartesian spin-echo images**
  - [d20201102_NV_brain.mat](https://drive.google.com/file/d/1yU42YylEXz8YH_UbEsBLQceD0ywTYFfb/view?usp=sharing)

Once datasets are downloaded, organize datasets in this structure:
 
     |---path-to-dataset (e.g., D:\lowfield\NHLBI\data\20201102_NV_brain)
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
    data_directory = 'path-to-dataset';

Run `demo_batch_maxgirf_spiral_se_human_axial.m`.

Run `figures5&6\demo_generate_figure5.m` and `figures5&6\demo_generate_figure6.m` to generate Figures 5 and 6, respectively.

## Figures 7 and 8 (human sagittal spiral spin-echo imaging)

Update variables `src_directory` and `ismrmrd_directory` in `demo_maxgirf_spiral_se_human.m`.

    src_directory = 'path-to-this-package';
    ismrmrd_directory = 'path-to-ISMRMRD-package';

Update variables `B0map_fullpath` and `data_directory` in `demo_batch_maxgirf_spiral_se_human_sagittal.m`.

    B0map_fullpath = 'path-to-this-package\B0map_nlinv_min1.0e-06_sagittal.mat';
    data_directory = 'path-to-dataset';

Run `demo_batch_maxgirf_spiral_se_human_sagittal.m`.

Run `figures7&8\demo_generate_figure7.m` and `figures7&8\demo_generate_figure8.m` to generate Figures 7 and 8, respectively.
