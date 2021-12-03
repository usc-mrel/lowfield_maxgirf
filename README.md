# lowfield_maxgirf

This repository contains code and datasets for

"MaxGIRF: Image Reconstruction Incorporating Concomitant
Field and Gradient Impulse Response Function Effects"

This code is distributed under the BSD license.

Nam Gyun Lee, University of Southern California, Dec 2021.

Download ISMRMD datasets:

## Figure 3 (numerical simulation)

1. Update variables `src_directory` and `mida_directory` in `demo_2d_spiral_simulation.m`.

    src_directory = 'E:\lowfield_maxgirf';
    mida_directory = 'E:\lowfield_maxgirf\thirdparty\MIDA_v1.0\MIDA_v1_voxels';
 
2. Run `demo_batch_2d_spiral_simulation.m`.

3. Run `demo_generate_figure3.m`.

 
## Figures 5 and 6 (human axial spiral spin-echo imaging)

* Run `demo_batch_maxgirf_spiral_se_human_axial.m`

## Figures 7 and 8 (human sagittal spiral spin-echo imaging)

* Run `demo_batch_maxgirf_spiral_se_human_sagittal.m`
