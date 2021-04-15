% demo_siemens_coordinate_systems.m
% Written by Namgyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 05/18/2020, Last modified: 04/11/2021

%% Clean slate
close all; clear; clc;

%% Set directory names
%--------------------------------------------------------------------------
% Package directory
%--------------------------------------------------------------------------
src_directory = 'E:\NHLBI-USC-share';

%--------------------------------------------------------------------------
% ISMRMRD directory
%--------------------------------------------------------------------------
ismrmrd_directory = 'D:\ismrmrd\ismrmrd';

%% Add paths
addpath(genpath(src_directory));
addpath(genpath(ismrmrd_directory));

%% Notation
% We denote X, Y, Z as the horizontal, vertical, and through-magnet axes,
% respectively in the physical coordinate system, and the corresponding
% coordinates as x, y, and z. We also denote U, V, and W as the readout,
% phase encoding, and slice axes, respectively in the logical coordinate
% system, and the corresponding coordinates as u, v, and w.

%% Define .dat filenames
% data_directory = 'D:\lowfield\NHLBI\data\20200506_NIST_phantom_rot_HFS';
% output_filename = '20200506_NIST_phantom_rot_HFS';
% % sagittal
% dat_list{1} = 'meas_MID00095_FID69440_crt_ssfp_sag_0deg';
% dat_list{2} = 'meas_MID00096_FID69441_crt_ssfp_sag_15deg';
% dat_list{3} = 'meas_MID00097_FID69442_crt_ssfp_sag_30deg';
% dat_list{4} = 'meas_MID00112_FID69457_crt_ssfp_sag_44deg';
% dat_list{5} = 'meas_MID00098_FID69443_crt_ssfp_sag_45deg_peSwitch';
% dat_list{6} = 'meas_MID00099_FID69444_crt_ssfp_sag_m15deg';
% dat_list{7} = 'meas_MID00100_FID69445_crt_ssfp_sag_m30deg';
% dat_list{8} = 'meas_MID00113_FID69458_crt_ssfp_sag_m44deg';
% dat_list{9} = 'meas_MID00101_FID69446_crt_ssfp_sag_m45deg_peSwitch';
% 
% % coronal
% dat_list{10} = 'meas_MID00102_FID69447_crt_ssfp_cor_0deg';
% dat_list{11} = 'meas_MID00104_FID69449_crt_ssfp_cor_15deg';
% dat_list{12} = 'meas_MID00105_FID69450_crt_ssfp_cor_30deg';
% dat_list{13} = 'meas_MID00107_FID69452_crt_ssfp_cor_44deg';
% dat_list{14} = 'meas_MID00106_FID69451_crt_ssfp_cor_45deg_peSwitch';
% dat_list{15} = 'meas_MID00108_FID69453_crt_ssfp_cor_m15deg';
% dat_list{16} = 'meas_MID00109_FID69454_crt_ssfp_cor_m30deg';
% dat_list{17} = 'meas_MID00110_FID69455_crt_ssfp_cor_m44deg';
% dat_list{18} = 'meas_MID00111_FID69456_crt_ssfp_cor_m45deg_peSwitch';
% 
% % axial
% dat_list{19} = 'meas_MID00087_FID69432_crt_ssfp_ax_0deg';
% dat_list{20} = 'meas_MID00088_FID69433_crt_ssfp_ax_15deg';
% dat_list{21} = 'meas_MID00089_FID69434_crt_ssfp_ax_30deg';
% dat_list{22} = 'meas_MID00091_FID69436_crt_ssfp_ax_45deg';
% dat_list{23} = 'meas_MID00092_FID69437_crt_ssfp_ax_m15deg';
% dat_list{24} = 'meas_MID00093_FID69438_crt_ssfp_ax_m30deg';
% dat_list{25} = 'meas_MID00094_FID69439_crt_ssfp_ax_m45deg';

%%
% data_directory = 'D:\lowfield\NHLBI\data\20200528_NIST_phantom_rot_ipshift';
% output_filename = '20200528_NIST_phantom_rot_ipshift';
% % axial
% dat_list{1} = 'meas_MID00423_FID69751_crt_ssfp_ax_0deg';
% dat_list{2} = 'meas_MID00424_FID69752_crt_ssfp_ax_180deg';
% dat_list{3} = 'meas_MID00425_FID69753_crt_ax_L3_P6_H0_0deg';
% dat_list{4} = 'meas_MID00426_FID69754_crt_ax_L3_P6_H0_30deg';
% dat_list{5} = 'meas_MID00427_FID69755_crt_ax_L3_P6_H0_m30deg';
% 
% % sagittal
% dat_list{6} = 'meas_MID00428_FID69762_crt_sag_L0_P6_H5_0deg';
% dat_list{7} = 'meas_MID00429_FID69763_crt_sag_L0_P6_H5_30deg';
% dat_list{8} = 'meas_MID00430_FID69764_crt_sag_L0_P6_H5_m30deg';
% 
% % coronal
% dat_list{9} = 'meas_MID00431_FID69771_crt_cor_L3_P0_H5_0deg';
% dat_list{10} = 'meas_MID00432_FID69772_crt_cor_L3_P0_H5_30deg';
% dat_list{11} = 'meas_MID00433_FID69773_crt_cor_L3_P0_H5_m30deg';

%%
% data_directory = 'D:\lowfield\NHLBI\data\20200506_NIST_phantom_rot_FFP';
% output_filename = '20200506_NIST_phantom_rot_FFP';
% % sagittal
% dat_list{1} = 'meas_MID00153_FID69498_crt_ssfp_sag_0deg';
% dat_list{2} = 'meas_MID00154_FID69499_crt_ssfp_sag_30deg';
% dat_list{3} = 'meas_MID00155_FID69500_crt_ssfp_sag_m30deg';
% 
% % coronal
% dat_list{4} = 'meas_MID00156_FID69501_crt_ssfp_cor_0deg';
% dat_list{5} = 'meas_MID00157_FID69502_crt_ssfp_cor_30deg';
% dat_list{6} = 'meas_MID00158_FID69503_crt_ssfp_cor_m30deg';
% 
% % axial
% dat_list{7} = 'meas_MID00150_FID69495_crt_ssfp_ax_0deg';
% dat_list{8} = 'meas_MID00151_FID69496_crt_ssfp_ax_30deg';
% dat_list{9} = 'meas_MID00152_FID69497_crt_ssfp_ax_m30deg';

%%
data_directory = 'D:\lowfield\NHLBI\data\20200608_NIST_Phantom_rot2pi';
output_filename = '20200608_NIST_Phantom_rot2pi';
dat_list{1} = 'meas_MID00236_FID70326_crt_ssfp_ax_0deg';
dat_list{2} = 'meas_MID00237_FID70327_crt_ax_L3_P6_H0_0deg_pe_AP';
dat_list{3} = 'meas_MID00238_FID70328_crt_ax_L3_P6_H0_PE_RL';
dat_list{4} = 'meas_MID00239_FID70329_crt_ax_L3_P6_H0_90deg_pe_RL';
dat_list{5} = 'meas_MID00240_FID70330_crt_ax_L3_P6_H0_45deg_pe_AP';
dat_list{6} = 'meas_MID00241_FID70331_crt_ax_L3_P6_H0_46deg_peSwitch_RL';
dat_list{7} = 'meas_MID00242_FID70332_crt_ax_L3_P6_H0_134deg_peSwitch_RL';
dat_list{8} = 'meas_MID00243_FID70333_crt_ax_L3_P6_H0_135deg_peSwitch_PA';
dat_list{9} = 'meas_MID00244_FID70334_crt_ax_L3_P6_H0_180deg_peSwitch_PA';
dat_list{10} = 'meas_MID00245_FID70335_crt_ax_L3_P6_H0_m45deg_pe_AP';
dat_list{11} = 'meas_MID00246_FID70336_crt_ax_L3_P6_H0_m46deg_peSwitch_LR';
dat_list{12} = 'meas_MID00247_FID70337_crt_ax_L3_P6_H0_m134deg_peSwitch_LR';
dat_list{13} = 'meas_MID00248_FID70338_crt_ax_L3_P6_H0_m134deg_peSwitch_PA';
dat_list{14} = 'meas_MID00249_FID70339_crt_ax_L3_P6_H0_m180deg_peSwitch_PA';
dat_list{15} = 'meas_MID00250_FID70340_crt_sag_L0_P6_H5_0deg_pe_AP';
dat_list{16} = 'meas_MID00251_FID70341_crt_sag_L0_P6_H5_44deg_pe_AP';
dat_list{17} = 'meas_MID00252_FID70342_crt_sag_L0_P6_H5_45deg_peSwitch_HF';
dat_list{18} = 'meas_MID00253_FID70343_crt_sag_L0_P6_H5_135deg_peSwitch_HF';
dat_list{19} = 'meas_MID00254_FID70344_crt_sag_L0_P6_H5_136deg_peSwitch_PA';
dat_list{20} = 'meas_MID00255_FID70345_crt_sag_L0_P6_H5_180deg_peSwitch_PA';
dat_list{21} = 'meas_MID00256_FID70346_crt_sag_L0_P6_H5_m44deg_pe_AP';
dat_list{22} = 'meas_MID00257_FID70347_crt_sag_L0_P6_H5_m45deg_peSwitch_FH';
dat_list{23} = 'meas_MID00260_FID70350_crt_sag_L0_P6_H5_m135deg_peSwitch_FH';
dat_list{24} = 'meas_MID00261_FID70351_crt_sag_L0_P6_H5_m136deg_peSwitch_PA';
dat_list{25} = 'meas_MID00262_FID70352_crt_sag_L0_P6_H5_m180deg_peSwitch_PA';
dat_list{26} = 'meas_MID00263_FID70353_crt_cor_L3_P0_H5_0deg_pe_RL';
dat_list{27} = 'meas_MID00264_FID70354_crt_cor_L3_P0_H5_44deg_pe_RL';
dat_list{28} = 'meas_MID00265_FID70355_crt_cor_L3_P0_H5_45deg_peSwitch_FH';
dat_list{29} = 'meas_MID00266_FID70356_crt_cor_L3_P0_H5_135deg_peSwitch_FH';
dat_list{30} = 'meas_MID00267_FID70357_crt_cor_L3_P0_H5_136deg_peSwitch_LR';
dat_list{31} = 'meas_MID00268_FID70358_crt_cor_L3_P0_H5_180deg_peSwitch_LR';
dat_list{32} = 'meas_MID00269_FID70359_crt_cor_L3_P0_H5_m44deg_pe_RL';
dat_list{33} = 'meas_MID00270_FID70360_crt_cor_L3_P0_H5_m45deg_peSwitch_HF';
dat_list{34} = 'meas_MID00271_FID70361_crt_cor_L3_P0_H5_m135deg_peSwitch_HF';
dat_list{35} = 'meas_MID00272_FID70362_crt_cor_L3_P0_H5_m136deg_peSwitch_LR';
dat_list{36} = 'meas_MID00273_FID70363_crt_cor_L3_P0_H5_m180deg_peSwitch_LR';

nr_files = length(dat_list);

%% Define directory names
h5_directory    = fullfile(data_directory, 'h5');
noise_directory = fullfile(data_directory, 'noise');
dicom_directory = fullfile(data_directory, 'dicom');
dicom_directory_info = dir(dicom_directory);

%% Define an output directory
output_directory = fullfile(pwd, output_filename);
mkdir(output_directory);

%% Define constants
SAGITTAL   = 0; % Patient axis perpendicular to the sagittal plane
CORONAL    = 1; % Patient axis perpendicular to the coroal plane
TRANSVERSE = 2; % Patient axis perpendicular to the transvers plane

%% Process datasets
for k = 1:nr_files
    dat_filename = dat_list{k};
    if isempty(dat_filename)
        continue;
    end

    under_loc = strfind(dat_filename, '_');
    deg_loc   = strfind(dat_filename, 'deg');
    prefix = dat_filename(under_loc(3)+1:end);

    %% Read a Siemens dat file
    twix = mapVBVD(fullfile(data_directory, dat_filename));
    if length(twix) > 1
        twix = twix{end};
    end

    %% Get a slice normal vector
    %----------------------------------------------------------------------
    % dNormalSag: Sagittal component of a slice normal vector (in PCS)
    %----------------------------------------------------------------------
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dSag')
        dNormalSag = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dSag;
    else
        dNormalSag = 0;
    end

    %----------------------------------------------------------------------
    % dNormalCor: Coronal component of a slice normal vector (in PCS)
    %----------------------------------------------------------------------
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dCor')
        dNormalCor = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dCor;
    else
        dNormalCor = 0;
    end

    %----------------------------------------------------------------------
    % dNormalTra: Transverse component of a slice normal vector (in PCS)
    %----------------------------------------------------------------------
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dTra')
        dNormalTra = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dTra;
    else
        dNormalTra = 0;
    end

    %----------------------------------------------------------------------
    % dRotAngle: Slice rotation angle ("swap Fre/Pha")
    %----------------------------------------------------------------------
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}, 'dInPlaneRot')
        dRotAngle = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dInPlaneRot; % [rad]
    else
        dRotAngle = 0; % [rad]
    end

    %% Get a slice offset in PCS
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}, 'sPosition')
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition, 'dSag')
            sag_offset = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dSag; % [mm]
        else
            sag_offset = 0; % [mm]
        end
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition, 'dCor')
            cor_offset = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dCor; % [mm]
        else
            cor_offset = 0; % [mm]
        end
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition, 'dTra')
            tra_offset = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dTra; % [mm]
        else
            tra_offset = 0; % [mm]
        end
    else
        sag_offset = 0; % [mm]
        cor_offset = 0; % [mm]
        tra_offset = 0; % [mm]
    end
    PCS_offset = [sag_offset; cor_offset; tra_offset] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Calculate the main orientation of a lice
    main_orientation = fGSLClassOri(dNormalSag, dNormalCor, dNormalTra);
    switch main_orientation
        case SAGITTAL
            image_ori = 'sag';
        case CORONAL
            image_ori = 'cor';
        case TRANSVERSE
            image_ori = 'tra';
        otherwise
    end
    
    %% Read an ismrmrd file
    dset = ismrmrd.Dataset(fullfile(h5_directory, sprintf('%s.h5', dat_filename)), 'dataset');

    %% Read some fields from the XML header
    header = ismrmrd.xml.deserialize(dset.readxml);

    Nkx = header.encoding.encodedSpace.matrixSize.x;
    Nky = header.encoding.encodedSpace.matrixSize.y;
    Nkz = header.encoding.encodedSpace.matrixSize.z;
    Nx  = header.encoding.reconSpace.matrixSize.x;
    Ny  = header.encoding.reconSpace.matrixSize.y;
    Nz  = header.encoding.reconSpace.matrixSize.z;
    Nc  = header.acquisitionSystemInformation.receiverChannels;

    patient_position  = header.measurementInformation.patientPosition;
    fieldOfView_mm(1) = header.encoding.reconSpace.fieldOfView_mm.x; % RO
    fieldOfView_mm(2) = header.encoding.reconSpace.fieldOfView_mm.y; % PE
    fieldOfView_mm(3) = header.encoding.reconSpace.fieldOfView_mm.z; % SS

    %% Parse the ISMRMRD header
    tic; fprintf('Parsing the ISMRMRD header... ');
    raw_data = dset.readAcquisition(); % read all the acquisitions
    fprintf('done! (%6.4f sec)\n', toc);

    %% Calculate a scaling matrix [m]
    scaling_matrix = [fieldOfView_mm(1) / Nx 0 0; 0 fieldOfView_mm(2) / Ny 0; 0 0 fieldOfView_mm(3) / Nz] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Calculate a rotation matrix from RCS to GCS [r,c,s] <=> [PE,RO,SS] 
    rotMatrixRCSToGCS = [0    1    0 ; % [PE]   [0 1 0] * [r]
                         1    0    0 ; % [RO] = [1 0 0] * [c]
                         0    0    1]; % [SS]   [0 0 1] * [s]

    %% Calculate a rotation matrix from GCS to PCS
    [rotMatrixGCSToPCS,PE_sign,RO_sign] = calcMatrixGCSToPCS(dNormalSag, dNormalCor, dNormalTra, dRotAngle);
    if strcmp(header.encoding.trajectory, 'spiral')
        rotMatrixGCSToPCS(:,1) = PE_sign * rotMatrixGCSToPCS(:,1);
        rotMatrixGCSToPCS(:,2) = RO_sign * rotMatrixGCSToPCS(:,2);
    end
    
    %% Calculate a rotation matrix from PCS to DCS
    rotMatrixPCSToDCS = calcMatrixPCSToDCS(patient_position);

    %% Calculate a rotation matrix from GCS to DCS
    rotMatrixGCSToDCS = rotMatrixPCSToDCS * rotMatrixGCSToPCS;

    %% Get a rotation matrix from GCS to PCS (ISMRMRD format)
    rotMatrixGCSToPCS_ismrmrd = double([raw_data.head.phase_dir(:,1) raw_data.head.read_dir(:,1) raw_data.head.slice_dir(:,1)]);

    %% Get a slice offset from ISMRMRD format
    sag_offset_ismrmrd = raw_data.head.position(1,1); % [mm]
    cor_offset_ismrmrd = raw_data.head.position(2,1); % [mm]
    tra_offset_ismrmrd = raw_data.head.position(3,1); % [mm]

    %% Display slice information
    fprintf('======================= SLICE INFORMATION =======================\n');
    fprintf('---------------------- From Siemens raw data --------------------\n');
    fprintf('                   [sag]   %10.5f [mm]\n', sag_offset);
    fprintf('slice offset     : [cor] = %10.5f [mm]\n', cor_offset);
    fprintf('                   [tra]   %10.5f [mm]\n', tra_offset);
    fprintf('-----------------------------------------------------------------\n');
    fprintf('                   [sag]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToPCS(1,1), rotMatrixGCSToPCS(1,2), rotMatrixGCSToPCS(1,3));
    fprintf('rotMatrixGCSToPCS: [cor] = [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToPCS(2,1), rotMatrixGCSToPCS(2,2), rotMatrixGCSToPCS(2,3));
    fprintf('                   [tra]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToPCS(3,1), rotMatrixGCSToPCS(3,2), rotMatrixGCSToPCS(3,3));
    fprintf('-----------------------------------------------------------------\n');
    fprintf('                   [sag]   [%10.5f %10.5f %10.5f]\n', rotMatrixPCSToDCS(1,1), rotMatrixPCSToDCS(1,2), rotMatrixPCSToDCS(1,3));
    fprintf('rotMatrixPCSToDCS: [cor] = [%10.5f %10.5f %10.5f]\n', rotMatrixPCSToDCS(2,1), rotMatrixPCSToDCS(2,2), rotMatrixPCSToDCS(2,3));
    fprintf('                   [tra]   [%10.5f %10.5f %10.5f]\n', rotMatrixPCSToDCS(3,1), rotMatrixPCSToDCS(3,2), rotMatrixPCSToDCS(3,3));
    fprintf('-----------------------------------------------------------------\n');
    fprintf('                   [sag]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToDCS(1,1), rotMatrixGCSToDCS(1,2), rotMatrixGCSToDCS(1,3));
    fprintf('rotMatrixGCSToDCS: [cor] = [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToDCS(2,1), rotMatrixGCSToDCS(2,2), rotMatrixGCSToDCS(2,3));
    fprintf('                   [tra]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToDCS(3,1), rotMatrixGCSToDCS(3,2), rotMatrixGCSToDCS(3,3));
    fprintf('---------------------- From ISMRMRD header ----------------------\n');
    fprintf('                   [sag]   %10.5f [mm]\n', sag_offset_ismrmrd);
    fprintf('slice offset     : [cor] = %10.5f [mm]\n', cor_offset_ismrmrd);
    fprintf('                   [tra]   %10.5f [mm]\n', tra_offset_ismrmrd);
    fprintf('-----------------------------------------------------------------\n');
    fprintf('                   [sag]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToPCS_ismrmrd(1,1), rotMatrixGCSToPCS_ismrmrd(1,2), rotMatrixGCSToPCS_ismrmrd(1,3));
    fprintf('rotMatrixGCSToPCS: [cor] = [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToPCS_ismrmrd(2,1), rotMatrixGCSToPCS_ismrmrd(2,2), rotMatrixGCSToPCS_ismrmrd(2,3));
    fprintf('                   [tra]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToPCS_ismrmrd(3,1), rotMatrixGCSToPCS_ismrmrd(3,2), rotMatrixGCSToPCS_ismrmrd(3,3));
    fprintf('=================================================================\n');

    %% Calculate spatial coordinates in PCS (LPH DICOM coordinates)
    r_range = (-floor(Nx/2):ceil(Nx/2)-1).';
    c_range = (-floor(Ny/2):ceil(Ny/2)-1).';
    s_range = (-floor(Nz/2):ceil(Nz/2)-1).';
    [r,c,s] = ndgrid(r_range, c_range, s_range); % 3 x 1

    xyz = rotMatrixGCSToPCS * rotMatrixRCSToGCS * scaling_matrix * cat(2, r(:), c(:), s(:)).';
    x = reshape(xyz(1,:), [Nx Ny]) + PCS_offset(1); % [m]
    y = reshape(xyz(2,:), [Nx Ny]) + PCS_offset(2); % [m]
    z = reshape(xyz(3,:), [Nx Ny]) + PCS_offset(3); % [m]

    %% Read k-space data
    %----------------------------------------------------------------------
    % Read k-space data of a target image
    %----------------------------------------------------------------------
    raw_data = dset.readAcquisition();
    is_noise = raw_data.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
    first_scan = find(is_noise == 0, 1, 'first');
    meas = raw_data.select(first_scan:raw_data.getNumber);

    %----------------------------------------------------------------------
    % Prepare k-space data
    %----------------------------------------------------------------------
    kspace = zeros(Nkx, Nky, Nc, 'double');
    acqs = find((meas.head.idx.slice == 0));
    for p = 1:length(acqs)
        ky = meas.head.idx.kspace_encode_step_1(acqs(p)) + 1;
        kspace(:,ky,:) = meas.data{acqs(p)};
    end

    %% Flip k-space
    if RO_sign == -1
        kspace = flip(kspace,1);
    end
    if PE_sign == -1
        kspace = flip(kspace,2);
    end

    %% Reconstruct in x
    imc = fftshift(fft(ifftshift(kspace, 1), [], 1), 1);

    %% Remove oversampling
    idx1_range = (-floor(Nx/2):ceil(Nx/2)-1).' + floor(Nkx/2) + 1;
    imc = imc(idx1_range,:,:);

    %% Reconstruct in y
    imc = fftshift(fft(ifftshift(imc, 2), [], 2), 2);

    %% Calculate coil sensitivity maps
    %----------------------------------------------------------------------
    % IFFT to k-space (k-space <=> image-space)
    %----------------------------------------------------------------------
    kspace_cartesian = ifft2c(imc); % Nx x Ny x Nc

    %----------------------------------------------------------------------
    % Calculate the calibration region of k-space
    %----------------------------------------------------------------------
    cal_shape = [32 32];
    cal_data = crop(kspace_cartesian, [cal_shape Nc]);
    cal_data = bsxfun(@times, cal_data, hamming(cal_shape(1)) * hamming(cal_shape(2)).');

    %----------------------------------------------------------------------
    % Calculate coil sensitivity maps (Nx x Ny x Nc)
    %----------------------------------------------------------------------
    cal_im = fft2c(zpad(cal_data, [Nx Ny Nc]));
    csm = ismrm_estimate_csm_walsh(cal_im);

    %% Perform optimal coil combination
    im = sum(imc .* conj(csm), 3);

    %% Sort out dicom folders (without GNC)
    % Dicoms include ND (non-distortion correction) and Distortion corrected images (no label)
    for idx1 = 1:length(dicom_directory_info)
        if ~isempty(strfind(dicom_directory_info(idx1).name, prefix)) && ...
            dicom_directory_info(idx1).isdir && ...
            isempty(strfind(dicom_directory_info(idx1).name, 'ND')) && ...
            isempty(strfind(dicom_directory_info(idx1).name, 'RR'))
            sub_dicom_directory = fullfile(dicom_directory_info(idx1).folder, dicom_directory_info(idx1).name);
            file_info = dir(sub_dicom_directory);
            for idx2 = 1:length(file_info)
                if strfind(file_info(idx2).name, 'IMA')
                    dicom_filename = fullfile(file_info(idx2).folder, file_info(idx2).name);
                end
            end
        end
    end

    %% Read a dicom file
    dcm_info = dicominfo(dicom_filename);
    im_dcm = double(dicomread(dcm_info));

    %% Calculate spatial coordinates for a dicom image in LPH
    iop             = dcm_info.ImageOrientationPatient;
    ipp             = dcm_info.ImagePositionPatient; % [mm]
    pixel_spacing   = dcm_info.PixelSpacing; % [mm]
    slice_thickness = dcm_info.SliceThickness; % [mm]

    scaling_dicom = [pixel_spacing(1) 0 0; 0 pixel_spacing(2) 0; 0 0 slice_thickness] * 1e-3; % [mm] * [m/1e3mm] => [m]

    R_rcs2lph = cat(2, iop(4:6), iop(1:3), cross(iop(1:3), iop(4:6)));

    r_range = (0:Nx-1).';
    c_range = (0:Ny-1).';
    s_range = (0:Nz-1).';
    [r,c,s] = ndgrid(r_range, c_range, s_range); % 3 x 1

    xyz_dcm = (repmat(ipp * 1e-3, [1 Nx*Ny*Nz]) +  R_rcs2lph * scaling_dicom * cat(2, r(:), c(:), s(:)).');
    x_dcm   = reshape(xyz_dcm(1,:), [Nx Ny]); % [m]
    y_dcm   = reshape(xyz_dcm(2,:), [Nx Ny]); % [m]
    z_dcm   = reshape(xyz_dcm(3,:), [Nx Ny]); % [m]

    %% Display images in PCS
    FontSize = 14;

    N = 1000;
    t = (0:N-1).' * (2 * pi / N);
    r = 100 * 1e-3;              % [mm] * [m/1e3mm] => [m]
    x_circle_offset = 0;         % [mm] * [m/1e3mm] => [m]
    y_circle_offset = 15 * 1e-3; % [mm] * [m/1e3mm] => [m]
    z_circle_offset = 0;         % [mm] * [m/1e3mm] => [m]

    if strcmp(image_ori, 'sag')
        az = 90;
        el = 0;
        x_circle = zeros(N,1) + x_circle_offset;
        y_circle = r * cos(t) + y_circle_offset;
        z_circle = r * sin(t) + z_circle_offset;
        x_text = 0;
        y_text = 0;
        z_text = max(z(:));
    elseif strcmp(image_ori, 'cor')
        az = 0;
        el = 0;
        x_circle = r * cos(t) + x_circle_offset;
        y_circle = zeros(N,1) + 0;
        z_circle = r * sin(t) + z_circle_offset;
        x_text = 0;
        y_text = 0;
        z_text = max(z(:));
    elseif strcmp(image_ori, 'tra')
        az = 0;
        el = -90;
        x_circle = r * cos(t) + x_circle_offset;
        y_circle = r * sin(t) + y_circle_offset;
        z_circle = zeros(N,1) + z_circle_offset;
        x_text = 0;
        y_text = min(y(:));
        z_text = 0;
    end

    %----------------------------------------------------------------------
    % MATLAB recon
    %----------------------------------------------------------------------
    figure('Color', 'w', 'Position', [1 1 1600 823]);
    color_order = get(gca, 'colororder');
    subplot(1,2,1); hold on;
    surf(x*1e3, y*1e3, z*1e3, abs(im), 'EdgeColor', 'none'); axis image; colormap(gray);
    plot3(x_circle*1e3, y_circle*1e3, z_circle*1e3, '--', 'LineWidth', 2, 'Color', [192 0 0]/255);
    set(gca, 'Box', 'On', 'TickLabelInterpreter', 'latex', 'FontSize', FontSize);
    xlabel('x [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    ylabel('y [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    zlabel('z [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    title(sprintf('ISMRMRD: %s, %s, %g deg', patient_position, image_ori, dRotAngle*180/pi), 'FontSize', FontSize);
    xlim([-230 230]);
    ylim([-230 230]);
    zlim([-230 230]);
    view(az,el);
    text(x_text*1e3, y_text*1e3, z_text*1e3, {'LPH DICOM coordinate system (=PCS)', ...
        sprintf('min(x) = %7.3f, max(x) = %7.3f', min(x(:))*1e3, max(x(:))*1e3), ...
        sprintf('min(y) = %7.3f, max(y) = %7.3f', min(y(:))*1e3, max(y(:))*1e3), ...
        sprintf('min(z) = %7.3f, max(z) = %7.3f', min(z(:))*1e3, max(z(:))*1e3), ...
        }, 'Color', color_order(3,:), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
    subplot(1,2,2); hold on;
    surf(x*1e3, y*1e3, z*1e3, abs(im), 'EdgeColor', 'none'); axis image; colormap(gray);
    plot3(0, 0, 0, '.', 'Color', 'r', 'MarkerSize', 10);
    set(gca, 'Box', 'On', 'TickLabelInterpreter', 'latex', 'FontSize', FontSize);
    xlabel('x [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    ylabel('y [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    zlabel('z [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    title(sprintf('ISMRMRD: %s, %s, %g deg', patient_position, image_ori, dRotAngle*180/pi), 'FontSize', FontSize);
    xlim([-25 25]);
    ylim([-25 25]);
    zlim([-25 25]);
    view(az,el);
    export_fig(fullfile(output_directory, sprintf('%s_ismrmrd_LPH', prefix)), '-m1', '-tif');
    close;

    %----------------------------------------------------------------------
    % DICOM
    %----------------------------------------------------------------------
    figure('Color', 'w', 'Position', [1 1 1600 823]);
    color_order = get(gca, 'colororder');
    subplot(1,2,1); hold on;
    surf(x_dcm*1e3, y_dcm*1e3, z_dcm*1e3, im_dcm, 'EdgeColor', 'none'); axis image; colormap(gray);
    plot3(x_circle*1e3, y_circle*1e3, z_circle*1e3, '--', 'LineWidth', 2, 'Color', [192 0 0]/255);
    set(gca, 'Box', 'On', 'TickLabelInterpreter', 'latex', 'FontSize', FontSize);
    xlabel('x [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    ylabel('y [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    zlabel('z [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    title(sprintf('DICOM: %s, %s, %g deg', patient_position, image_ori, dRotAngle*180/pi), 'FontSize', FontSize);
    xlim([-230 230]);
    ylim([-230 230]);
    zlim([-230 230]);
    view(az,el);
    text(x_text*1e3, y_text*1e3, z_text*1e3, {'LPH DICOM coordinate system (=PCS)', ...
        sprintf('min(x) = %7.3f, max(x) = %7.3f', min(x_dcm(:))*1e3, max(x_dcm(:))*1e3), ...
        sprintf('min(y) = %7.3f, max(y) = %7.3f', min(y_dcm(:))*1e3, max(y_dcm(:))*1e3), ...
        sprintf('min(z) = %7.3f, max(z) = %7.3f', min(z_dcm(:))*1e3, max(z_dcm(:))*1e3), ...
        }, 'Color', color_order(3,:), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
    subplot(1,2,2); hold on;
    surf(x_dcm*1e3, y_dcm*1e3, z_dcm*1e3, im_dcm, 'EdgeColor', 'none'); axis image; colormap(gray);
    plot3(0, 0, 0, '.', 'Color', 'r', 'MarkerSize', 10);
    set(gca, 'Box', 'On', 'TickLabelInterpreter', 'latex', 'FontSize', FontSize);
    xlabel('x [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    ylabel('y [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    zlabel('z [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    title(sprintf('DICOM: %s, %s, %g deg', patient_position, image_ori, dRotAngle*180/pi), 'FontSize', FontSize);
    xlim([-25 25]);
    ylim([-25 25]);
    zlim([-25 25]);
    view(az,el);
    export_fig(fullfile(output_directory, sprintf('%s_dicom_LPH', prefix)), '-m1', '-tif');
    close;
end
