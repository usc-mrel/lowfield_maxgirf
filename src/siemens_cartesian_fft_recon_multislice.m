function [im_multislice,header,r_dcs_multislice] = siemens_cartesian_fft_recon_multislice(ismrmrd_noise_fullpath, ismrmrd_data_fullpath, siemens_dat_fullpath, user_opts)
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/16/2022, Last modified: 02/05/2022

%% Read k-space data (ISMRMRD format)
start_time = tic;
tic; fprintf('Reading an ISMRMRD file: %s... ', ismrmrd_data_fullpath);
if exist(ismrmrd_data_fullpath, 'file')
    dset = ismrmrd.Dataset(ismrmrd_data_fullpath, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , ismrmrd_data_fullpath);
end

%% Get imaging parameters from the XML header
header = ismrmrd.xml.deserialize(dset.readxml);

%--------------------------------------------------------------------------
% measurement information
%--------------------------------------------------------------------------
patient_position = header.measurementInformation.patientPosition;

%--------------------------------------------------------------------------
% Encoding
%--------------------------------------------------------------------------
encoded_fov(1) = header.encoding.encodedSpace.fieldOfView_mm.x; % RO
encoded_fov(2) = header.encoding.encodedSpace.fieldOfView_mm.y; % PE
encoded_fov(3) = header.encoding.encodedSpace.fieldOfView_mm.z; % SL

recon_fov(1) = header.encoding.reconSpace.fieldOfView_mm.x; % RO
recon_fov(2) = header.encoding.reconSpace.fieldOfView_mm.y; % PE
recon_fov(3) = header.encoding.reconSpace.fieldOfView_mm.z; % SL

Nkx = header.encoding.encodedSpace.matrixSize.x; % number of readout samples in k-space
Nky = header.encoding.encodedSpace.matrixSize.y; % number of phase encodes in k-space
Nkz = header.encoding.encodedSpace.matrixSize.z; % number of slice encodes in k-space
Nx  = header.encoding.reconSpace.matrixSize.x;   % number of samples in image-space (RO)
Ny  = header.encoding.reconSpace.matrixSize.y;   % number of samples in image-space (PE)
Nz  = header.encoding.reconSpace.matrixSize.z;   % number of samples in image-space (SL)
Nc  = header.acquisitionSystemInformation.receiverChannels;

encoded_resolution = encoded_fov ./ [Nkx Nky Nkz]; % [mm]
recon_resolution = recon_fov ./ [Nx Ny Nz]; % [mm]

if isfield(user_opts, 'remove_oversampling')
    remove_oversampling = user_opts.remove_oversampling;
else
    remove_oversampling = 1;
end

if remove_oversampling
    osf = encoded_fov(1) / recon_fov(1); % oversampling factor in the x direction
else
    osf = 1;
end

zpad_matrix_size(1) = floor(encoded_fov(1) / osf / recon_resolution(1));
zpad_matrix_size(2) = round(encoded_fov(2) / recon_resolution(2));
zpad_matrix_size(3) = round(encoded_fov(3) / recon_resolution(3));

%--------------------------------------------------------------------------
% Calculate reconstruction parameters
%--------------------------------------------------------------------------
Nk = Nkx / osf; % number of readout samples
N1 = Nk;
N2 = Nky;
N3 = Nkz;
N = N1 * N2 * N3;

%% Parse the ISMRMRD header
tic; fprintf('Parsing the ISMRMRD header... ');
raw_data = dset.readAcquisition(); % read all the acquisitions
fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

%--------------------------------------------------------------------------
% ISMRMRD header
%--------------------------------------------------------------------------
% uint16_t version;                                    /**< First unsigned int indicates the version */
% uint64_t flags;                                      /**< bit field with flags */
% uint32_t measurement_uid;                            /**< Unique ID for the measurement */
% uint32_t scan_counter;                               /**< Current acquisition number in the measurement */
% uint32_t acquisition_time_stamp;                     /**< Acquisition clock */
% uint32_t physiology_time_stamp[ISMRMRD_PHYS_STAMPS]; /**< Physiology time stamps, e.g. ecg, breating, etc. */
% uint16_t number_of_samples;                          /**< Number of samples acquired */
% uint16_t available_channels;                         /**< Available coils */
% uint16_t active_channels;                            /**< Active coils on current acquisiton */
% uint64_t channel_mask[ISMRMRD_CHANNEL_MASKS];        /**< Mask to indicate which channels are active. Support for 1024 channels */
% uint16_t discard_pre;                                /**< Samples to be discarded at the beginning of  acquisition */
% uint16_t discard_post;                               /**< Samples to be discarded at the end of acquisition */
% uint16_t center_sample;                              /**< Sample at the center of k-space */
% uint16_t encoding_space_ref;                         /**< Reference to an encoding space, typically only one per acquisition */
% uint16_t trajectory_dimensions;                      /**< Indicates the dimensionality of the trajectory vector (0 means no trajectory) */
% float sample_time_us;                                /**< Time between samples in micro seconds, sampling BW */
% float position[3];                                   /**< Three-dimensional spatial offsets from isocenter */
% float read_dir[3];                                   /**< Directional cosines of the readout/frequency encoding */
% float phase_dir[3];                                  /**< Directional cosines of the phase */
% float slice_dir[3];                                  /**< Directional cosines of the slice direction */
% float patient_table_position[3];                     /**< Patient table off-center */
% ISMRMRD_EncodingCounters idx;                        /**< Encoding loop counters, see above */
% int32_t user_int[ISMRMRD_USER_INTS];                 /**< Free user parameters */
% float user_float[ISMRMRD_USER_FLOATS];               /**< Free user parameters */
%--------------------------------------------------------------------------
% Where EncodingCounters are defined as:
% uint16_t kspace_encode_step_1;    /**< e.g. phase encoding line number */
% uint16_t kspace_encode_step_2;    /**< e.g. partition encoding number */
% uint16_t average;                 /**< e.g. signal average number */
% uint16_t slice;                   /**< e.g. imaging slice number */
% uint16_t contrast;                /**< e.g. echo number in multi-echo */
% uint16_t phase;                   /**< e.g. cardiac phase number */
% uint16_t repetition;              /**< e.g. dynamic number for dynamic scanning */
% uint16_t set;                     /**< e.g. flow encoding set */
% uint16_t segment;                 /**< e.g. segment number for segmented acquisition */
% uint16_t user[ISMRMRD_USER_INTS]; /**< Free user parameters */
%--------------------------------------------------------------------------
number_of_samples  = double(max(raw_data.head.number_of_samples));
discard_pre        = double(max(raw_data.head.discard_pre));
discard_post       = double(max(raw_data.head.discard_post));
center_sample      = double(max(raw_data.head.center_sample));
nr_channels        = double(max(raw_data.head.active_channels));
nr_phase_encodings = double(max(raw_data.head.idx.kspace_encode_step_1)) + 1; % nr_interleaves for spiral imaging
nr_slice_encodings = double(max(raw_data.head.idx.kspace_encode_step_2)) + 1;
nr_averages        = double(max(raw_data.head.idx.average)) + 1;
nr_slices          = double(max(raw_data.head.idx.slice)) + 1;
nr_contrasts       = double(max(raw_data.head.idx.contrast)) + 1;
nr_phases          = double(max(raw_data.head.idx.phase)) + 1;
nr_repetitions     = double(max(raw_data.head.idx.repetition)) + 1;
nr_sets            = double(max(raw_data.head.idx.set)) + 1;
nr_segments        = double(max(raw_data.head.idx.segment)) + 1;
nr_samples         = number_of_samples - discard_pre - discard_post;

%--------------------------------------------------------------------------
% Get the dimensionality of the trajectory vector (0 means no trajectory)
%--------------------------------------------------------------------------
trajectory_dimensions = double(max(raw_data.head.trajectory_dimensions));

%--------------------------------------------------------------------------
% Get the dwell time in [sec]
%--------------------------------------------------------------------------
dt = double(max(raw_data.head.sample_time_us)) * 1e-6; % [usec] * [sec/1e-6 usec] => [sec]

%--------------------------------------------------------------------------
% Calculate the readout duration [sec]
%--------------------------------------------------------------------------
T = number_of_samples * dt; % readout duration [sec]

%% Display ISMRMRD header
fprintf('========================= ISMRMRD header ========================\n');
fprintf('encoded_fov        = %8.4f %8.4f %8.4f\n', encoded_fov(1), encoded_fov(2), encoded_fov(3));
fprintf('Nkx Nky Nkz        = %d      %d        %d\n', Nkx, Nky, Nkz);
fprintf('encoded_resolution = %8.4f %8.4f %8.4f\n', encoded_resolution(1), encoded_resolution(2), encoded_resolution(3));
fprintf('-----------------------------------------------------------------\n');
fprintf('recon_fov          = %8.4f %8.4f %8.4f\n', recon_fov(1), recon_fov(2), recon_fov(3));
fprintf('Nx Ny Nz           = %d      %d        %d\n', Nx, Ny, Nz);
fprintf('recon_resolution   = %8.4f %8.4f %8.4f\n', recon_resolution(1), recon_resolution(2), recon_resolution(3));
fprintf('-----------------------------------------------------------------\n');
fprintf('trajectory         = %s\n', header.encoding.trajectory);
fprintf('number_of_samples  = %d\n', number_of_samples);
fprintf('discard_pre        = %d\n', discard_pre);
fprintf('discard_post       = %d\n', discard_post);
fprintf('center_sample      = %d\n', center_sample);
fprintf('nr_channels        = %d\n', nr_channels);
fprintf('nr_phase_encodings = %d\n', nr_phase_encodings);
fprintf('nr_slice_encodings = %d\n', nr_slice_encodings);
fprintf('nr_averages        = %d\n', nr_averages);
fprintf('nr_slices          = %d\n', nr_slices);
fprintf('nr_contrasts       = %d\n', nr_contrasts);
fprintf('nr_phases          = %d\n', nr_phases);
fprintf('nr_repetitions     = %d\n', nr_repetitions);
fprintf('nr_sets            = %d\n', nr_sets);
fprintf('nr_segments        = %d\n', nr_segments);
fprintf('dt                 = %5.2f [usec]\n', dt * 1e6);
fprintf('readout duration   = %5.2f [msec]\n', T * 1e3);
fprintf('=================================================================\n');

%% Calculate the receiver noise matrix
[Psi,inv_L] = calculate_noise_decorrelation_matrix(ismrmrd_noise_fullpath);

%% Read a Siemens .dat file
if exist(siemens_dat_fullpath, 'file')
    fprintf('Reading a Siemens .dat file: %s\n', siemens_dat_fullpath);
    twix = mapVBVD(siemens_dat_fullpath);
    if length(twix) > 1
        twix = twix{end};
    end
end

%% Perform FFT reconstruction per slice
nr_recons = nr_slices * nr_contrasts * nr_phases * nr_repetitions * nr_sets;
im_multislice = complex(zeros(N1, N2, nr_slices, 'double'));
r_dcs_multislice = zeros(N, 3, nr_slices, 'double');

for idx = 1:nr_recons
    %% Get information about the current slice
    [slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr] = ind2sub([nr_slices nr_contrasts nr_phases nr_repetitions nr_sets], idx);
    fprintf('(%2d/%2d): Reconstructing slice (slice = %2d, contrast = %2d, phase = %2d, repetition = %2d, set = %2d)\n', idx, nr_recons, slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr);

    %% Calculate the actual slice number for Siemens interleaved multislice imaging
    if nr_slices > 1 % multi-slice
        if mod(nr_slices,2) == 0 % even
            offset1 = 0;
            offset2 = 1;
        else % odd
            offset1 = 1;
            offset2 = 0;
        end
        if slice_nr <= ceil(nr_slices / 2)
            actual_slice_nr = 2 * slice_nr - offset1;
        else
            actual_slice_nr = 2 * (slice_nr - ceil(nr_slices/2)) - offset2;
        end
    else
        actual_slice_nr = slice_nr;
    end

    %% Get a list of profiles for all segments
    profile_list = find((raw_data.head.idx.slice      == (slice_nr - 1))      & ...
                        (raw_data.head.idx.contrast   == (contrast_nr - 1))   & ...
                        (raw_data.head.idx.phase      == (phase_nr - 1))      & ...
                        (raw_data.head.idx.repetition == (repetition_nr - 1)) & ...
                        (raw_data.head.idx.set        == (set_nr - 1)));

    %% Get a slice offset in PCS from ISMRMRD format
    sag_offset_ismrmrd = double(raw_data.head.position(1,profile_list(1))); % [mm]
    cor_offset_ismrmrd = double(raw_data.head.position(2,profile_list(1))); % [mm]
    tra_offset_ismrmrd = double(raw_data.head.position(3,profile_list(1))); % [mm]

    %% Get a rotation matrix from GCS to PCS (ISMRMRD format)
    phase_dir = double(raw_data.head.phase_dir(:,profile_list(1)));
    read_dir  = double(raw_data.head.read_dir(:,profile_list(1)));
    slice_dir = double(raw_data.head.slice_dir(:,profile_list(1)));
    R_gcs2pcs_ismrmrd = [phase_dir read_dir slice_dir];

    %% Get information from Siemens TWIX format
    if exist(siemens_dat_fullpath, 'file')
        %% Get a slice normal vector from Siemens TWIX format
        %------------------------------------------------------------------
        % dNormalSag: Sagittal component of a slice normal vector (in PCS)
        %------------------------------------------------------------------
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal, 'dSag')
            dNormalSag = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal.dSag;
        else
            dNormalSag = 0;
        end

        %------------------------------------------------------------------
        % dNormalCor: Coronal component of a slice normal vector (in PCS)
        %------------------------------------------------------------------
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal, 'dCor')
            dNormalCor = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal.dCor;
        else
            dNormalCor = 0;
        end

        %------------------------------------------------------------------
        % dNormalTra: Transverse component of a slice normal vector (in PCS)
        %------------------------------------------------------------------
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal, 'dTra')
            dNormalTra = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal.dTra;
        else
            dNormalTra = 0;
        end

        %------------------------------------------------------------------
        % dRotAngle: Slice rotation angle ("swap Fre/Pha")
        %------------------------------------------------------------------
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}, 'dInPlaneRot')
            dRotAngle = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.dInPlaneRot; % [rad]
        else
            dRotAngle = 0; % [rad]
        end

        %% Get a slice offset in PCS from Siemens TWIX format
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}, 'sPosition')
            if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition, 'dSag')
                sag_offset_twix = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition.dSag; % [mm]
            else
                sag_offset_twix = 0; % [mm]
            end
            if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition, 'dCor')
                cor_offset_twix = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition.dCor; % [mm]
            else
                cor_offset_twix = 0; % [mm]
            end
            if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition, 'dTra')
                tra_offset_twix = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition.dTra; % [mm]
            else
                tra_offset_twix = 0; % [mm]
            end
        else
            sag_offset_twix = 0; % [mm]
            cor_offset_twix = 0; % [mm]
            tra_offset_twix = 0; % [mm]
        end
        pcs_offset_twix = [sag_offset_twix; cor_offset_twix; tra_offset_twix] * 1e-3; % [mm] * [m/1e3mm] => [m]
    end

    %% Use a slice offset in PCS from Siemens TWIX format
    pcs_offset = [sag_offset_twix; cor_offset_twix; tra_offset_twix] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Calculate spatial coordinates in DCS [m]
    %----------------------------------------------------------------------
    % Calculate a transformation matrix from RCS to GCS [r,c,s] <=> [PE,RO,SL]
    %----------------------------------------------------------------------
    R_rcs2gcs = [0    1    0 ; % [PE]   [0 1 0] * [r]
                 1    0    0 ; % [RO] = [1 0 0] * [c]
                 0    0    1]; % [SL]   [0 0 1] * [s]

    %----------------------------------------------------------------------
    % Calculate a rotation matrix from GCS to PCS
    %----------------------------------------------------------------------
    if exist(siemens_dat_fullpath, 'file')
        [R_gcs2pcs,phase_sign,read_sign] = siemens_calculate_matrix_gcs_to_pcs(dNormalSag, dNormalCor, dNormalTra, dRotAngle);
    else
        phase_sign = user_opts.phase_sign;
        read_sign = user_opts.read_sign;
        R_gcs2pcs = [phase_sign * phase_dir read_sign * read_dir slice_dir];
    end

    %----------------------------------------------------------------------
    % Calculate a rotation matrix from PCS to DCS
    %----------------------------------------------------------------------
    R_pcs2dcs = siemens_calculate_matrix_pcs_to_dcs(patient_position);

    %----------------------------------------------------------------------
    % Calculate a rotation matrix from GCS to DCS
    %----------------------------------------------------------------------
    R_gcs2dcs = R_pcs2dcs * R_gcs2pcs;

    %----------------------------------------------------------------------
    % Calculate a scaling matrix
    %----------------------------------------------------------------------
    scaling_matrix = diag(encoded_resolution) * 1e-3; % [mm] * [m/1e3mm] => [m]

    %----------------------------------------------------------------------
    % Calculate a slice offset in DCS [m]
    %----------------------------------------------------------------------
    dcs_offset = R_pcs2dcs * pcs_offset; % 3 x 1

    %----------------------------------------------------------------------
    % Calculate spatial coordinates in RCS [m]
    %----------------------------------------------------------------------
    [I1,I2,I3] = ndgrid((1:N1).', (1:N2).', (1:N3).');
    r_rcs = (scaling_matrix * cat(2, I1(:) - (floor(N1/2) + 1), I2(:) - (floor(N2/2) + 1), I3(:) - (floor(N3/2) + 1)).').'; % N x 3

    %----------------------------------------------------------------------
    % Calculate spatial coordinates in GCS [m] [PE,RO,SL]
    %----------------------------------------------------------------------
    r_gcs = (R_rcs2gcs * r_rcs.').'; % N x 3

    %----------------------------------------------------------------------
    % Calculate spatial coordinates in DCS [m]
    %----------------------------------------------------------------------
    r_dcs = (repmat(dcs_offset, [1 N]) + R_pcs2dcs * R_gcs2pcs * r_gcs.').'; % N x 3
    r_dcs_multislice(:,:,idx) = r_dcs;

    %% Display slice information
    fprintf('======================= SLICE INFORMATION ========================\n');
    fprintf('slice_nr = %d, actual_slice_nr = %d\n', slice_nr, actual_slice_nr);
    fprintf('phase_sign = %+g, read_sign = %+g\n', phase_sign, read_sign);
    fprintf('---------------------- From Siemens TWIX format ------------------\n');
    fprintf('                   [sag]   %10.5f [mm]\n', sag_offset_twix);
    fprintf('slice offset(PCS): [cor] = %10.5f [mm]\n', cor_offset_twix);
    fprintf('                   [tra]   %10.5f [mm]\n', tra_offset_twix);
    fprintf('---------------------- From ISMRMRD format -----------------------\n');
    fprintf('                   [sag]   %10.5f [mm]\n', sag_offset_ismrmrd);
    fprintf('slice offset(PCS): [cor] = %10.5f [mm]\n', cor_offset_ismrmrd);
    fprintf('                   [tra]   %10.5f [mm]\n', tra_offset_ismrmrd);
    fprintf('---------------------- From Siemens TWIX format ------------------\n');
    fprintf('                   [sag]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2pcs(1,1), R_gcs2pcs(1,2), R_gcs2pcs(1,3));
    fprintf('R_gcs2pcs        : [cor] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2pcs(2,1), R_gcs2pcs(2,2), R_gcs2pcs(2,3));
    fprintf('                   [tra]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2pcs(3,1), R_gcs2pcs(3,2), R_gcs2pcs(3,3));
    fprintf('---------------------- From ISMRMRD format (incorrect!)-----------\n');
    fprintf('                   [sag]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2pcs_ismrmrd(1,1), R_gcs2pcs_ismrmrd(1,2), R_gcs2pcs_ismrmrd(1,3));
    fprintf('R_gcs2pcs        : [cor] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2pcs_ismrmrd(2,1), R_gcs2pcs_ismrmrd(2,2), R_gcs2pcs_ismrmrd(2,3));
    fprintf('                   [tra]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2pcs_ismrmrd(3,1), R_gcs2pcs_ismrmrd(3,2), R_gcs2pcs_ismrmrd(3,3));
    fprintf('------------------------------------------------------------------\n');
    fprintf('                   [ x ]   [%10.5f %10.5f %10.5f][sag]\n', R_pcs2dcs(1,1), R_pcs2dcs(1,2), R_pcs2dcs(1,3));
    fprintf('R_pcs2dcs        : [ y ] = [%10.5f %10.5f %10.5f][cor]\n', R_pcs2dcs(2,1), R_pcs2dcs(2,2), R_pcs2dcs(2,3));
    fprintf('                   [ z ]   [%10.5f %10.5f %10.5f][tra]\n', R_pcs2dcs(3,1), R_pcs2dcs(3,2), R_pcs2dcs(3,3));
    fprintf('------------------------------------------------------------------\n');
    fprintf('                   [ x ]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2dcs(1,1), R_gcs2dcs(1,2), R_gcs2dcs(1,3));
    fprintf('R_gcs2dcs        : [ y ] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2dcs(2,1), R_gcs2dcs(2,2), R_gcs2dcs(2,3));
    fprintf('                   [ z ]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2dcs(3,1), R_gcs2dcs(3,2), R_gcs2dcs(3,3));
    fprintf('==================================================================\n');

    %% Get information about k-space encoding
    kspace_encoding_step_1_center = header.encoding.encodingLimits.kspace_encoding_step_1.center;
    kspace_encoding_step_2_center = header.encoding.encodingLimits.kspace_encoding_step_2.center;
    kspace_encoding_step_1_maximum = header.encoding.encodingLimits.kspace_encoding_step_1.maximum;
    kspace_encoding_step_2_maximum = header.encoding.encodingLimits.kspace_encoding_step_2.maximum;

    %% Calculate a pattern showing the number of averages at each k-space location (Nky x Nkz)
    average_pattern = zeros(Nky, Nkz, 'double');
    for segment_nr = 1:nr_segments
        %------------------------------------------------------------------
        % Get a list of profiles per segment
        %------------------------------------------------------------------
        profile_list_segment = find((raw_data.head.idx.slice      == (slice_nr - 1))      & ...
                                    (raw_data.head.idx.contrast   == (contrast_nr - 1))   & ...
                                    (raw_data.head.idx.phase      == (phase_nr - 1))      & ...
                                    (raw_data.head.idx.repetition == (repetition_nr - 1)) & ...
                                    (raw_data.head.idx.set        == (set_nr - 1))        & ...
                                    (raw_data.head.idx.segment    == (segment_nr - 1)));

        %------------------------------------------------------------------
        % Calculate a sampling pattern (Nky x Nkz)
        %------------------------------------------------------------------
        for idx1 = 1:length(profile_list_segment)
            kspace_encode_step_1 = double(raw_data.head.idx.kspace_encode_step_1(profile_list_segment(idx1)));
            kspace_encode_step_2 = double(raw_data.head.idx.kspace_encode_step_2(profile_list_segment(idx1)));
            k2_index = kspace_encode_step_1 - (kspace_encoding_step_1_maximum - kspace_encoding_step_1_center + 1) + floor(Nky/2) + 1;
            k3_index = kspace_encode_step_2 - (kspace_encoding_step_2_maximum - kspace_encoding_step_2_center + 1) + floor(Nkz/2) + 1;
            if k3_index == 0, k3_index = k3_index + 1; end % For 2D imaging
            average_pattern(k2_index,k3_index) = average_pattern(k2_index,k3_index) + 1;
        end
    end

    %% Get k-space data (Nkx x Nky x Nkz x Nc)
    %--------------------------------------------------------------------------
    % Calculate the index range (1st index) of a profile
    %--------------------------------------------------------------------------
    k1_range = ((discard_pre + 1):(number_of_samples - discard_post)).' + (Nkx - nr_samples); % need to confirm this with partial Fourier datasets

    kspace = complex(zeros(Nkx, Nky, Nkz, Nc, 'double'));
    for idx1 = 1:length(profile_list)
        tstart = tic; fprintf('(%2d/%2d): Reading k-space data (%d/%d)... ', idx, nr_recons, idx1, length(profile_list));
        %------------------------------------------------------------------
        % Calculate the (2nd,3rd) matrix index of a profile
        %------------------------------------------------------------------
        kspace_encode_step_1 = double(raw_data.head.idx.kspace_encode_step_1(profile_list(idx1)));
        kspace_encode_step_2 = double(raw_data.head.idx.kspace_encode_step_2(profile_list(idx1)));
        k2_index = kspace_encode_step_1 - (kspace_encoding_step_1_maximum - kspace_encoding_step_1_center + 1) + floor(Nky/2) + 1;
        k3_index = kspace_encode_step_2 - (kspace_encoding_step_2_maximum - kspace_encoding_step_2_center + 1) + floor(Nkz/2) + 1;
        if k3_index == 0, k3_index = k3_index + 1; end % For 2D imaging

        %------------------------------------------------------------------
        % Prewhiten k-space data
        %------------------------------------------------------------------
        profile = raw_data.data{profile_list(idx1)}; % number_of_samples x nr_channels
        profile = (inv_L * profile.').';

        %------------------------------------------------------------------
        % Calculate the average of k-space data
        %------------------------------------------------------------------
        profile = profile / average_pattern(k2_index,k3_index);

        %------------------------------------------------------------------
        % Accumulate k-space
        %------------------------------------------------------------------
        kspace(k1_range,k2_index,k3_index,:) = kspace(k1_range,k2_index,k3_index,:) + reshape(profile, [nr_samples 1 1 Nc]);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %% Flip k-space
    tstart = tic; fprintf('(%2d/%2d): Flipping k-space... ', idx, nr_recons);
    if read_sign == -1
        kspace = flip(kspace,1);
    end
    if phase_sign == -1
        kspace = flip(kspace,2);
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Remove readout oversampling: (Nkx x Nky x Nkz x Nc) => (Nk x Nky x Nkz x Nc)
    %----------------------------------------------------------------------
    % FFT along the readout direction (k-space <=> image-space)
    %----------------------------------------------------------------------
    tstart = tic; fprintf('(%2d/%2d): Removing readout oversampling... ', idx, nr_recons);
    projection = 1 / sqrt(Nkx) * fftshift(fft(ifftshift(kspace, 1), [], 1), 1);
    idx1_range = (-floor(Nk/2):ceil(Nk/2)-1).' + floor(Nkx/2) + 1;
    kspace = sqrt(Nk / osf) * fftshift(ifft(ifftshift(projection(idx1_range,:,:,:), 1), [], 1), 1);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    clear projection;

    %% Calculate coil sensitivity maps
    tstart = tic; fprintf('(%2d/%2d): Calculating coil sensitivity maps with Walsh method... ', idx, nr_recons);
    %----------------------------------------------------------------------
    % Calculate the calibration region of k-space
    %----------------------------------------------------------------------
    cal_shape = [32 32];
    cal_data = crop(reshape(kspace, [N1 N2 Nc]), [cal_shape Nc]);
    cal_data = bsxfun(@times, cal_data, hamming(cal_shape(1)) * hamming(cal_shape(2)).');

    %----------------------------------------------------------------------
    % Calculate coil sensitivity maps
    %----------------------------------------------------------------------
    cal_im = zpad(cal_data, [N1 N2 Nc]);
    for dim = 1:2
        cal_im = 1 / sqrt(size(cal_im,dim)) * fftshift(fft(ifftshift(cal_im, dim), [], dim), dim);
    end
    csm = ismrm_estimate_csm_walsh(cal_im);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Perform optimal coil combination
    imc = reshape(kspace, [N1 N2 Nc]);
    for dim = 1:2
        imc = 1 / sqrt(size(imc,dim)) * fftshift(fft(ifftshift(imc, dim), [], dim), dim);
    end
    im_multislice(:,:,idx) = sum(bsxfun(@times, conj(csm), imc), 3);
end
end