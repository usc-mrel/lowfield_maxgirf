function [im,header,x,y,z,imc,csm] = cartesian_recon(noise_fullpath, cartesian_fullpath, dat_fullpath, user_opts_cartesian)
% Inputs
%   noise_fullpath       full file path to noise only h5 file
%   cartesian_fullpath   full file path to cartesian h5 file
%   dat_fullpath         full file path to Siemens dat file

%% Define constants
gamma = 2.67522212e+8;    % gyromagnetic ratio for 1H [rad/sec/T]

%% Read an ismrmrd file
start_time = tic;
tic; fprintf('Reading an ismrmrd file: %s... ', cartesian_fullpath);
if exist(cartesian_fullpath, 'file')
    dset = ismrmrd.Dataset(cartesian_fullpath, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , cartesian_fullpath);
end

%% Get imaging parameters from the XML header
header = ismrmrd.xml.deserialize(dset.readxml);

%--------------------------------------------------------------------------
% measurement information
%--------------------------------------------------------------------------
patient_position = header.measurementInformation.patientPosition;

%--------------------------------------------------------------------------
% Sequence parameters
%--------------------------------------------------------------------------
TR         = header.sequenceParameters.TR * 1e-3;     % [msec] * [sec/1e3msec] => [sec]
TE         = header.sequenceParameters.TE * 1e-3;     % [msec] * [sec/1e3msec] => [sec]
flip_angle = header.sequenceParameters.flipAngle_deg; % [degrees]

%--------------------------------------------------------------------------
% Experimental conditions
%--------------------------------------------------------------------------
B0 = header.experimentalConditions.H1resonanceFrequency_Hz * (2 * pi / gamma); % [Hz] * [2pi rad/cycle] / [rad/sec/T] => [T]

%--------------------------------------------------------------------------
% Encoding
%--------------------------------------------------------------------------
field_of_view_mm(1) = header.encoding.reconSpace.fieldOfView_mm.x; % RO
field_of_view_mm(2) = header.encoding.reconSpace.fieldOfView_mm.y; % PE
field_of_view_mm(3) = header.encoding.reconSpace.fieldOfView_mm.z; % SS

Nkx = header.encoding.encodedSpace.matrixSize.x; % number of samples in k-space (RO,row)
Nky = header.encoding.encodedSpace.matrixSize.y; % number of samples in k-space (PE,col)
Nkz = header.encoding.encodedSpace.matrixSize.z; % number of samples in k-space (SS,slice)
Nx  = header.encoding.reconSpace.matrixSize.x;   % number of samples in image-space (RO,row)
Ny  = header.encoding.reconSpace.matrixSize.y;   % number of samples in image-space (PE,col)
Nz  = header.encoding.reconSpace.matrixSize.z;   % number of samples in image-space (SS,slice)

%--------------------------------------------------------------------------
% Zeropad in kspace
%--------------------------------------------------------------------------
zpad_factor = user_opts_cartesian.zpad_factor;
Nx = Nx * zpad_factor;
Ny = Ny * zpad_factor;
N  = Nx * Ny * Nz; % total number of voxels in image-space

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
nr_samples         = number_of_samples - discard_pre - discard_post;
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

%--------------------------------------------------------------------------
% Get the dimensionality of the trajectory vector (0 means no trajectory)
%--------------------------------------------------------------------------
trajectory_dimensions = double(max(raw_data.head.trajectory_dimensions));

%--------------------------------------------------------------------------
% Get the dwell time in [sec]
%--------------------------------------------------------------------------
dt = double(max(raw_data.head.sample_time_us)) * 1e-6; % [usec] * [sec/1e-6 usec] => [sec]

fprintf('trajectory         = %s\n', header.encoding.trajectory);
fprintf('number_of_samples  = %d\n', number_of_samples);
fprintf('discard_pre        = %d\n', discard_pre);
fprintf('discard_post       = %d\n', discard_post);
fprintf('nr_samples         = %d\n', nr_samples);
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
fprintf('dt                 = %g [usec]\n', dt * 1e6);
fprintf('readout duration   = %g [msec]\n', nr_samples * dt * 1e3);

%% Get all nominal k-space trajectories and k-space data
tic; fprintf('Sorting out all k-space data... ');

data = complex(zeros(nr_samples, nr_phase_encodings, nr_slice_encodings, nr_averages, nr_slices, nr_contrasts, nr_phases, nr_repetitions, nr_sets, nr_segments, nr_channels, 'double'));

nr_profiles = length(raw_data.data);
for idx = 1:nr_profiles
    %----------------------------------------------------------------------
    % Get information about the current profile
    %----------------------------------------------------------------------
    phase_encoding_nr = double(raw_data.head.idx.kspace_encode_step_1(idx)) + 1;
    slice_encoding_nr = double(raw_data.head.idx.kspace_encode_step_2(idx)) + 1;
    average_nr        = double(raw_data.head.idx.average(idx)) + 1;
    slice_nr          = double(raw_data.head.idx.slice(idx)) + 1;
    contrast_nr       = double(raw_data.head.idx.contrast(idx)) + 1;
    phase_nr          = double(raw_data.head.idx.phase(idx)) + 1;
    repetition_nr     = double(raw_data.head.idx.repetition(idx)) + 1;
    set_nr            = double(raw_data.head.idx.set(idx)) + 1;
    segment_nr        = double(raw_data.head.idx.segment(idx)) + 1;

    %----------------------------------------------------------------------
    % Set the index range of a profile
    %----------------------------------------------------------------------
    number_of_samples = double(raw_data.head.number_of_samples(idx));
    discard_pre       = double(raw_data.head.discard_pre(idx));
    discard_post      = double(raw_data.head.discard_post(idx));
    start_index       = discard_pre + 1;
    end_index         = number_of_samples - discard_post;
    index_range       = (start_index:end_index).';

    %----------------------------------------------------------------------
    % Get k-space data
    %----------------------------------------------------------------------
    full_profile = raw_data.data{idx}; % number_of_samples x nr_channels
    data(index_range, phase_encoding_nr, slice_encoding_nr, average_nr, slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr, segment_nr, :) = full_profile(index_range,:); % nr_samples x nr_channels
end
fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

%% Process noise only ismrmrd data
if ~isempty(noise_fullpath)
    [Psi,inv_L] = process_noise_ismrmrd_data(noise_fullpath);
else
    Psi = eye(nr_channels);
    inv_L = eye(nr_channels);
end

%% Read a Siemens dat file
if ~isempty(dat_fullpath)
    fprintf('Reading a Siemens .dat file: %s\n', dat_fullpath);
    twix = mapVBVD(dat_fullpath);
    if length(twix) > 1
        twix = twix{end};
    end
else
    fprintf('\nFile %s does not exist. Skip reading this file.\n' , noise_fullpath);
end

%% Define parameters used in reconstruction
Nc = nr_channels; % number of coils

%% Process a particular dataset at a time
im  = complex(zeros(Nx, Ny, Nz, nr_averages, nr_slices, nr_contrasts, nr_phases, nr_repetitions, nr_sets, nr_segments, 'double'));
imc = complex(zeros(Nx, Ny, Nz, nr_averages, nr_slices, nr_contrasts, nr_phases, nr_repetitions, nr_sets, nr_segments, nr_channels, 'double'));
csm = complex(zeros(Nx, Ny, Nz, nr_averages, nr_slices, nr_contrasts, nr_phases, nr_repetitions, nr_sets, nr_segments, nr_channels, 'double'));
x = zeros(Nx, Ny, Nz, nr_averages, nr_slices, nr_contrasts, nr_phases, nr_repetitions, nr_sets, nr_segments, 'double');
y = zeros(Nx, Ny, Nz, nr_averages, nr_slices, nr_contrasts, nr_phases, nr_repetitions, nr_sets, nr_segments, 'double');
z = zeros(Nx, Ny, Nz, nr_averages, nr_slices, nr_contrasts, nr_phases, nr_repetitions, nr_sets, nr_segments, 'double');
nr_recons = nr_averages * nr_slices * nr_contrasts * nr_phases * nr_repetitions * nr_sets * nr_segments;

for idx = 1:nr_recons
    %% Get information about the current slice
    [average_nr, slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr, segment_nr] = ind2sub([nr_averages nr_slices nr_contrasts nr_phases nr_repetitions nr_sets nr_segments], idx);
    tic; fprintf('(%2d/%2d): Reconstructing a slice (average = %d, slice = %d, contrast = %d, phase = %d, repetition = %d, set = %d, segment = %d)\n', idx, nr_recons, average_nr, slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr, segment_nr);

    %% Calculate the actual slice number
    %----------------------------------------------------------------------
    % my slice number    : 1 2 3 4 5  6 7 8 9 10
    % actual slice number: 2 4 6 8 10 1 3 5 7 9
    % (bottom to top)
    %----------------------------------------------------------------------
    % ISMRMRD (slice_nr):
    %   1,   2, 3,  4,  5,   6,     7,     8,    9,   10,   11
    % -70, -35, 0, 35, 70, 105, -52.5, -17.5, 17.5, 52.5, 87.5
    % twix (actual_slice_nr):
    % -70, -52.5, -35, -17.5, 0, 17.5, 35, 52.5, 70, 87.5, 105
    %   1,     2,   3,     4, 5,    6,  7,    8,  9,   10,  11
    % slice_nr => actual_slice_nr:
    % 1  => 1
    % 2  => 3
    % 3  => 5
    % 4  => 7
    % 5  => 9
    % 6  => 11
    % 7  => 2
    % 8  => 4
    % 9  => 6
    % 10 => 8
    % 11 => 10
    %----------------------------------------------------------------------
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

    %% Get a slice normal vector from Siemens raw data
    %----------------------------------------------------------------------
    % dNormalSag: Sagittal component of a slice normal vector (in PCS)
    %----------------------------------------------------------------------
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal, 'dSag')
        dNormalSag = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal.dSag;
    else
        dNormalSag = 0;
    end

    %----------------------------------------------------------------------
    % dNormalCor: Coronal component of a slice normal vector (in PCS)
    %----------------------------------------------------------------------
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal, 'dCor')
        dNormalCor = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal.dCor;
    else
        dNormalCor = 0;
    end

    %----------------------------------------------------------------------
    % dNormalTra: Transverse component of a slice normal vector (in PCS)
    %----------------------------------------------------------------------
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal, 'dTra')
        dNormalTra = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal.dTra;
    else
        dNormalTra = 0;
    end

    %----------------------------------------------------------------------
    % dRotAngle: Slice rotation angle ("swap Fre/Pha")
    %----------------------------------------------------------------------
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}, 'dInPlaneRot')
        dRotAngle = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.dInPlaneRot; % [rad]
    else
        dRotAngle = 0; % [rad]
    end

    %% Get a slice offset in PCS from Siemens raw data
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}, 'sPosition')
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition, 'dSag')
            sag_offset = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition.dSag; % [mm]
        else
            sag_offset = 0; % [mm]
        end
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition, 'dCor')
            cor_offset = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition.dCor; % [mm]
        else
            cor_offset = 0; % [mm]
        end
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition, 'dTra')
            tra_offset = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition.dTra; % [mm]
        else
            tra_offset = 0; % [mm]
        end
    else
        sag_offset = 0; % [mm]
        cor_offset = 0; % [mm]
        tra_offset = 0; % [mm]
    end
    PCS_offset = [sag_offset; cor_offset; tra_offset] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Calculate spatial coordinates in DCS [mm]
    %----------------------------------------------------------------------
    % Calculate a scaling matrix
    %----------------------------------------------------------------------
    scaling_matrix = diag(field_of_view_mm ./ [Nx Ny Nz] * 1e-3); % [mm] * [m/1e3mm] => [m]

    %----------------------------------------------------------------------
    % Calculate a transformation matrix from RCS to GCS [r,c,s] <=> [PE,RO,SS] 
    %----------------------------------------------------------------------
    rotMatrixRCSToGCS = [0    1    0 ; % [PE]   [0 1 0] * [r]
                         1    0    0 ; % [RO] = [1 0 0] * [c]
                         0    0    1]; % [SS]   [0 0 1] * [s]

    %----------------------------------------------------------------------
    % Calculate a rotation matrix from GCS to PCS from Siemens raw data
    %----------------------------------------------------------------------
    [rotMatrixGCSToPCS,PE_sign,RO_sign] = calcMatrixGCSToPCS(dNormalSag, dNormalCor, dNormalTra, dRotAngle);

    %----------------------------------------------------------------------
    % Calculate a rotation matrix from PCS to DCS
    %----------------------------------------------------------------------
    rotMatrixPCSToDCS = calcMatrixPCSToDCS(patient_position);

    %----------------------------------------------------------------------
    % Calculate a rotation matrix from GCS to DCS
    %----------------------------------------------------------------------
    rotMatrixGCSToDCS = rotMatrixPCSToDCS * rotMatrixGCSToPCS;

    %----------------------------------------------------------------------
    % Calculate a slice offset in DCS
    %----------------------------------------------------------------------
    DCS_offset = rotMatrixPCSToDCS * PCS_offset; % 3 x 1

    %----------------------------------------------------------------------
    % Calculate spatial coordinates in DCS
    %----------------------------------------------------------------------
    r_range = (-floor(Nx/2):ceil(Nx/2)-1).';
    c_range = (-floor(Ny/2):ceil(Ny/2)-1).';
    s_range = (-floor(Nz/2):ceil(Nz/2)-1).';
    [r,c,s] = ndgrid(r_range, c_range, s_range); % Nx x Ny x Nz

    xyz = rotMatrixPCSToDCS * rotMatrixGCSToPCS * rotMatrixRCSToGCS * scaling_matrix * cat(2, r(:), c(:), s(:)).'; % 3 x N
    x(:, :, :, average_nr, actual_slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr, segment_nr) = reshape(xyz(1,:), [Nx Ny Nz]) + DCS_offset(1); % [m]
    y(:, :, :, average_nr, actual_slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr, segment_nr) = reshape(xyz(2,:), [Nx Ny Nz]) + DCS_offset(2); % [m]
    z(:, :, :, average_nr, actual_slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr, segment_nr) = reshape(xyz(3,:), [Nx Ny Nz]) + DCS_offset(3); % [m]

    %% Get a slice offset in PCS from ISMRMRD format
    %----------------------------------------------------------------------
    % Get a list of profiles in the current slice
    %----------------------------------------------------------------------
    index = find((raw_data.head.idx.average    == (average_nr - 1))    & ...
                 (raw_data.head.idx.slice      == (slice_nr - 1))      & ...
                 (raw_data.head.idx.contrast   == (contrast_nr - 1))   & ...
                 (raw_data.head.idx.phase      == (phase_nr - 1))      & ...
                 (raw_data.head.idx.repetition == (repetition_nr - 1)) & ...
                 (raw_data.head.idx.set        == (set_nr - 1))        & ...
                 (raw_data.head.idx.segment    == (segment_nr - 1)));

    sag_offset_ismrmrd = raw_data.head.position(1,index(1)); % [mm]
    cor_offset_ismrmrd = raw_data.head.position(2,index(1)); % [mm]
    tra_offset_ismrmrd = raw_data.head.position(3,index(1)); % [mm]

    %% Calculate a rotation matrix from GCS to PCS from ISMRMRD format
    rotMatrixGCSToPCS_ismrmrd = double([raw_data.head.phase_dir(:,index(1)) raw_data.head.read_dir(:,index(1)) raw_data.head.slice_dir(:,index(1))]);

    %% Calculate a rotation matrix from GCS to DCS
    rotMatrixGCSToDCS = rotMatrixPCSToDCS * rotMatrixGCSToPCS;

    %% Display slice information
    fprintf('======================= SLICE INFORMATION =======================\n');
    fprintf('slice_nr = %d, actual_slice_nr = %d\n', slice_nr, actual_slice_nr);
    fprintf('dNormalSag = %+g \ndNormalCor = %+g \ndNormalTra = %+g \ndRotAngle = %g [rad]\n', dNormalSag, dNormalCor, dNormalTra, dRotAngle);
    fprintf('---------------------- From Siemens raw data --------------------\n');
    fprintf('                   [sag]   %10.5f [mm]\n', sag_offset);
    fprintf('slice offset     : [cor] = %10.5f [mm]\n', cor_offset);
    fprintf('                   [tra]   %10.5f [mm]\n', tra_offset);
    fprintf('---------------------- From ISMRMRD header ----------------------\n');
    fprintf('                   [sag]   %10.5f [mm]\n', sag_offset_ismrmrd);
    fprintf('slice offset     : [cor] = %10.5f [mm]\n', cor_offset_ismrmrd);
    fprintf('                   [tra]   %10.5f [mm]\n', tra_offset_ismrmrd);
    fprintf('---------------------- From Siemens raw data --------------------\n');
    fprintf('                   [sag]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToPCS(1,1), rotMatrixGCSToPCS(1,2), rotMatrixGCSToPCS(1,3));
    fprintf('rotMatrixGCSToPCS: [cor] = [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToPCS(2,1), rotMatrixGCSToPCS(2,2), rotMatrixGCSToPCS(2,3));
    fprintf('                   [tra]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToPCS(3,1), rotMatrixGCSToPCS(3,2), rotMatrixGCSToPCS(3,3));
    fprintf('---------------------- From ISMRMRD header ----------------------\n');
    fprintf('                   [sag]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToPCS_ismrmrd(1,1), rotMatrixGCSToPCS_ismrmrd(1,2), rotMatrixGCSToPCS_ismrmrd(1,3));
    fprintf('rotMatrixGCSToPCS: [cor] = [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToPCS_ismrmrd(2,1), rotMatrixGCSToPCS_ismrmrd(2,2), rotMatrixGCSToPCS_ismrmrd(2,3));
    fprintf('                   [tra]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToPCS_ismrmrd(3,1), rotMatrixGCSToPCS_ismrmrd(3,2), rotMatrixGCSToPCS_ismrmrd(3,3));
    fprintf('-----------------------------------------------------------------\n');
    fprintf('                   [sag]   [%10.5f %10.5f %10.5f]\n', rotMatrixPCSToDCS(1,1), rotMatrixPCSToDCS(1,2), rotMatrixPCSToDCS(1,3));
    fprintf('rotMatrixPCSToDCS: [cor] = [%10.5f %10.5f %10.5f]\n', rotMatrixPCSToDCS(2,1), rotMatrixPCSToDCS(2,2), rotMatrixPCSToDCS(2,3));
    fprintf('                   [tra]   [%10.5f %10.5f %10.5f]\n', rotMatrixPCSToDCS(3,1), rotMatrixPCSToDCS(3,2), rotMatrixPCSToDCS(3,3));
    fprintf('-----------------------------------------------------------------\n');
    fprintf('                   [sag]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToDCS(1,1), rotMatrixGCSToDCS(1,2), rotMatrixGCSToDCS(1,3));
    fprintf('rotMatrixGCSToDCS: [cor] = [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToDCS(2,1), rotMatrixGCSToDCS(2,2), rotMatrixGCSToDCS(2,3));
    fprintf('                   [tra]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToDCS(3,1), rotMatrixGCSToDCS(3,2), rotMatrixGCSToDCS(3,3));
    fprintf('=================================================================\n');

    %% Get k-space data (Nkx x Nky x Nkz x Nc)
    kspace = reshape(data(:, :, :, average_nr, slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr, segment_nr, :), [Nkx Nky Nkz Nc]);

    %% Apply the noise prewhitening matrix
    kspace = ipermute(reshape(inv_L * reshape(permute(kspace,[4 1 2 3]), [Nc Nkx*Nky*Nkz]), [Nc Nkx Nky Nkz]), [4 1 2 3]);

    %% Perform Cartesian reconstruction
    %----------------------------------------------------------------------
    % Flip k-space
    %----------------------------------------------------------------------
    if RO_sign == -1
        kspace = flip(kspace,1);
    end
    if PE_sign == -1
        kspace = flip(kspace,2);
    end

    %----------------------------------------------------------------------
    % Reconstruct in x
    % Siemens use forward FFT to move from k-space to image-space
    % (k-space <=> image-space)
    %----------------------------------------------------------------------
    imc_ = fftshift(fft(ifftshift(kspace, 1), [], 1), 1); % Nkx x Nky x Nkz x Nc

    %----------------------------------------------------------------------
    % Remove oversampling
    %----------------------------------------------------------------------
    No = Nx / zpad_factor; % original matrix size in image-domain
    idx1_range = (-floor(No/2):ceil(No/2)-1).' + floor(Nkx/2) + 1;
    imc_ = imc_(idx1_range,:,:,:);  % No x Nky x Nkz x Nc

    %----------------------------------------------------------------------
    % Back to k-space (k-space <=> image-space)
    %----------------------------------------------------------------------
    kspace = fftshift(ifft(ifftshift(imc_, 1), [], 1), 1); % No x Nky x Nkz x Nc
    
    %----------------------------------------------------------------------
    % Zeropad in k-space
    %----------------------------------------------------------------------
    kspace_zpad = zeros(Nx, Ny, Nz, Nc, 'double');
    idx1_range = (-floor(No/2):ceil(No/2)-1).' + floor(Nx/2) + 1;
    idx2_range = (-floor(Nky/2):ceil(Nky/2)-1).' + floor(Ny/2) + 1;
    idx3_range = (-floor(Nkz/2):ceil(Nkz/2)-1).' + floor(Nz/2) + 1;
    kspace_zpad(idx1_range,idx2_range,idx3_range,:) = kspace;

    %----------------------------------------------------------------------
    % Reconstruct in x
    % Siemens use forward FFT to move from k-space to image-space
    % (k-space <=> image-space)
    %----------------------------------------------------------------------
    imc_ = fftshift(fft(ifftshift(kspace_zpad, 1), [], 1), 1); % Nx x Ny x Nz x Nc

    %----------------------------------------------------------------------
    % Reconstruct in y
    % Siemens use forward FFT to move from k-space to image-space
    %----------------------------------------------------------------------
    imc_ = fftshift(fft(ifftshift(imc_, 2), [], 2), 2);
    imc_ = squeeze(imc_);

    %----------------------------------------------------------------------
    % Save coil images
    %----------------------------------------------------------------------
    imc(:, :, :, average_nr, actual_slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr, segment_nr, :) = imc_;

    %----------------------------------------------------------------------
    % IFFT to k-space (k-space <=> image-space)
    %----------------------------------------------------------------------
    kspace = ifft2c(imc_);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

    %% Calculate coil sensitivity maps
    tic; fprintf('(%2d/%2d): Calculating coil sensitivity maps with Walsh method... ', idx, nr_recons);
    %----------------------------------------------------------------------
    % Calculate the calibration region of k-space
    %----------------------------------------------------------------------
    cal_shape = [32 32];
    cal_data = crop(kspace, [cal_shape Nc]);
    cal_data = bsxfun(@times, cal_data, hamming(cal_shape(1)) * hamming(cal_shape(2)).');

    %----------------------------------------------------------------------
    % Calculate coil sensitivity maps (Nx x Ny x Nc)
    %----------------------------------------------------------------------
    cal_im = fft2c(zpad(cal_data, [Nx Ny Nc]));
    csm_ = ismrm_estimate_csm_walsh(cal_im);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

    %% Save CSM
    csm(:, :, :, average_nr, actual_slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr, segment_nr, :) = csm_;

    %% Perform optimal coil combination
    im(:, :, :, average_nr, actual_slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr, segment_nr) = sum(conj(csm_) .* imc_, 3) ./ sum(abs(csm_).^2, 3);
end

end
