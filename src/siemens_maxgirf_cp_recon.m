function [im_cpr_multislice,header,r_dcs_multislice,output] = siemens_maxgirf_cp_recon(ismrmrd_noise_fullpath, ismrmrd_data_fullpath, siemens_dat_fullpath, B0map, user_opts)
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/16/2022, Last modified: 01/16/2022

%% Define constants
gamma = 4257.59 * (1e4 * 2 * pi); % gyromagnetic ratio for 1H [rad/sec/T]

%% Read k-space data (ISMRMRD format)
start_time = tic;
tstart = tic; fprintf('Reading an ISMRMRD file: %s... ', ismrmrd_data_fullpath);
if exist(ismrmrd_data_fullpath, 'file')
    dset = ismrmrd.Dataset(ismrmrd_data_fullpath, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
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
% Experimental conditions
%--------------------------------------------------------------------------
B0 = header.experimentalConditions.H1resonanceFrequency_Hz * (2 * pi / gamma); % [Hz] * [2pi rad/cycle] / [rad/sec/T] => [T]

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

encoded_resolution = encoded_fov ./ [Nkx Nky Nkz]; % [mm]
recon_resolution = recon_fov ./ [Nx Ny Nz]; % [mm]

%--------------------------------------------------------------------------
% Trajectory description (HargreavesVDS2000)
%--------------------------------------------------------------------------
if strcmp(header.encoding.trajectory, 'spiral')
    Gmax              = header.encoding.trajectoryDescription.userParameterDouble(1).value;      % [G/cm]
    Smax              = header.encoding.trajectoryDescription.userParameterDouble(2).value;      % [G/cm/sec]
    FOV               = header.encoding.trajectoryDescription.userParameterDouble(3).value;      % [cm]
    krmax             = header.encoding.trajectoryDescription.userParameterDouble(4).value;      % [cycle/cm]
    interleaves       = header.encoding.trajectoryDescription.userParameterLong(1).value;        % interleaves
    fov_coefficients  = header.encoding.trajectoryDescription.userParameterLong(2).value;        % fov_coefficients
    sampling_time     = header.encoding.trajectoryDescription.userParameterLong(3).value * 1e-9; % SamplingTime_ns [sec]
    fprintf('==================== Trajectory description =====================\n');
    fprintf('Gmax          = %5.3f [G/cm]\n'    , Gmax);
    fprintf('Smax          = %5.3f [G/cm/sec]\n', Smax);
    fprintf('FOV           = %5.2f [cm]\n'      , FOV);
    fprintf('interleaves   = %d\n'              , interleaves);
    fprintf('sampling_time = %5.3f [sec]\n'     , sampling_time);
    fprintf('=================================================================\n');
end

%% Parse the ISMRMRD header
tstart = tic; fprintf('Parsing the ISMRMRD header... ');
raw_data = dset.readAcquisition(); % read all the acquisitions
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

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

%% Define parameters for reconstruction
%--------------------------------------------------------------------------
% Set the number of samples in the row (r) direction of the RCS
%--------------------------------------------------------------------------
if isfield(user_opts, 'N1')
    N1 = user_opts.N1;
else
    N1 = Nkx;
end

%--------------------------------------------------------------------------
% Set the number of samples in the column (c) direction of the RCS
%--------------------------------------------------------------------------
if isfield(user_opts, 'N2')
    N2 = user_opts.N2;
else
    N2 = Nky;
end

N3 = Nkz;                % number of samples in the slice (s) direction of the RCS
Ni = nr_phase_encodings; % number of spiral interleaves
N  = N1 * N2;            % total number of voxels in image-space

%--------------------------------------------------------------------------
% Set the number of spatial basis functions
%--------------------------------------------------------------------------
if isfield(user_opts, 'Nl')
    Nl = user_opts.Nl;
else
    Nl = 19;
end

%--------------------------------------------------------------------------
% Calculate the number of samples per spiral arm
%--------------------------------------------------------------------------
if ~isempty(user_opts.discard_pre)
    discard_pre = user_opts.discard_pre;
end
if ~isempty(user_opts.discard_post)
    discard_post = user_opts.discard_post;
end
Nk = number_of_samples - discard_pre - discard_post;

%--------------------------------------------------------------------------
% Select the channels
%--------------------------------------------------------------------------
if isfield(user_opts, 'selected_channels')
    selected_channels = user_opts.selected_channels;
else
    selected_channels = (1:nr_channels).';
end
Nc = length(selected_channels); % number of coils

%--------------------------------------------------------------------------
% Set static_B0_correction
%--------------------------------------------------------------------------
if isfield(user_opts, 'static_B0_correction')
    static_B0_correction = user_opts.static_B0_correction;
else
    static_B0_correction = 1;
end

%--------------------------------------------------------------------------
% Set the maximum rank of the SVD approximation of a higher-order encoding matrix
%--------------------------------------------------------------------------
if isfield(user_opts, 'Lmax')
    Lmax = user_opts.Lmax;
else
    Lmax = 50;
end

%--------------------------------------------------------------------------
% Set the oversampling parameter for randomized SVD
%--------------------------------------------------------------------------
os = 5; 

%% Prepare a static off-resonance map [Hz]
if isempty(B0map)
    B0map = complex(zeros(N1, N2, nr_slices, 'double'));
end

%% Read a Siemens .dat file
if exist(siemens_dat_fullpath, 'file')
    fprintf('Reading a Siemens .dat file: %s\n', siemens_dat_fullpath);
    twix = mapVBVD(siemens_dat_fullpath);
    if length(twix) > 1
        twix = twix{end};
    end
end

%% Get a slice normal vector from Siemens TWIX format
if exist(siemens_dat_fullpath, 'file')
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
end

%% Calculate a scaling matrix
scaling_matrix = diag(encoded_resolution) * 1e-3; % [mm] * [m/1e3mm] => [m]

%% Calculate a transformation matrix from RCS to GCS [r,c,s] <=> [PE,RO,SL]
R_rcs2gcs = [0    1    0 ; % [PE]   [0 1 0] * [r]
             1    0    0 ; % [RO] = [1 0 0] * [c]
             0    0    1]; % [SL]   [0 0 1] * [s]

%% Get a rotation matrix from GCS to PCS (ISMRMRD format)
phase_dir = double(raw_data.head.phase_dir(:,1));
read_dir  = double(raw_data.head.read_dir(:,1));
slice_dir = double(raw_data.head.slice_dir(:,1));
R_gcs2pcs_ismrmrd = [phase_dir read_dir slice_dir];

%% Calculate a rotation matrix from GCS to PCS
if exist(siemens_dat_fullpath, 'file')
    [R_gcs2pcs,phase_sign,read_sign] = siemens_calculate_matrix_gcs_to_pcs(dNormalSag, dNormalCor, dNormalTra, dRotAngle);
else
    phase_sign = user_opts.phase_sign;
    read_sign = user_opts.read_sign;
    R_gcs2pcs = [phase_sign * phase_dir read_sign * read_dir slice_dir];
end

%% Calculate a rotation matrix from PCS to DCS
R_pcs2dcs = siemens_calculate_matrix_pcs_to_dcs(patient_position);

%% Calculate a rotation matrix from GCS to DCS
R_gcs2dcs = R_pcs2dcs * R_gcs2pcs;

%% Calculate nominal gradient waveforms in GCS [G/cm]
%--------------------------------------------------------------------------
% Calculate one spiral interleaf: k in [cycle/cm], g in [G/cm]
%--------------------------------------------------------------------------
Fcoeff = [FOV -FOV * (1 - user_opts.vds_factor / 100)]; % FOV decreases linearly from 30 to 15cm
krmax  = 1 / (2 * (FOV / Nkx));                         % [cycle/cm]
res    = 1 / (2 * krmax) * 10;                          % resolution [mm]
[k_spiral_arm,g_spiral_arm,s_spiral_arm,time] = vdsmex(interleaves, Fcoeff, res, Gmax, Smax, sampling_time, 10000000);
g_spiral_arm = -g_spiral_arm(1:Nk,:); % Nk x 2

%--------------------------------------------------------------------------
% Rotate the spiral interleaf by 360/Ni degrees every TR
% (Re{g_spiral} + 1j * Im{g_spiral}) * (cos(arg) - 1j * sin(arg))
% RO: real =>  Re{g_spiral} * cos(arg) + Im{g_spiral} * sin(arg)
% PE: imag => -Re{g_spiral} * sin(arg) + Im{g_spiral} * cos(arg)
%--------------------------------------------------------------------------
gu_nominal = zeros(Nk, Ni, 'double'); % Nk x Ni [G/cm] PE (gu)
gv_nominal = zeros(Nk, Ni, 'double'); % Nk x Ni [G/cm] RO (gv)
gw_nominal = zeros(Nk, Ni, 'double'); % Nk x Ni [G/cm] SL (gw)
for i = 1:Ni
    arg = 2 * pi / Ni * (i - 1); % [rad]
    gv_nominal(:,i) =  g_spiral_arm(:,1) * cos(arg) + g_spiral_arm(:,2) * sin(arg); % RO (gv)
    gu_nominal(:,i) = -g_spiral_arm(:,1) * sin(arg) + g_spiral_arm(:,2) * cos(arg); % PE (gu)
end

%% Calculate GIRF-predicted gradient waveforms in GCS [G/cm]
g_nominal = cat(3, gu_nominal, gv_nominal);
tRR = 0; % custom clock-shift
sR.R = R_gcs2dcs;
sR.T = header.acquisitionSystemInformation.systemFieldStrength_T;
[~,g_predicted] = apply_GIRF(g_nominal, dt, sR, tRR); % k:[cycle/cm] and g:[G/cm]

%% Change the sign of GIRF-predicted gradient waveforms in GCS [G/cm] [PE,RO,SL]
g_gcs = zeros(Nk, 3, Ni, 'double');
g_gcs(:,1,:) = phase_sign * g_predicted(:,:,1); % [G/cm] PE (gu)
g_gcs(:,2,:) = read_sign  * g_predicted(:,:,2); % [G/cm] RO (gv)
g_gcs(:,3,:) = g_predicted(:,:,3);              % [G/cm] SL (gw)

%% Calculate GIRF-predicted gradient waveforms in DCS [G/cm] [x,y,z]
g_dcs = zeros(Nk, 3, Ni, 'double');
for i = 1:Ni
    g_dcs(:,:,i) = (R_gcs2dcs * g_gcs(:,:,i).').'; % Nk x 3
end

%% Calculate the time courses of phase coefficients (Nk x Nl x Ni) [rad/m], [rad/m^2], [rad/m^3]
tstart = tic; fprintf('Calculating the time courses of phase coefficients... ');
k = calculate_concomitant_field_coefficients(reshape(g_dcs(:,1,:), [Nk Ni]), reshape(g_dcs(:,2,:), [Nk Ni]), reshape(g_dcs(:,3,:), [Nk Ni]), Nl, B0, gamma, dt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Calculate a time vector [sec]
t = (0:Nk-1).' * dt; % Nk x 1 [sec]

%% Calculate GIRF-predicted k-space trajectories in GCS [rad/m] [PE,RO,SL]
%--------------------------------------------------------------------------
% Numerically integrate the coefficients
% [rad/sec/T] * [G/cm] * [T/1e4G] * [1e2cm/m] * [sec] => [rad/m]
%--------------------------------------------------------------------------
k_gcs = cumsum(gamma * g_gcs * 1e-2 * dt); % [rad/m]

%% Calculate GIRF-predicted k-space trajectories in DCS [rad/m] [x,y,z]
k_dcs = zeros(Nk, 3, Ni, 'double');
for i = 1:Ni
    k_dcs(:,:,i) = (R_gcs2dcs * k_gcs(:,:,i).').'; % Nk x 3
end

%% Calculate GIRF-predicted k-space trajectories in RCS [rad/m] [R,C,S]
k_rcs = zeros(Nk, 3, Ni, 'double');
for i = 1:Ni
    k_rcs(:,:,i) = (R_rcs2gcs.' * k_gcs(:,:,i).').'; % Nk x 3
end

%% Initialize structure for NUFFT for all interleaves (NUFFT reconstruction)
tstart = tic; fprintf('Initializing structure for NUFFT for all interleaves... ');
% scaled to [-0.5,0.5] and then [-pi,pi]
% [rad/m] / ([cycle/cm] * [2pi rad/cycle] * [1e2cm/m]) => [rad/m] / [rad/m] = [unitless] ([-0.5,0.5]) * 2pi => [-pi,pi]
% The definition of FFT is opposite in NUFFT
om = -cat(2, reshape(k_rcs(:,1,:), [Nk*Ni 1]), reshape(k_rcs(:,2,:), [Nk*Ni 1])) / (2 * krmax * 1e2); % Nk*Ni x 2
Nd = [N1 N2]; % matrix size
Jd = [6 6];   % kernel size
Kd = Nd * 2;  % oversampled matrix size
nufft_st = nufft_init(om, Nd, Jd, Kd, Nd/2, 'minmax:kb');
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Initialize structure for NUFFT per interleaf (MaxGIRF reconstruction)
st = cell(Ni,1);
for i = 1:Ni
    tstart = tic; fprintf('(%2d/%2d): Initializing structure for NUFFT per interleaf... ', i, Ni);
    % scaled to [-0.5,0.5] and then [-pi,pi]
    % [rad/m] / ([cycle/cm] * [2pi rad/cycle] * [1e2cm/m]) => [rad/m] / [rad/m] = [unitless]
    om = -cat(2, k_rcs(:,1,i), k_rcs(:,2,i)) / (2 * krmax * 1e2); % Nk x 2
    Nd = [N1 N2];   % matrix size
    Jd = [6 6];     % kernel size
    Kd = Nd * 2;    % oversampled matrix size
    st{i} = nufft_init(om, Nd, Jd, Kd, Nd/2, 'minmax:kb');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Calculate a density compensation function
tstart = tic; fprintf('Calculating a density compensation function based on sdc3_MAT.c... ');
% [rad/m] / ([cycle/cm] * [2pi rad/cycle] * [1e2cm/m]) => [unitless]
coords = permute(k_rcs, [2 1 3]) / (2 * krmax * 2 * pi * 1e2); % Nk x 3 x Ni => 3 x Nk x Ni
coords(3,:,:) = 0;
numIter = 25;
effMtx  = Nkx;
verbose = 0;
DCF = sdc3_MAT(coords, numIter, effMtx, verbose, 2.1);
w = DCF / max(DCF(:));
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Calculate the receiver noise matrix
[Psi,inv_L] = calculate_noise_decorrelation_matrix(ismrmrd_noise_fullpath);

%% Perform NUFFT reconstruction per slice
nr_recons = nr_slices * nr_contrasts * nr_phases * nr_repetitions * nr_sets * nr_segments;
im_cpr_multislice = complex(zeros(N1, N2, nr_slices, Lmax, 'double'));
r_dcs_multislice = zeros(N, 3, nr_slices, 'double');
output = repmat(struct('s', [], 'computation_time_svd', [], 'computation_time_cpr', []), [nr_recons 1]);

for idx = 1:nr_recons
    %% Get information about the current slice
    [slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr, segment_nr] = ind2sub([nr_slices nr_contrasts nr_phases nr_repetitions nr_sets nr_segments], idx);
    fprintf('(%2d/%2d): Reconstructing slice (slice = %2d, contrast = %2d, phase = %2d, repetition = %2d, set = %2d, segment = %2d)\n', idx, nr_recons, slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr, segment_nr);

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
                        (raw_data.head.idx.set        == (set_nr - 1))        & ...
                        (raw_data.head.idx.segment    == (segment_nr - 1)));

    %% Get a slice offset in PCS from ISMRMRD format
    sag_offset_ismrmrd = double(raw_data.head.position(1,slice_nr)); % [mm]
    cor_offset_ismrmrd = double(raw_data.head.position(2,slice_nr)); % [mm]
    tra_offset_ismrmrd = double(raw_data.head.position(3,slice_nr)); % [mm]

    %% Get a slice offset in PCS from Siemens TWIX format
    if exist(siemens_dat_fullpath, 'file')
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
    end

    %% Use a slice offset in PCS from Siemens TWIX format
    pcs_offset = [sag_offset_twix; cor_offset_twix; tra_offset_twix] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Calculate spatial coordinates in DCS [m]
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

    %% Get k-space data (Nk x Ni x Nc)
    %----------------------------------------------------------------------
    % Calculate the index range of readout samples for reconstruction
    %----------------------------------------------------------------------
    index_range = ((discard_pre + 1):(number_of_samples - discard_post)).';

    kspace = complex(zeros(Nk, Ni, Nc, 'double'));
    for idx1 = 1:length(profile_list)
        tstart = tic; fprintf('(%2d/%2d): Reading k-space data (%d/%d)... ', idx, nr_recons, idx1, length(profile_list));
        %------------------------------------------------------------------
        % Determine the interleaf number
        %------------------------------------------------------------------
        interleaf_nr = raw_data.head.idx.kspace_encode_step_1(profile_list(idx1)) + 1;

        %------------------------------------------------------------------
        % Prewhiten k-space data
        %------------------------------------------------------------------
        profile = raw_data.data{profile_list(idx1)}; % number_of_samples x nr_channels
        profile = (inv_L * profile.').';

        %------------------------------------------------------------------
        % Calculate the average of k-space data
        %------------------------------------------------------------------
        profile = profile / nr_averages;

        %------------------------------------------------------------------
        % Accumulate k-space
        %------------------------------------------------------------------
        kspace(:,interleaf_nr,:) = kspace(:,interleaf_nr,:) + reshape(profile(index_range,selected_channels), [Nk 1 Nc]);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %% Calculate concomitant field basis functions (N x Nl) [m], [m^2], [m^3]
    tstart = tic; fprintf('(%2d/%2d): Calculating concomitant field basis functions... ', idx, nr_recons);
    p = calculate_concomitant_field_basis(r_dcs(:,1), r_dcs(:,2), r_dcs(:,3), Nl);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Perform NUFFT reconstruction
    imc_nufft = complex(zeros(N1, N2, Nc, 'double'));
    scale_factor = 1 / sqrt(prod(Nd));
    for c = 1:Nc
        tstart = tic; fprintf('(%2d/%2d): NUFFT reconstruction (c=%2d/%2d)... ', idx, nr_recons, c, Nc);
        imc_nufft(:,:,c) = nufft_adj(reshape(kspace(:,:,c) .* w, [Nk*Ni 1]), nufft_st) * scale_factor;
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %% Calculate coil sensitivity maps
    tstart = tic; fprintf('(%2d/%2d): Calculating coil sensitivity maps with Walsh method... ', idx, nr_recons);
    %----------------------------------------------------------------------
    % IFFT to k-space (k-space <=> image-space)
    %----------------------------------------------------------------------
    kspace_gridded = imc_nufft;
    for dim = 1:2
        kspace_gridded = sqrt(size(kspace_gridded,dim)) * fftshift(ifft(ifftshift(kspace_gridded, dim), [], dim), dim);
    end

    %----------------------------------------------------------------------
    % Calculate the calibration region of k-space
    %----------------------------------------------------------------------
    cal_shape = [32 32];
    cal_data = crop(reshape(kspace_gridded, [N1 N2 Nc]), [cal_shape Nc]);
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

    %% Calculate the SVD of the higher-order encoding matrix (Nk x N)
    start_time_svd = tic;
    U = complex(zeros(Nk, Lmax, Ni, 'double'));
    V = complex(zeros(N, Lmax, Ni, 'double'));
    s = zeros(Lmax, Ni, 'double');
    for i = 1:Ni
        tstart = tic; fprintf('(%2d/%2d): Calculating randomized SVD (i=%2d/%2d)... ', idx, nr_recons, i, Ni);
        [U_,S_,V_] = calculate_rsvd_higher_order_encoding_matrix(k(:,4:end,i), p(:,4:end), Lmax, os, reshape(B0map(:,:,actual_slice_nr), [N 1]), t, static_B0_correction);
        U(:,:,i) = U_(:,1:Lmax); % U: Nk x Lmax+os => Nk x Lmax
        V(:,:,i) = V_(:,1:Lmax) * S_(1:Lmax,1:Lmax)'; % V: N x Lmax+os => N x Lmax
        s(:,i) = diag(S_(1:Lmax,1:Lmax));
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
    computation_time_svd = toc(start_time_svd);

    %% Perform CP-based MaxGIRF reconstruction (conjugate phase reconstruction)
    start_time_cpr = tic;
    b = complex(zeros(N, Lmax, 'double'));
    for i = 1:Ni
        tstart = tic; fprintf('(%d/%d): Performing CP-based MaxGIRF reconstruction (i=%d/%d)... ', idx, nr_recons, i, Ni);

        for c = 1:Nc
            %--------------------------------------------------------------
            % Caclulate d_{i,c}
            %--------------------------------------------------------------
            d = kspace(:,i,c); % kspace: Nk x Ni x Nc

            %--------------------------------------------------------------
            % Calculate sum_{ell=1}^L diag(V(ell,i)) * Fi^H * diag(conj(U(ell,i)))
            %--------------------------------------------------------------
            AHd = complex(zeros(N, Lmax, 'double'));
            AHd_ = complex(zeros(N, 1, 'double'));
            scale_factor = 1 / sqrt(prod(st{i}.Nd));
            for ell = 1:Lmax
                % Preconditioning with density compensation
                FHDuHd = nufft_adj((conj(U(:,ell,i)) .* d) .* w(:,i), st{i}) * scale_factor;
                AHd_ = AHd_ + V(:,ell,i) .* reshape(FHDuHd, [N 1]);
                AHd(:,ell) = AHd_;
            end

            %--------------------------------------------------------------
            % Calculate Sc^H * Ei^H * d_{i,c}
            %--------------------------------------------------------------
            for ell = 1:Lmax
                AHd(:,ell) = reshape(conj(csm(:,:,c)), [N 1]) .* AHd(:,ell);
            end

            %--------------------------------------------------------------
            % Calculate b (N x 1)
            %--------------------------------------------------------------
            for ell = 1:Lmax
                b(:,ell) = b(:,ell) + AHd(:,ell);
            end
        end
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
    computation_time_cpr = toc(start_time_cpr);

    %% Collect the output
    im_cpr_multislice(:,:,idx,:) = reshape(b, [N1 N2 1 Lmax]);
    output(idx).s = s;
    output(idx).computation_time_svd = computation_time_svd;
    output(idx).computation_time_cpr = computation_time_cpr;
end
end