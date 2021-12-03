% demo_maxgirf_spiral_se_human.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 05/23/2020, Last modified: 09/18/2021

%% Notation
%--------------------------------------------------------------------------
% We denote X, Y, Z as the horizontal, vertical, and through-magnet axes,
% respectively in the physical coordinate system, and the corresponding
% coordinates as x, y, and z. We also denote U, V, and W as the phase encoding,
% readout, and slice axes, respectively in the logical coordinate
% system, and the corresponding coordinates as u, v, and w.
%
% HFS position:
%
%                     Anterior
%                        ^ (+y)
%                        |
%           Right        |            Left
%                        |
%                  Head  +--------> (+x)
%                       /
%                      /
%                     /
%               Foot v (+z)
%                      Posterior
%
%--------------------------------------------------------------------------

%% Set directory names
computer_type = computer;
if strcmp(computer_type, 'PCWIN64')
    src_directory = 'E:\lowfield_maxgirf';
    ismrmrd_directory = 'D:\ismrmrd\ismrmrd';
elseif strcmp(computer_type, 'GLNXA64')
    src_directory = '/server/home/nlee/lowfield_maxgirf';
    ismrmrd_directory = '/server/home/nlee/ismrmrd';
end

%% Add paths
addpath(genpath(src_directory));
addpath(genpath(ismrmrd_directory));
warning off;

%% Define a fullpath to each filename
spiral_data_fullpath  = fullfile(data_directory, 'h5'   , sprintf('%s.h5', spiral_filename));
spiral_noise_fullpath = fullfile(data_directory, 'noise', sprintf('noise_%s.h5', spiral_filename));
spiral_dat_fullpath   = fullfile(data_directory, sprintf('%s.dat', spiral_filename));
spiral_dat_fullpath = '';

%% Define an output directory
[filepath,data_filename] = fileparts(data_directory);
underscore_loc = strfind(spiral_filename, '_');
output_filename = sprintf('%s_%s_osf%d_B0_correction%d_concomitant_correction%d_Lmax%d_L%d', data_filename, spiral_filename(underscore_loc(3)+1:end), osf, static_B0_correction, concomitant_field_correction, Lmax, L);
output_directory = fullfile(src_directory, output_directory_filename, output_filename);
mkdir(output_directory);

%% Define constants
gamma = 4257.59 * (1e4 * 2 * pi); % gyromagnetic ratio for 1H [rad/sec/T]

SAGITTAL   = 0; % Patient axis perpendicular to the sagittal plane
CORONAL    = 1; % Patient axis perpendicular to the coronal plane
TRANSVERSE = 2; % Patient axis perpendicular to the transverse plane

%% Load a smooth static off-resonance map [Hz]
load(B0map_fullpath);
B0map_orig = B0map_nlinv;

%% Read spiral k-space data (ISMRMRD format)
start_time = tic;
tic; fprintf('Reading an ISMRMRD file: %s... ', spiral_data_fullpath);
if exist(spiral_data_fullpath, 'file')
    dset = ismrmrd.Dataset(spiral_data_fullpath, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , spiral_data_fullpath);
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
field_of_view_mm(1) = header.encoding.encodedSpace.fieldOfView_mm.x; % RO
field_of_view_mm(2) = header.encoding.encodedSpace.fieldOfView_mm.y; % PE
field_of_view_mm(3) = header.encoding.encodedSpace.fieldOfView_mm.z; % SL
matrix_size(1)      = header.encoding.encodedSpace.matrixSize.x;     % number of samples in idx1 (r)
matrix_size(2)      = header.encoding.encodedSpace.matrixSize.y;     % number of samples in idx2 (c)
matrix_size(3)      = header.encoding.encodedSpace.matrixSize.z;     % number of samples in idx3 (s)

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
    fprintf('Gmax          = %5.3f [G/cm]\n' , Gmax);
    fprintf('Smax          = %d [G/cm/sec]\n', Smax);
    fprintf('FOV           = %4.2f [cm]\n'   , FOV);
    fprintf('interleaves   = %d\n'           , interleaves);
    fprintf('sampling_time = %g [sec]\n'     , sampling_time);
    fprintf('=================================================================\n');
end

%% Parse the ISMRMRD header
tstart = tic; fprintf('Parsing the ISMRMRD header... ');
raw_data = dset.readAcquisition(); % read all the acquisitions
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% ISMRMRD header
%--------------------------------------------------------------------------
nr_samples         = double(max(raw_data.head.number_of_samples));
discard_pre        = double(max(raw_data.head.discard_pre));
discard_post       = double(max(raw_data.head.discard_post));
nr_channels        = double(max(raw_data.head.active_channels));
nr_phase_encodings = double(max(raw_data.head.idx.kspace_encode_step_1)) + 1;
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
dt = single(max(raw_data.head.sample_time_us)) * 1e-6; % [usec] * [sec/1e-6 usec] => [sec]

%--------------------------------------------------------------------------
% Calculate the readout duration [sec]
%--------------------------------------------------------------------------
T = nr_samples * dt; % readout duration [sec]

%% Display ISMRMRD header
fprintf('========================= ISMRMRD header ========================\n');
fprintf('nr_samples         = %d\n', nr_samples);
fprintf('discard_pre        = %d\n', discard_pre);
fprintf('discard_post       = %d\n', discard_post);
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
fprintf('readout duration   = %g [msec]\n', T * 1e3);
fprintf('=================================================================\n');

%% Process noise only ismrmrd data
[Psi,inv_L] = calculate_receiver_noise_matrix(spiral_noise_fullpath);

%% Read a Siemens dat file
if exist(spiral_dat_fullpath, 'file')
    fprintf('Reading a Siemens .dat file: %s\n', spiral_dat_fullpath);
    twix = mapVBVD(spiral_dat_fullpath);
    if length(twix) > 1
        twix = twix{end};
    end
else    
end

%% Define parameters used in reconstruction
N1 = osf * matrix_size(1); % number of samples in RCS (r)
N2 = osf * matrix_size(2); % number of samples in RCS (c)
N3 = matrix_size(3);       % number of samples in RCS (s)
Ni = nr_phase_encodings;   % number of spiral interleaves

%--------------------------------------------------------------------------
% Calculate the number of samples per spiral arm
%--------------------------------------------------------------------------
if ~isempty(user_opts.discard_pre)
    discard_pre = user_opts.discard_pre;
end

if ~isempty(user_opts.discard_post)
    discard_post = user_opts.discard_post;
end

Nk = nr_samples - discard_pre - discard_post; % number of samples in one spiral interleaf
Np = Nk * Ni;  % total number of fully-sampled k-space samples per image
N  = N1 * N2;  % total number of voxels in image-space
Nl = 19;       % number of spatial basis functions

%--------------------------------------------------------------------------
% Calculate the range of k-space indices for reconstruction
%--------------------------------------------------------------------------
start_index = discard_pre + 1;
end_index   = nr_samples - discard_post;
index_range = (start_index:end_index).';

%--------------------------------------------------------------------------
% Select the channels
%--------------------------------------------------------------------------
if isempty(channel_range)
    channel_range = (1:nr_channels).';
end
Nc = length(channel_range); % number of coils

%% Zeropad a static off-resonance map
if osf > 1
    [N1_,N2_,N3_] = size(B0map_orig);
    idx1_range = (-floor(N1_/2):ceil(N1_/2)-1).' + floor(N1/2) + 1;
    idx2_range = (-floor(N2_/2):ceil(N2_/2)-1).' + floor(N2/2) + 1;
    idx3_range = (1:N3_).';
    B0map = zeros(N1, N2, N3_, 'double');
    B0map(idx1_range,idx2_range,idx3_range) = B0map_orig;
else
    B0map = B0map_orig;
end

%% Reconstruct images per slice
nr_recons = nr_slices * nr_contrasts * nr_phases * nr_repetitions * nr_sets * nr_segments;
% 9 for sagittal
%actual: 1 2 3 4 5 6 7 8 9 10 11
%slice : 1,3,5,7,9,11,2,4,6,8,10
%10 for axial

for idx = 10%1:nr_recons
    %% Get information about the current slice
    [slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr, segment_nr] = ind2sub([nr_slices nr_contrasts nr_phases nr_repetitions nr_sets nr_segments], idx);
    fprintf('(%2d/%2d): Reconstructing slice (slice = %2d, contrast = %2d, phase = %2d, repetition = %2d, set = %2d, segment = %2d)\n', idx, nr_recons, slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr, segment_nr);

    %% Get a list of profiles in the current slice (all averages)
    profile_list = find((raw_data.head.idx.slice      == (slice_nr - 1))      & ...
                        (raw_data.head.idx.contrast   == (contrast_nr - 1))   & ...
                        (raw_data.head.idx.phase      == (phase_nr - 1))      & ...
                        (raw_data.head.idx.repetition == (repetition_nr - 1)) & ...
                        (raw_data.head.idx.set        == (set_nr - 1))        & ...
                        (raw_data.head.idx.segment    == (segment_nr - 1)));

    %% Get k-space data (Nk x Ni x Nc)
    kspace = complex(zeros(Nk, Ni, Nc, 'single'));
    for idx1 = 1:length(profile_list)
        interleaf_nr = raw_data.head.idx.kspace_encode_step_1(profile_list(idx1)) + 1;
        kspace(:,interleaf_nr,:) = kspace(:,interleaf_nr,:) + reshape(raw_data.data{profile_list(idx1)}(index_range,channel_range), [Nk 1 Nc]);
    end

    %% Prewhiten k-space data
    tstart = tic; fprintf('Prewhitening k-space data... ');
    kspace = ipermute(reshape(inv_L(channel_range,channel_range) * reshape(permute(kspace, [3 1 2]), [Nc Nk*Ni]), [Nc Nk Ni]), [3 1 2]);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

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

    if exist(spiral_dat_fullpath, 'file')
        %% Get a slice normal vector from Siemens raw data
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

        % It seems that there is a bug in the pulse sequence
        main_orientation = fGSLClassOri(dNormalSag, dNormalCor, dNormalTra);
        if main_orientation == SAGITTAL
            %PCS_offset(3) = 0;
            PCS_offset(3) = PCS_offset(3) - 0.013; %% 0.013 is the offset of the center slice
        end
    else
        %% Get a slice offset in PCS from ISMRMRD format
        %------------------------------------------------------------------
        % Get a list of profiles in the current slice
        %------------------------------------------------------------------
        sag_offset = double(raw_data.head.position(1,slice_nr)); % [mm]
        cor_offset = double(raw_data.head.position(2,slice_nr)); % [mm]
        tra_offset = double(raw_data.head.position(3,slice_nr)); % [mm]
        PCS_offset = [sag_offset; cor_offset; tra_offset] * 1e-3; % [mm] * [m/1e3mm] => [m]

        % It seems that there is a bug in the pulse sequence
        if main_orientation == SAGITTAL
            PCS_offset(3) = PCS_offset(3) - 0.013; %% 0.013 is the offset of the center slice
        end
    end

    %% Calculate spatial coordinates in DCS [m]
    %----------------------------------------------------------------------
    % Calculate a transformation matrix from RCS to GCS [r,c,s] <=> [PE,RO,SL] 
    %----------------------------------------------------------------------
    rotMatrixRCSToGCS = [0    1    0 ; % [PE]   [0 1 0] * [r]
                         1    0    0 ; % [RO] = [1 0 0] * [c]
                         0    0    1]; % [SL]   [0 0 1] * [s]

    if exist(spiral_dat_fullpath, 'file')
        %------------------------------------------------------------------
        % Calculate a rotation matrix from GCS to PCS from Siemens raw data
        %------------------------------------------------------------------
        [rotMatrixGCSToPCS,PE_sign,RO_sign,main_orientation] = calcMatrixGCSToPCS(dNormalSag, dNormalCor, dNormalTra, dRotAngle);
        
        
    else
        %------------------------------------------------------------------
        % Get a rotation matrix from GCS to PCS (ISMRMRD format)
        %------------------------------------------------------------------
        rotMatrixGCSToPCS = double([raw_data.head.phase_dir(:,slice_nr) raw_data.head.read_dir(:,slice_nr) raw_data.head.slice_dir(:,slice_nr)]);
        rotMatrixGCSToPCS(:,1) = PE_sign * rotMatrixGCSToPCS(:,1);
        rotMatrixGCSToPCS(:,2) = RO_sign * rotMatrixGCSToPCS(:,2);
    end

    %----------------------------------------------------------------------
    % Calculate a rotation matrix from PCS to DCS
    %----------------------------------------------------------------------
    rotMatrixPCSToDCS = calcMatrixPCSToDCS(patient_position);

    %----------------------------------------------------------------------
    % Calculate a rotation matrix from GCS to DCS
    %----------------------------------------------------------------------
    rotMatrixGCSToDCS = rotMatrixPCSToDCS * rotMatrixGCSToPCS;

    %----------------------------------------------------------------------
    % Calculate a scaling matrix
    %----------------------------------------------------------------------
    scaling_matrix = diag(field_of_view_mm ./ matrix_size * 1e-3); % [mm] * [m/1e3mm] => [m]

    %----------------------------------------------------------------------
    % Calculate a slice offset in DCS
    %----------------------------------------------------------------------
    DCS_offset = rotMatrixPCSToDCS * PCS_offset; % 3 x 1

    %----------------------------------------------------------------------
    % Calculate spatial coordinates in DCS
    %----------------------------------------------------------------------
    idx1_range = (-floor(N1/2):ceil(N1/2)-1).';
    idx2_range = (-floor(N2/2):ceil(N2/2)-1).';
    idx3_range = (-floor(N3/2):ceil(N3/2)-1).';
    [I1,I2,I3] = ndgrid(idx1_range, idx2_range, idx3_range);
    r_dcs = (repmat(DCS_offset, [1 N]) + rotMatrixPCSToDCS * rotMatrixGCSToPCS * rotMatrixRCSToGCS * scaling_matrix * cat(2, I1(:), I2(:), I3(:)).').'; % N x 3
    x = r_dcs(:,1); % N x 1 [m]
    y = r_dcs(:,2); % N x 1 [m]
    z = r_dcs(:,3); % N x 1 [m]

    %% Display slice information
    fprintf('======================= SLICE INFORMATION =======================\n');
    fprintf('main_orientation = %d (SAGITTAL/CORONAL/TRANSVERSE = 0/1/2)\n', main_orientation);
    fprintf('slice_nr = %d, actual_slice_nr = %d\n', slice_nr, actual_slice_nr);
    %fprintf('dNormalSag = %+g \ndNormalCor = %+g \ndNormalTra = %+g \ndRotAngle = %g [rad]\n', dNormalSag, dNormalCor, dNormalTra, dRotAngle);
    fprintf('---------------------- From Siemens raw data --------------------\n');
    fprintf('                   [Sag]   %10.5f [mm]\n', sag_offset);
    fprintf('slice offset(PCS): [Cor] = %10.5f [mm]\n', cor_offset);
    fprintf('                   [Tra]   %10.5f [mm]\n', tra_offset);
    fprintf('---------------------- From Siemens raw data --------------------\n');
    fprintf('                   [Sag]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToPCS(1,1), rotMatrixGCSToPCS(1,2), rotMatrixGCSToPCS(1,3));
    fprintf('rotMatrixGCSToPCS: [Cor] = [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToPCS(2,1), rotMatrixGCSToPCS(2,2), rotMatrixGCSToPCS(2,3));
    fprintf('                   [Tra]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToPCS(3,1), rotMatrixGCSToPCS(3,2), rotMatrixGCSToPCS(3,3));
    fprintf('-----------------------------------------------------------------\n');
    fprintf('                   [Sag]   [%10.5f %10.5f %10.5f]\n', rotMatrixPCSToDCS(1,1), rotMatrixPCSToDCS(1,2), rotMatrixPCSToDCS(1,3));
    fprintf('rotMatrixPCSToDCS: [Cor] = [%10.5f %10.5f %10.5f]\n', rotMatrixPCSToDCS(2,1), rotMatrixPCSToDCS(2,2), rotMatrixPCSToDCS(2,3));
    fprintf('                   [Tra]   [%10.5f %10.5f %10.5f]\n', rotMatrixPCSToDCS(3,1), rotMatrixPCSToDCS(3,2), rotMatrixPCSToDCS(3,3));
    fprintf('-----------------------------------------------------------------\n');
    fprintf('                   [Sag]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToDCS(1,1), rotMatrixGCSToDCS(1,2), rotMatrixGCSToDCS(1,3));
    fprintf('rotMatrixGCSToDCS: [Cor] = [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToDCS(2,1), rotMatrixGCSToDCS(2,2), rotMatrixGCSToDCS(2,3));
    fprintf('                   [Tra]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToDCS(3,1), rotMatrixGCSToDCS(3,2), rotMatrixGCSToDCS(3,3));
    fprintf('-----------------------------------------------------------------\n');
    fprintf('                   [ x ]   %10.5f [mm]\n', DCS_offset(1) * 1e3);
    fprintf('slice offset(DCS): [ y ] = %10.5f [mm]\n', DCS_offset(2) * 1e3);
    fprintf('                   [ z ]   %10.5f [mm]\n', DCS_offset(3) * 1e3);
    fprintf('=================================================================\n');

    %% Calculate nominal gradient waveforms in GCS [G/cm]
    %----------------------------------------------------------------------
    % Calculate one spiral interleaf: k in [cycle/cm], g in [G/cm]
    %----------------------------------------------------------------------
    Fcoeff = [FOV -FOV * (1 - user_opts.vds_factor / 100)]; % FOV decreases linearly from 30 to 15cm
    krmax  = 1 / (2 * (FOV / matrix_size(1)));              % [cycle/cm]
    res    = 1 / (2 * krmax) * 10;                          % resolution [mm]
    [k_spiral_arm,g_spiral_arm,s_spiral_arm,time] = vdsmex(interleaves, Fcoeff, res, Gmax, Smax, sampling_time, 10000000);
    g_spiral_arm = -g_spiral_arm(1:Nk,:); % Nk x 2

    %----------------------------------------------------------------------
    % Rotate the spiral interleaf by 360/Ni degrees every TR
    % (Re{g_spiral} + 1j * Im{g_spiral}) * (cos(arg) - 1j * sin(arg))
    % RO: real =>  Re{g_spiral} * cos(arg) + Im{g_spiral} * sin(arg)
    % PE: imag => -Re{g_spiral} * sin(arg) + Im{g_spiral} * cos(arg)
    %----------------------------------------------------------------------
    gu_nominal = zeros(Nk, Ni, 'double'); % Nk x Ni [G/cm] PE
    gv_nominal = zeros(Nk, Ni, 'double'); % Nk x Ni [G/cm] RO
    gw_nominal = zeros(Nk, Ni, 'double'); % Nk x Ni [G/cm] SL
    for i = 1:Ni
        arg = 2 * pi / Ni * (i - 1); % [rad]
        gv_nominal(:,i) =  g_spiral_arm(:,1) * cos(arg) + g_spiral_arm(:,2) * sin(arg);
        gu_nominal(:,i) = -g_spiral_arm(:,1) * sin(arg) + g_spiral_arm(:,2) * cos(arg);
    end

    %% Calculate GIRF-predicted gradient waveforms in GCS [G/cm]
    g_nominal = cat(3, gu_nominal, gv_nominal);
    tRR = 0; % custom clock-shift
    sR.R = rotMatrixGCSToDCS;
    sR.T = header.acquisitionSystemInformation.systemFieldStrength_T;
    [~,g_predicted] = apply_GIRF(g_nominal, dt, sR, tRR); % k:[cycle/cm] and g:[G/cm]

    %% Change the sign of GIRF-predicted gradient waveforms in GCS [G/cm] [PE,RO,SL]
    g_gcs = zeros(Nk, 3, Ni, 'double');
    g_gcs(:,1,:) = PE_sign * g_predicted(:,:,1); % [G/cm] PE (gu)
    g_gcs(:,2,:) = RO_sign * g_predicted(:,:,2); % [G/cm] RO (gv)
    g_gcs(:,3,:) = g_predicted(:,:,3);           % [G/cm] SL (gw)

    %% Calculate GIRF-predicted gradient waveforms in DCS [G/cm] [x,y,z]
    g_dcs = zeros(Nk, 3, Ni, 'double');
    for i = 1:Ni
        g_dcs(:,:,i) = (rotMatrixGCSToDCS * g_gcs(:,:,i).').'; % Nk x 3
    end

    %% Calculate GIRF-predicted k-space trajectories in GCS [rad/m] [PE,RO,SL]
    %----------------------------------------------------------------------
    % Numerically integrate the coefficients
    % [rad/sec/T] * [G/cm] * [T/1e4G] * [1e2cm/m] * [sec] => [rad/m]
    %----------------------------------------------------------------------
    k_gcs = cumsum(gamma * g_gcs * 1e-2 * double(dt)); % [rad/m]

    %% Calculate GIRF-predicted k-space trajectories in DCS [rad/m] [x,y,z]
    k_dcs = zeros(Nk, 3, Ni, 'double');
    for i = 1:Ni
        k_dcs(:,:,i) = (rotMatrixGCSToDCS * k_gcs(:,:,i).').'; % Nk x 3
    end

    %% Calculate GIRF-predicted k-space trajectories in RCS [rad/m] [R,C,S]
    k_rcs = zeros(Nk, 3, Ni, 'double');
    for i = 1:Ni
        k_rcs(:,:,i) = (rotMatrixRCSToGCS.' * k_gcs(:,:,i).').'; % Nk x 3
    end

    %% Calculate concomitant field basis functions (N x Nl) [m], [m^2], [m^3]
    tstart = tic; fprintf('Calculating concomitant field basis functions... ');
    p = calculate_concomitant_field_basis(x, y, z, Nl);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Calculate the time courses of phase coefficients (Nk x Nl x Ni) [rad/m], [rad/m^2], [rad/m^3]
    tstart = tic; fprintf('Calculating the time courses of phase coefficients... ');
    k = calculate_concomitant_field_coefficients(reshape(g_dcs(:,1,:), [Nk Ni]), reshape(g_dcs(:,2,:), [Nk Ni]), reshape(g_dcs(:,3,:), [Nk Ni]), Nl, B0, gamma, dt);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Set concomitant fields zero
    if concomitant_field_correction == 0
        p(:,4:end) = 0;
        k(:,4:end) = 0;
    end

    %% Calculate a density compensation function
    tstart = tic; fprintf('Cacluating a density compensation function based on sdc3_MAT.c... ');
    % [rad/m] / ([cycle/cm] * [2pi rad/cycle] * [1e2cm/m]) => [unitless]
    coords = permute(k_rcs, [2 1 3]) / (2 * krmax * 2 * pi * 1e2); % Nk x 3 x Ni => 3 x Nk x Ni
    coords(3,:,:) = 0;
    numIter = 25;
    effMtx  = matrix_size(1);
    verbose = 0;
    DCF = sdc3_MAT(coords, numIter, effMtx, verbose, 2.1);
    w = DCF / max(DCF(:));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Calculate a circular mask (image support constraint)
    switch main_orientation % (S/C/T = 0/1/2)
        case SAGITTAL
            mask = reshape((y.^2 + z.^2 < (field_of_view_mm(1) * osf * 1e-3 / 2)^2), [N1 N2]); % [mm] * [m/1e3mm] => [m]
        case CORONAL
            mask = reshape((x.^2 + z.^2 < (field_of_view_mm(1) * osf * 1e-3 / 2)^2), [N1 N2]); % [mm] * [m/1e3mm] => [m]
        case TRANSVERSE
            mask = reshape((x.^2 + y.^2 < (field_of_view_mm(1) * osf * 1e-3 / 2)^2), [N1 N2]); % [mm] * [m/1e3mm] => [m]
        otherwise
    end

    %% Demodulate k-space data (Nk x Ni x Nc)
    for i = 1:Ni
        kspace(:,i,:) = bsxfun(@times, kspace(:,i,:), exp(-1j * k_dcs(:,:,i) * DCS_offset)); % (Nk x 3) * (3 x 1) => Nk x 1
    end

    %% Perform NUFFT reconstruction
    % scaled to [-0.5,0.5] and then [-pi,pi]
    % [rad/m] / ([cycle/cm] * [2pi rad/cycle] * [1e2cm/m]) => [rad/m] / [rad/m] = [unitless] ([-0.5,0.5]) * 2pi => [-pi,pi]
    % The definition of FFT is opposite in NUFFT
    om = -cat(2, reshape(k_rcs(:,1,:), [Nk*Ni 1]), reshape(k_rcs(:,2,:), [Nk*Ni 1])) / (2 * krmax * 1e2); % Nk*Ni x 2
    Nd = [N1 N2]; % matrix size
    Jd = [6 6];   % kernel size
    Kd = Nd * 2;  % oversampled matrix size
    nufft_st = nufft_init(om, Nd, Jd, Kd, Nd/2, 'minmax:kb');

    imc_nufft = complex(zeros(N1, N2, Nc, 'double'));
    for c = 1:Nc
        tstart = tic; fprintf('(%2d/%2d): NUFFT reconstruction (c=%2d/%2d)... ', idx, nr_recons, c, Nc);
        %------------------------------------------------------------------
        % Apply the adjoint of 2D NUFFT
        %------------------------------------------------------------------
        imc_nufft(:,:,c) = nufft_adj(reshape(double(kspace(:,:,c)) .* w, [Nk*Ni 1]), nufft_st) / sqrt(prod(Nd));
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %% Calculate coil sensitivity maps
    tstart = tic; fprintf('(%2d/%2d): Calculating coil sensitivity maps with Walsh method... ', idx, nr_recons);
    %----------------------------------------------------------------------
    % IFFT to k-space (k-space <=> image-space)
    %----------------------------------------------------------------------
    kspace_gridded = ifft2c(imc_nufft);

    %----------------------------------------------------------------------
    % Calculate the calibration region of k-space
    %----------------------------------------------------------------------
    cal_shape = [32 32];
    cal_data = crop(kspace_gridded, [cal_shape Nc]);
    cal_data = bsxfun(@times, cal_data, hamming(cal_shape(1)) * hamming(cal_shape(2)).');

    %----------------------------------------------------------------------
    % Calculate coil sensitivity maps
    %----------------------------------------------------------------------
    cal_im = fft2c(zpad(cal_data, [N1 N2 Nc]));
    csm = ismrm_estimate_csm_walsh(cal_im);

    %----------------------------------------------------------------------
    % Apply a circular mask to CSM
    %----------------------------------------------------------------------
    csm = bsxfun(@times, csm, mask);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Perform optimal coil combination for NUFFT reconstruction
    im_nufft = sum(imc_nufft .* conj(csm), 3);
    
    %% Perform King's method for concomitant field correction
    start_time_king = tic; fprintf('(%2d/%2d): Performing King''s method...\n', idx, nr_recons);
    [im_king,im_fs] = perform_deblurring_king_method(double(kspace), nufft_st, w, csm, ...
	   reshape(g_dcs(:,1,:), [Nk Ni]), reshape(g_dcs(:,2,:), [Nk Ni]), reshape(g_dcs(:,3,:), [Nk Ni]), ...
	   x, y, z, rotMatrixRCSToGCS, rotMatrixGCSToPCS, rotMatrixPCSToDCS, field_of_view_mm, DCS_offset, gamma, B0, dt);
    computation_time_king = toc(start_time_king);
    fprintf('done! (%6.4f/%6.4f sec)\n', computation_time_king, toc(start_time));

    %% Perform King's method with a standard B0 correction
    if main_orientation == TRANSVERSE
        tstart = tic; fprintf('(%2d/%2d): Performing King''s method with a standard B0 correction...\n', idx, nr_recons);
        [im_king_with_B0_correction,im_fs] = perform_deblurring_king_method_with_B0_correction(double(kspace), nufft_st, w, csm, ...
            reshape(g_dcs(:,1,:), [Nk Ni]), reshape(g_dcs(:,2,:), [Nk Ni]), reshape(g_dcs(:,3,:), [Nk Ni]), ...
            x, y, z, rotMatrixRCSToGCS, rotMatrixGCSToPCS, rotMatrixPCSToDCS, field_of_view_mm, DCS_offset, gamma, B0, dt, B0map(:,:,actual_slice_nr));
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %% Perform CG-SENSE reconstruction
    tstart = tic; fprintf('(%2d/%2d): Performing CG-SENSE reconstruction...\n', idx, nr_recons);
    b = reshape(sum(conj(csm).* imc_nufft, 3), [N 1]);
    max_iterations = 30;
    limit = 1e-5;
    E = @(x,tr) encoding_sense(x, csm, w, nufft_st, tr);
    [m_sense, flag, relres, iter, resvec, lsvec] = lsqr(E, b, limit, max_iterations, [], [], []); % N x 1
    im_sense = reshape(m_sense, [N1 N2]);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Estimate a mask for limited spatial support using ESPIRiT
    if support_constraint
        tstart = tic; fprintf('(%2d/%2d): Calculating coil sensitivity maps with ESPIRiT... ', idx, nr_recons);
        %------------------------------------------------------------------
        % Set parameters
        %------------------------------------------------------------------
        ncalib = 24;    % use 24 calibration lines to compute compression
        ksize  = [6 6]; % kernel size

        % Threshold for picking singular vercors of the calibration matrix
        % (relative to largest singlular value.
        eigThresh_1 = 0.02;

        % threshold of eigen vector decomposition in image space.
        eigThresh_2 = 0.95;

        %------------------------------------------------------------------
        % Calculate "clean" k-space data
        %------------------------------------------------------------------
        kspace_sense = fft2c(bsxfun(@times, csm, im_sense));

        %------------------------------------------------------------------
        % crop a calibration area
        %------------------------------------------------------------------
        calib = crop(kspace_sense, [ncalib ncalib Nc]);

        %------------------------------------------------------------------
        % Compute ESPIRiT EigenVectors
        % We perform calibration in k-space followed by an eigen-decomposition
        % in image space to produce the EigenMaps.
        %------------------------------------------------------------------
        % compute the calibration matrix, perform SVD, and convert singular
        % vectors into k-space kernels
        [kernel,S] = dat2Kernel(calib, ksize);
        index = max(find(S >= S(1) * eigThresh_1));

        %------------------------------------------------------------------
        % Compute eigen-value decomposition in image space to get CSM maps
        %------------------------------------------------------------------
        [V,W] = kernelEig(kernel(:,:,:,1:index), [N1 N2]);

        %------------------------------------------------------------------
        % Calculate a "support" mask based on an eigenvalue threshold
        %------------------------------------------------------------------
        mask_support = (W(:,:,end) > eigThresh_2);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

        %% Apply a support mask to CSM
        csm = bsxfun(@times, csm, mask_support);
    else
        mask_support = true(N1, N2);
    end

    %% Calculate a time vector [sec]
    t = (0:Nk-1).' * dt; % Nk x 1 [sec]

    %% Calculate the SVD of the higher-order encoding matrix (Nk x N)
    start_time_svd = tic;
    os = 5; % oversampling parameter for randomized SVD
    u = complex(zeros(Nk, Lmax, Ni, 'double'));
    v = complex(zeros(N , Lmax, Ni, 'double'));
    s = complex(zeros(Lmax, Lmax, Ni, 'double'));
    for i = 1:Ni
        tstart = tic; fprintf('(%2d/%2d): Calculating the randomized SVD (i=%2d/%2d)... \n', idx, nr_recons, i, Ni);
        [U,S,V] = calculate_rsvd_higher_order_encoding_matrix(k(:,4:end,i), p(:,4:end), Lmax, os, vec(B0map(:,:,actual_slice_nr)), t, static_B0_correction);
        u(:,:,i) = U(:,1:Lmax); % Nk x Lmax
        v(:,:,i) = V(:,1:Lmax) * S(1:Lmax,1:Lmax)'; % N x Lmax
        s(:,:,i) = S(1:Lmax,1:Lmax);
        fprintf('(%2d/%2d): Calculating the randomized SVD (i=%2d/%2d)... done! (%6.4f/%6.4f sec)\n', idx, nr_recons, i, Ni, toc(tstart), toc(start_time));
    end
    computation_time_svd = toc(start_time_svd);

    %% Calculate NUFFT structures for lowrank MaxGIRF reconstruction
    st = cell(Ni,1);
    for i = 1:Ni
        tstart = tic; fprintf('(%2d/%2d): Calculating NUFFT structure... ', i, Ni);
        % scaled to [-0.5,0.5] and then [-pi,pi]
        % [rad/m] / ([cycle/cm] * [2pi rad/cycle] * [1e2cm/m]) => [rad/m] / [rad/m] = [unitless]
        om = -cat(2, k_rcs(:,1,i), k_rcs(:,2,i)) / (2 * krmax * 1e2); % Nk x 2
        Nd = [N1 N2];   % matrix size
        Jd = [6 6];     % kernel size
        Kd = Nd * 2;    % oversampled matrix size
        st{i} = nufft_init(om, Nd, Jd, Kd, Nd/2, 'minmax:kb');
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %% Perform conjugate phase reconstruction using the MaxGIRF encoding model
    start_time_cpr = tic;
    b = complex(zeros(N, Lmax, 'double'));
    for i = 1:Ni
        Nd = st{i}.Nd;
        tstart = tic; fprintf('(%d/%d): Calculating conjugate phase reconstruction (i=%d/%d)... ', idx, nr_recons, i, Ni);

        for c = 1:Nc
            %--------------------------------------------------------------
            % Caclulate d_{i,c}
            %--------------------------------------------------------------
            d = double(kspace(:,i,c)); % kspace: Nk x Ni x Nc

            %--------------------------------------------------------------
            % Calculate sum_{ell=1}^L ...
            % diag(v_tilde(ell,i)) * Fc^H * diag(conj(u_tilde(ell,i)))
            %--------------------------------------------------------------
            AHd  = zeros(N, Lmax, 'double');
            AHd_ = zeros(N, 1, 'double');
            for ell = 1:Lmax
                % Preconditioning with density compensation
                FHDuHd = nufft_adj((conj(u(:,ell,i)) .* d) .* w(:,i), st{i}) / sqrt(prod(Nd));
                AHd_ = AHd_ + v(:,ell,i) .* reshape(FHDuHd, [N 1]);
                AHd(:,ell) = AHd_;
            end

            %--------------------------------------------------------------
            % Calculate Sc^H * Ei^H * d_{i,c}
            %--------------------------------------------------------------
            for ell = 1:Lmax
                AHd(:,ell) = reshape(conj(double(csm(:,:,c))), [N 1]) .* AHd(:,ell);
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
    im_maxgirf_cpr = reshape(b, [N1 N2 Lmax]);

    %% Perform lowrank MaxGIRF reconstruction
    start_time_maxgirf = tic; fprintf('(%d/%d): Performing lowrank MaxGIRF reconstruction...\n', idx, nr_recons);
    max_iterations = 15;
    limit = 1e-5;
    E = @(x,tr) encoding_lowrank_maxgirf(x, csm, u(:,1:L,:), v(:,1:L,:), w, st, tr);
    [m_lowrank, flag, relres, iter, resvec, lsvec] = lsqr(E, b(:,L), limit, max_iterations, [], [], []); % NL x 1
    im_maxgirf_lowrank = reshape(m_lowrank, [N1 N2]);
    computation_time_maxgirf = toc(start_time_maxgirf);
    fprintf('done! (%6.4f/%6.4f sec)\n', computation_time_maxgirf, toc(start_time));
    
    %% Save each reconstruction in a mat file
    %----------------------------------------------------------------------
    % Save k and p
    %----------------------------------------------------------------------
    output_fullpath = fullfile(output_directory, sprintf('KP_slice%d', idx));
    tstart = tic; fprintf('Saving results: %s... ', output_fullpath);
    save(output_fullpath, 'k', 'p', '-v7.3');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Save NUFFT reconstruction
    %----------------------------------------------------------------------
    output_fullpath = fullfile(output_directory, sprintf('nufft_slice%d', idx));
    tstart = tic; fprintf('Saving results: %s... ', output_fullpath);
    save(output_fullpath, 'im_nufft', 'slice_nr', 'actual_slice_nr', '-v7.3');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Save King's method without a B0 correction
    %----------------------------------------------------------------------
    output_fullpath = fullfile(output_directory, sprintf('king_slice%d', idx));
    tstart = tic; fprintf('Saving results: %s... ', output_fullpath);
    save(output_fullpath, 'im_king', 'slice_nr', 'actual_slice_nr', '-v7.3');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Save King's method with a B0 correction
    %----------------------------------------------------------------------
    if main_orientation == TRANSVERSE
        output_fullpath = fullfile(output_directory, sprintf('king_with_B0_correction_slice%d', idx));
        tstart = tic; fprintf('Saving results: %s... ', output_fullpath);
        save(output_fullpath, 'im_king_with_B0_correction', 'slice_nr', 'actual_slice_nr');
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %----------------------------------------------------------------------
    % Save conjugate phase reconstruction (MaxGIRF)
    %----------------------------------------------------------------------
    output_fullpath = fullfile(output_directory, sprintf('maxgirf_cpr_slice%d', idx));
    tstart = tic; fprintf('Saving results: %s... ', output_fullpath);
    save(output_fullpath, 'im_maxgirf_cpr', 'slice_nr', 'actual_slice_nr', 'computation_time_cpr', 's', 'computation_time_svd', 'Lmax', 'L', '-v7.3');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Save iterative MaxGIRF reconstruction with a lowrank approximation
    %----------------------------------------------------------------------
    output_fullpath = fullfile(output_directory, sprintf('maxgirf_lowrank_slice%d', idx));
    tstart = tic; fprintf('Saving results: %s... ', output_fullpath);
    save(output_fullpath, 'im_maxgirf_lowrank', 'computation_time_maxgirf', 'slice_nr', 'actual_slice_nr', '-v7.3');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end
