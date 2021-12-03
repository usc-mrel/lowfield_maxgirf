% demo_2d_spiral_simulation.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 09/24/2020, Last modified: 09/25/2021

%% Clean slate
%close all; clear all; clc;

%% Initialize a random number generator
rng('default');

%% Set directory names
computer_type = computer;
if strcmp(computer_type, 'PCWIN64')
    src_directory = 'E:\lowfield_maxgirf';
    mida_directory = 'E:\lowfield_maxgirf\thirdparty\MIDA_v1.0\MIDA_v1_voxels';
    mida_filename = 'MIDA_v1';
elseif strcmp(computer_type, 'GLNXA64')
    src_directory = '/server/home/nlee/lowfield_maxgirf';
    mida_directory = '/server/home/nlee/lowfield_maxgirf/thirdparty/MIDA_v1.0/MIDA_v1_voxels';
    mida_filename = 'MIDA_v1';
end

%% Add paths
addpath(genpath(src_directory));

%% Define constants
gamma = 4257.59 * (1e4 * 2 * pi); % gyromagnetic ratio for 1H [rad/sec/T]

SAGITTAL   = 0; % Patient axis perpendicular to the sagittal plane
CORONAL    = 1; % Patient axis perpendicular to the coroal plane
TRANSVERSE = 2; % Patient axis perpendicular to the transvers plane

%% Define imaging parameters
psdname   = 'FISP';     % pulse sequence type: FISP or FLASH
SNR       = inf;        % signal-to-noise ratio, SNR = s / sigma, sigma = the standard deviation of the noise
FOV1      = 24;         % field of view in the RO dimension (row)   [cm]
FOV2      = 24;         % field of view in the PE dimension (col)   [cm]
FOV3      = 24;         % field of view in the SL dimension (slice) [cm]
osf       = 1;          % reconstruction oversampling factor
N1        = 256;        % number of voxels in the RO dimension (row)
N2        = 256;        % number of voxels in the PE dimension (column)
N3        = 256;        % number of voxels in the SL dimension (slice)
N         = N1 * N2;    % number of voxels
Nl        = 19;         % number of comcomitant field basis

field_of_view_mm = [FOV1 FOV2 FOV3] * 1e1; % [cm] * [10mm/cm] => [mm]

% image_ori  = 'sagittal'; % orientation of a slice
% B0        = 0.55;        % main field strength [T]
% Lmax      = 80;          % maximum rank of the SVD approximation of a higher-order encoding matrix
% L         = 50;          % rank of the SVD approximation of a higher-order encoding matrix
% offset    = 0;           % offset from isocenter [mm]

%% Set the slice number of interest (2D single slice at a time)
switch image_ori
    case 'sagittal'
        slice_nr = 135;
    case 'coronal'
        slice_nr = 143;
    case 'axial'
        slice_nr = 174;
    otherwise
end

%% Define MR parameters (based on Brainweb acquired at 1.5T)
start_time = tic;
nr_classes = 13; % number of tissue labels for Brainweb
MR_parameters = repmat(struct('tissue_name' , [], ...
                              'tissue_label', [], ...
                              'T1'          , [], ... % spin-lattice relaxation time [sec]
                              'T2'          , [], ... % spin-spin relaxation time    [sec]
                              'T2s'         , [], ... % T2 star relaxation time      [sec]
                              'PD'          , []), ... 
                              [nr_classes 1]);

idx = 1;
MR_parameters(idx).tissue_name  = 'BACKGROUND';
MR_parameters(idx).tissue_label = 0;
MR_parameters(idx).T1           = 0;
MR_parameters(idx).T2           = 0;
MR_parameters(idx).T2s          = 0;
MR_parameters(idx).PD           = 0;

idx = idx + 1;
MR_parameters(idx).tissue_name  = 'CSF';
MR_parameters(idx).tissue_label = 1;
MR_parameters(idx).T1           = 2569 * 1e-3; % 2569
MR_parameters(idx).T2           = 329 * 1e-3; % 329
MR_parameters(idx).T2s          = 58 * 1e-3;
MR_parameters(idx).PD           = 1;

idx = idx + 1;
MR_parameters(idx).tissue_name  = 'GREY MATTER';
MR_parameters(idx).tissue_label = 2;
MR_parameters(idx).T1           = 833 * 1e-3;
MR_parameters(idx).T2           = 83 * 1e-3;
MR_parameters(idx).T2s          = 69 * 1e-3;
MR_parameters(idx).PD           = 0.86;

idx = idx + 1;
MR_parameters(idx).tissue_name  = 'WHITE MATTER';
MR_parameters(idx).tissue_label = 3;
MR_parameters(idx).T1           = 500 * 1e-3;
MR_parameters(idx).T2           = 70 * 1e-3;
MR_parameters(idx).T2s          = 61 * 1e-3;
MR_parameters(idx).PD           = 0.77;

idx = idx + 1;
MR_parameters(idx).tissue_name  = 'FAT';
MR_parameters(idx).tissue_label = 4;
MR_parameters(idx).T1           = 350 * 1e-3;
MR_parameters(idx).T2           = 70 * 1e-3;
MR_parameters(idx).T2s          = 58 * 1e-3;
MR_parameters(idx).PD           = 1;

idx = idx + 1;
MR_parameters(idx).tissue_name  = 'MUSCLE / SKIN';
MR_parameters(idx).tissue_label = 5;
MR_parameters(idx).T1           = 900 * 1e-3;
MR_parameters(idx).T2           = 47 * 1e-3;
MR_parameters(idx).T2s          = 30 * 1e-3;
MR_parameters(idx).PD           = 1;

idx = idx + 1;
MR_parameters(idx).tissue_name  = 'SKIN';
MR_parameters(idx).tissue_label = 6;
MR_parameters(idx).T1           = 2569 * 1e-3;
MR_parameters(idx).T2           = 329 * 1e-3;
MR_parameters(idx).T2s          = 58 * 1e-3;
MR_parameters(idx).PD           = 1;

idx = idx + 1;
MR_parameters(idx).tissue_name  = 'SKULL';
MR_parameters(idx).tissue_label = 7;
MR_parameters(idx).T1           = 0;
MR_parameters(idx).T2           = 0;
MR_parameters(idx).T2s          = 0;
MR_parameters(idx).PD           = 0;

idx = idx + 1;
MR_parameters(idx).tissue_name  = 'Glial MATTER';
MR_parameters(idx).tissue_label = 8;
MR_parameters(idx).T1           = 833 * 1e-3;
MR_parameters(idx).T2           = 83 * 1e-3;
MR_parameters(idx).T2s          = 69 * 1e-3;
MR_parameters(idx).PD           = 0.86;

idx = idx + 1;
MR_parameters(idx).tissue_name  = 'MEAT';
MR_parameters(idx).tissue_label = 9;
MR_parameters(idx).T1           = 500 * 1e-3;
MR_parameters(idx).T2           = 70 * 1e-3;
MR_parameters(idx).T2s          = 61 * 1e-3;
MR_parameters(idx).PD           = 0.77;

idx = idx + 1;
MR_parameters(idx).tissue_name  = 'Blood';
MR_parameters(idx).tissue_label = {};
MR_parameters(idx).T1           = 1387.5 * 1e-3;
MR_parameters(idx).T2           = 308.5 * 1e-3;
MR_parameters(idx).T2s          = 308.5 * 1e-3;
MR_parameters(idx).PD           = 1;

idx = idx + 1;
MR_parameters(idx).tissue_name  = 'Bone Marrow (Red)';
MR_parameters(idx).tissue_label = {};
MR_parameters(idx).T1           = 549.0 * 1e-3;
MR_parameters(idx).T2           = 49 * 1e-3;
MR_parameters(idx).T2s          = 49 * 1e-3;
MR_parameters(idx).PD           = 1;

idx = idx + 1;
MR_parameters(idx).tissue_name  = 'Cartilage';
MR_parameters(idx).tissue_label = {};
MR_parameters(idx).T1           = 1045.5 * 1e-3;
MR_parameters(idx).T2           = 37.3 * 1e-3;
MR_parameters(idx).T2s          = 37.3 * 1e-3;
MR_parameters(idx).PD           = 1;

%% Define tissue classes for the MIDA brain model
%--------------------------------------------------------------------------
% BACKGROUND: 0 in Brainweb
%--------------------------------------------------------------------------
tissue_background   = {'Air Internal - Ethmoidal Sinus' ; ...
                       'Air Internal - Frontal Sinus'   ; ...
                       'Air Internal - Maxillary Sinus' ; ...
                       'Air Internal - Sphenoidal Sinus'; ...
                       'Air Internal - Mastoid'         ; ...
                       'Air Internal - Nasal/Pharynx'   ; ...
                       'Background'                     ; ...
                       'Eye Lens'                       ; ...
                       'Ear Auditory Canal'             ; ...
                       'Air Internal - Oral Cavity'     ; ...
                      };

%--------------------------------------------------------------------------
% CSF: 1 in Brainweb
%--------------------------------------------------------------------------
tissue_csf          = {'CSF Ventricles'; ...
                       'CSF General'   ; ...
                       'Eye Vitreous'  ; ...
                       'Eye Cornea'    ; ...
                       'Eye Aqueous'   ; ...
                      };

%--------------------------------------------------------------------------
% GREY MATTER: 2 in Brainweb
%--------------------------------------------------------------------------
tissue_grey_matter  = {'Cerebellum Gray Matter'       ; ...
                       'Pineal Body'                  ; ...
                       'Amygdala'                     ; ...
                       'Hippocampus'                  ; ...
                       'Caudate Nucleus'              ; ...
                       'Putamen'                      ; ...
                       'Brain Gray Matter'            ; ...
                       'Nucleus Accumbens'            ; ...
                       'Globus Pallidus'              ; ...
                       'Hypophysis or Pituitary Gland'; ...
                       'Mammillary Body'              ; ...
                       'Hypothalamus'                 ; ...
                       'Ear Cochlea'                  ; ...
                       'Ear Semicircular Canals'      ; ...
                       'Submandibular Gland'          ; ...
                       'Parotid Gland'                ; ...
                       'Sublingual Gland'             ; ...
                       'Substantia Nigra'             ; ...
                       'Thalamus'                     ; ...
                      }; 

%--------------------------------------------------------------------------
% WHITE MATTER: 3 in Brainweb
%--------------------------------------------------------------------------
tissue_white_matter = {'Cerebellum White Matter'               ; ...
                       'Brainstem Midbrain'                    ; ...
                       'Brain White Matter'                    ; ...
                       'Spinal Cord'                           ; ...
                       'Brainstem Pons'                        ; ...
                       'Brainstem Medulla'                     ; ...
                       'Optic Tract'                           ; ...
                       'Commissura (Anterior)'                 ; ...
                       'Commissura (Posterior)'                ; ...
                       'Eye Retina/Choroid/Sclera'             ; ...
                       'Cerebral Peduncles'                    ; ...
                       'Optic Chiasm'                          ; ...
                       'Cranial Nerve I - Olfactory'           ; ...
                       'Cranial Nerve II - Optic'              ; ...
                       'Cranial Nerve III - Oculomotor'        ; ...
                       'Cranial Nerve IV - Trochlear'          ; ...
                       'Cranial Nerve V - Trigeminal'          ; ...
                       'Cranial Nerve V2 - Maxillary Division' ; ...
                       'Cranial Nerve V3 - Mandibular Division'; ...
                       'Cranial Nerve VI - Abducens'           ; ...
                       'Cranial Nerve VII - Facial'            ; ...
                       'Cranial Nerve VIII - Vestibulocochlear'; ...
                       'Cranial Nerve IX - Glossopharyngeal'   ; ...
                       'Cranial Nerve X - Vagus'               ; ...
                       'Cranial Nerve XI - Accessory'          ; ...
                       'Cranial Nerve XII - Hypoglossal'       ; ...
                      };

%--------------------------------------------------------------------------
% FAT: 4 in Brainweb
%--------------------------------------------------------------------------
tissue_fat          = {'Adipose Tissue'             ; ...
                       'Subcutaneous Adipose Tissue'; ...
                      };

%--------------------------------------------------------------------------
% MUSCLE / SKIN: 5 in Brainweb
%--------------------------------------------------------------------------
tissue_muscle_skin  = {'Dura'                                         ; ...
                       'Mucosa'                                       ; ...
                       'Muscle (General)'                             ; ...
                       'Tongue'                                       ; ...
                       'Muscle - Platysma'                            ; ...
                       'Muscle - Temporalis/Temporoparietalis'        ; ...
                       'Muscle - Occipitiofrontalis - Frontal Belly'  ; ...
                       'Muscle - Lateral Pterygoid'                   ; ...
                       'Muscle - Masseter'                            ; ...
                       'Muscle - Splenius Capitis'                    ; ...
                       'Muscle - Sternocleidomastoid'                 ; ...
                       'Muscle - Occipitiofrontalis - Occipital Belly'; ...
                       'Muscle - Trapezius'                           ; ...
                       'Muscle - Mentalis'                            ; ...
                       'Muscle - Depressor Anguli Oris'               ; ...
                       'Muscle - Depressor Labii'                     ; ...
                       'Muscle - Nasalis'                             ; ...
                       'Muscle - Orbicularis Oris'                    ; ...
                       'Muscles - Procerus'                           ; ...
                       'Muscle - Levator Labii Superioris'            ; ...
                       'Muscle - Zygomaticus Major'                   ; ...
                       'Muscle - Orbicularis Oculi'                   ; ...
                       'Muscle - Levator Scapulae'                    ; ...
                       'Muscle - Medial Pterygoid'                    ; ...
                       'Muscle - Zygomaticus Minor'                   ; ...
                       'Muscles - Risorius'                           ; ...
                       'Muscle - Buccinator'                          ; ...
                       'Ear Pharyngotympanic Tube'                    ; ...
                       'Muscle - Superior Rectus'                     ; ...
                       'Muscle - Medial Rectus'                       ; ...
                       'Muscle - Lateral Rectus'                      ; ...
                       'Muscle - Inferior Rectus'                     ; ...
                       'Muscle - Superior Oblique'                    ; ...
                       'Muscle - Inferior Oblique'                    ; ...
                      };

%--------------------------------------------------------------------------
% SKIN: 6 in Brainweb
%--------------------------------------------------------------------------
tissue_skin         = {'Epidermis/Dermis'};

%--------------------------------------------------------------------------
% SKULL: 7 in Brainweb
%--------------------------------------------------------------------------
tissue_skull        = {'Mandible'             ; ...
                       'Skull'                ; ...
                       'Teeth'                ; ...
                       'Skull Inner Table'    ; ...
                       'Skull Outer Table'    ; ...
                       };

%--------------------------------------------------------------------------
% Glial MATTER: 8 in Brainweb
%--------------------------------------------------------------------------
tissue_glial_matter = {'NONE'};

%--------------------------------------------------------------------------
% MEAT: 9 in Brainweb
%--------------------------------------------------------------------------
tissue_meat         = {'Tendon - Galea Aponeurotica';...
                       'Tendon - Temporalis Tendon' ; ...
                      };

%--------------------------------------------------------------------------
% Blood
%--------------------------------------------------------------------------
tissue_blood         = {'Blood Arteries'; ...
                        'Blood Veins'   ; ...
                       };

%--------------------------------------------------------------------------
% Bone Marrow (Red)
%--------------------------------------------------------------------------
tissue_bone_marrow  = {'Vertebra - C1 (atlas)'; ...
                       'Vertebra - C2 (axis)' ; ...
                       'Vertebra - C3'        ; ...
                       'Vertebra - C4'        ; ...
                       'Vertebra - C5'        ; ...
                       'Skull Diploë'         ; ...
                       'Hyoid Bone'           ; ...
                      };

%--------------------------------------------------------------------------
% Cartilage
%--------------------------------------------------------------------------
tissue_cartilage    = {'Ear Auricular Cartilage (Pinna)'; ...
                       'Nasal Septum (Cartilage)'       ; ...
                       'Intervertebral Discs'           ; ...
                      };

MR_parameters(1).tissue_label  = tissue_background;
MR_parameters(2).tissue_label  = tissue_csf;
MR_parameters(3).tissue_label  = tissue_grey_matter;
MR_parameters(4).tissue_label  = tissue_white_matter;
MR_parameters(5).tissue_label  = tissue_fat;
MR_parameters(6).tissue_label  = tissue_muscle_skin;
MR_parameters(7).tissue_label  = tissue_skin;
MR_parameters(8).tissue_label  = tissue_skull;
MR_parameters(9).tissue_label  = tissue_glial_matter;
MR_parameters(10).tissue_label = tissue_meat;
MR_parameters(11).tissue_label = tissue_blood;
MR_parameters(12).tissue_label = tissue_bone_marrow;
MR_parameters(13).tissue_label = tissue_cartilage;

%% Read a MIDA txt file
tic; fprintf('Reading a MIDA txt file... ');
fid = fopen(fullfile(mida_directory, sprintf('%s.txt', mida_filename)));
C = textscan(fid, '%d %f %f %f %s', 'delimiter', sprintf('\t'));
fclose(fid);
fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

%% Read a MIDA brain model
% Original MIDA matrix size (number of cells)
n1 = 480;
n2 = 480;
n3 = 350;

%--------------------------------------------------------------------------
% Zeropad in image-space
%--------------------------------------------------------------------------
n = max([n1 n2 n3]);
idx1_range = (-floor(n1/2):ceil(n1/2)-1).' + floor(n/2) + 1;
idx2_range = (-floor(n2/2):ceil(n2/2)-1).' + floor(n/2) + 1;
idx3_range = (-floor(n3/2):ceil(n3/2)-1).' + floor(n/2) + 1;
mida_model = zeros(n, n, n, 'double');

%--------------------------------------------------------------------------
% Read in a MIDA head model
%--------------------------------------------------------------------------
tic; fprintf('Reading a MIDA head model... ');
fid = fopen(fullfile(mida_directory, sprintf('%s.raw', mida_filename)), 'rb');
mida_model(idx1_range, idx2_range, idx3_range) = reshape(fread(fid, 'uint8=>double'), [n1 n2 n3]); % in [0,9]
fclose(fid);
fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

%--------------------------------------------------------------------------
% Change the storage order of a MIDA brain model to match 'HFS' in PCS
%
%                                       ^ Tra (3)
%                                      /
%      R +-----------+              S +-----------+
%       /           /|               /           /|
%    L /  I     S  / |            I /  R     L  / |
%     +-----------+  |             +-----------+--|---> Sag (1)
%   A |           |  |   ====>   A |           |  |
%     |           |  |             |           |  |
%     |           | /              |           | /
%   P |           |/             P |           |/
%     +-----------+                +-----------+
%                                  |
%                                  v Cor (2)
%--------------------------------------------------------------------------
mida_model = flip(permute(mida_model, [1 3 2]), 2);

%% Calculate T1, T2, T2*, and M0 maps
T1map_highres  = zeros(n, n, n, 'double');
T2map_highres  = zeros(n, n, n, 'double');
T2smap_highres = zeros(n, n, n, 'double');
M0map_highres  = zeros(n, n, n, 'double');

nr_mida_labels = length(C{5});
for idx1 = 1:nr_classes
    tic; fprintf('(%2d/%2d) Assigning %s class... ', idx1, nr_classes, MR_parameters(idx1).tissue_name);
    %----------------------------------------------------------------------
    % Get tissue labels
    %----------------------------------------------------------------------
    tissue_label = MR_parameters(idx1).tissue_label;
    nr_labels = length(tissue_label);

    for idx2 = 1:nr_labels
        for idx3 = 1:nr_mida_labels
            if strcmp(tissue_label{idx2}, C{5}(idx3))
                %----------------------------------------------------------
                % Find the indices of voxels
                %----------------------------------------------------------
                index = find(mida_model == C{1}(idx3));

                %----------------------------------------------------------
                % Assign T1, T2, T2*, and M0
                %----------------------------------------------------------
                T1map_highres(index)  = MR_parameters(idx1).T1;  % [sec]
                T2map_highres(index)  = MR_parameters(idx1).T2;  % [sec]
                T2smap_highres(index) = MR_parameters(idx1).T2s; % [sec]
                M0map_highres(index)  = MR_parameters(idx1).PD;
            end
        end
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
end

%% Set background to zero
mask_highres = mida_model;
mask_highres(mask_highres == 50) = -50;
mask_highres = mask_highres + 50;
clear mida_model;

%% Set a slice normal vector
%--------------------------------------------------------------------------
% dNormalSag: Sagittal component of a slice normal vector (in PCS)
% dNormalCor: Coronal component of a slice normal vector (in PCS)
% dNormalTra: Transverse component of a slice normal vector (in PCS)
% dRotAngle: Slice rotation angle ("swap Fre/Pha")
%--------------------------------------------------------------------------
if ~(exist('dNormalSag', 'var') && exist('dNormalCor', 'var') && exist('dNormalTra', 'var') && exist('dRotAngle', 'var'))
    switch image_ori
        case 'sagittal'
            dNormalSag = 1;
            dNormalCor = 0;
            dNormalTra = 0;
            dRotAngle  = 0; % left-hand rule about the slice normal vector [rad]
        case 'coronal'
            dNormalSag = 0;
            dNormalCor = 1;
            dNormalTra = 0;
            dRotAngle  = 0; % left-hand rule about the slice normal vector [rad]
        case 'axial'
            dNormalSag = 0;
            dNormalCor = 0;
            dNormalTra = 1;
            dRotAngle  = 0; % left-hand rule about the slice normal vector [rad]
        otherwise
    end
end

%% Set a slice offset of a target slice in PCS
switch image_ori
    % shift the center of a scan plane by -3 cm along the "tra" direction
    case 'sagittal'
        sag_offset = offset;  % offset from isocenter [mm]
        cor_offset = 0;       % offset from isocenter [mm]
        tra_offset = -30;     % offset from isocenter [mm]
    case 'coronal'
        sag_offset = 0;       % offset from isocenter [mm]
        cor_offset = -offset; % offset from isocenter [mm]
        tra_offset = -30;     % offset from isocenter [mm]
    case 'axial'
        sag_offset = 0;       % offset from isocenter [mm]
        cor_offset = 0;       % offset from isocenter [mm]
        tra_offset = -offset; % offset from isocenter [mm]
    otherwise
end
PCS_offset = [sag_offset; cor_offset; tra_offset] * 1e-3; % [mm] * [m/1e3mm] => [m]

%% Calculate spatial coordinates in PCS
% FOV = 24 cm, voxel size = 0.0005 [m]
cor_range = (-floor(n/2):ceil(n/2)-1).' * 0.0005; % [m]
sag_range = (-floor(n/2):ceil(n/2)-1).' * 0.0005; % [m]
tra_range = (-floor(n/2):ceil(n/2)-1).' * 0.0005; % [m]
[cor,sag,tra] = ndgrid(cor_range, sag_range, tra_range);

%% Interpolate the MIDA brain model
% We are going to use a trick here. We would like to find a rotation of spatial
% coordinates within PCS but Siemens code rotates the PCS coordinates first and 
% gives the resultant coordinates in GCS. Thus, we intially apply a rotation 
% matrix from PCS to GCS and re-apply a transform matrix from GCS to PCS 
% which is computed assuming no rotation.

%--------------------------------------------------------------------------
% Calculate a rotation matrix from PCS to GCS
%--------------------------------------------------------------------------
[rotMatrixGCSToPCS,PE_sign,RO_sign,main_orientation] = calcMatrixGCSToPCS(dNormalSag, dNormalCor, dNormalTra, -dRotAngle);
rotMatrixGCSToPCS(:,1) = PE_sign * rotMatrixGCSToPCS(:,1);
rotMatrixGCSToPCS(:,2) = RO_sign * rotMatrixGCSToPCS(:,2);
rotMatrixPCSToGCS = rotMatrixGCSToPCS.';

%--------------------------------------------------------------------------
% Calculate a rotation matrix from GCS to PCS
%--------------------------------------------------------------------------
fprintf('main_orientation = %d (SAGITTAL/CORONAL/TRANSVERSE = 0/1/2)\n', main_orientation);
switch main_orientation % (S/C/T = 0/1/2)
    case SAGITTAL
        [rotMatrixGCSToPCS,PE_sign,RO_sign] = calcMatrixGCSToPCS(1, 0, 0, 0);
        rotMatrixGCSToPCS(:,1) = PE_sign * rotMatrixGCSToPCS(:,1);
        rotMatrixGCSToPCS(:,2) = RO_sign * rotMatrixGCSToPCS(:,2);
    case CORONAL
        [rotMatrixGCSToPCS,PE_sign,RO_sign] = calcMatrixGCSToPCS(0, 1, 0, 0);
        rotMatrixGCSToPCS(:,1) = PE_sign * rotMatrixGCSToPCS(:,1);
        rotMatrixGCSToPCS(:,2) = RO_sign * rotMatrixGCSToPCS(:,2);
    case TRANSVERSE
        [rotMatrixGCSToPCS,PE_sign,RO_sign] = calcMatrixGCSToPCS(0, 0, 1, 0);
        rotMatrixGCSToPCS(:,1) = PE_sign * rotMatrixGCSToPCS(:,1);
        rotMatrixGCSToPCS(:,2) = RO_sign * rotMatrixGCSToPCS(:,2);
    otherwise
end

%--------------------------------------------------------------------------
% Calculate spatial coordinates in PCS for interpolation
%--------------------------------------------------------------------------
tic; fprintf('Calculating spatial coordinates in PCS for interpolation... ');
switch main_orientation
    case SAGITTAL
        % Sag == SL (FOV3), Cor == PE (FOV2), Tra == RO (FOV1)
        matrix_size = [N3 N2 N1];
        voxel_size = [FOV3 FOV2 FOV1] * 1e-2 ./ matrix_size;
    case CORONAL
        % Sag == PE (FOV2), Cor == SL (FOV3), Tra == RO (FOV1)
        matrix_size = [N2 N3 N1];
        voxel_size = [FOV2 FOV3 FOV1] * 1e-2 ./ matrix_size;
    case TRANSVERSE
        % Sag == RO (FOV1), Cor == PE (FOV2), Tra == SL (FOV3)
        matrix_size = [N1 N2 N3];
        voxel_size = [FOV1 FOV2 FOV3] * 1e-2 ./ matrix_size;
    otherwise
end
sag_range = (-floor(matrix_size(1)/2):ceil(matrix_size(1)/2)-1).' * voxel_size(1); % [m]
cor_range = (-floor(matrix_size(2)/2):ceil(matrix_size(2)/2)-1).' * voxel_size(2); % [m]
tra_range = (-floor(matrix_size(3)/2):ceil(matrix_size(3)/2)-1).' * voxel_size(3); % [m]
[cor_,sag_,tra_] = ndgrid(cor_range, sag_range, tra_range);

pcs = rotMatrixGCSToPCS * rotMatrixPCSToGCS * cat(2, sag_(:), cor_(:), tra_(:)).';
sag_interp = reshape(pcs(1,:), matrix_size); % 3D [m]
cor_interp = reshape(pcs(2,:), matrix_size); % 3D [m]
tra_interp = reshape(pcs(3,:), matrix_size); % 3D [m]
fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

%--------------------------------------------------------------------------
% Perform 3D linear interpolation
%--------------------------------------------------------------------------
tic; fprintf('Performing 3D cubic interpolation (mask_highres)... ');
mask_pcs = interp3(sag, cor, tra, mask_highres, sag_interp, cor_interp, tra_interp, 'linear');
mask_pcs(isnan(mask_pcs)) = 0;
fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

tic; fprintf('Performing 3D cubic interpolation (T1map_highres)... ');
T1map_pcs = interp3(sag, cor, tra, T1map_highres, sag_interp, cor_interp, tra_interp, 'linear');
T1map_pcs(isnan(T1map_pcs)) = 0;
fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

tic; fprintf('Performing 3D cubic interpolation (T2map_highres)... ');
T2map_pcs = interp3(sag, cor, tra, T2map_highres, sag_interp, cor_interp, tra_interp, 'linear');
T2map_pcs(isnan(T2map_pcs)) = 0;
fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

tic; fprintf('Performing 3D cubic interpolation (T2smap_highres)... ');
T2smap_pcs = interp3(sag, cor, tra, T2smap_highres, sag_interp, cor_interp, tra_interp, 'linear');
T2smap_pcs(isnan(T2smap_pcs)) = 0;
fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

tic; fprintf('Performing 3D cubic interpolation (M0map_highres)... ');
M0map_pcs = interp3(sag, cor, tra, M0map_highres, sag_interp, cor_interp, tra_interp, 'linear');
M0map_pcs(isnan(M0map_pcs)) = 0;
fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

%--------------------------------------------------------------------------
% Clear variables
%--------------------------------------------------------------------------
clear sag cor tra sag_ cor_ tra_ sag_interp cor_interp tra_interp;
clear mask_highres T1map_highres T2map_highres T2smap_highres M0map_highres;
clear pcs;

%% Reorient the data into RCS
switch image_ori
    case 'sagittal'
        mask_rcs   = permute(mask_pcs  , [3 1 2]);
        T1map_rcs  = permute(T1map_pcs , [3 1 2]);
        T2map_rcs  = permute(T2map_pcs , [3 1 2]);
        T2smap_rcs = permute(T2smap_pcs, [3 1 2]);
        M0map_rcs  = permute(M0map_pcs , [3 1 2]);
    case 'coronal'
        mask_rcs   = flip(permute(mask_pcs  , [3 2 1]), 1);
        T1map_rcs  = flip(permute(T1map_pcs , [3 2 1]), 1);
        T2map_rcs  = flip(permute(T2map_pcs , [3 2 1]), 1);
        T2smap_rcs = flip(permute(T2smap_pcs, [3 2 1]), 1);
        M0map_rcs  = flip(permute(M0map_pcs , [3 2 1]), 1);
    case 'axial'
        mask_rcs   = flip(permute(mask_pcs  , [2 1 3]), 1);
        T1map_rcs  = flip(permute(T1map_pcs , [2 1 3]), 1);
        T2map_rcs  = flip(permute(T2map_pcs , [2 1 3]), 1);
        T2smap_rcs = flip(permute(T2smap_pcs, [2 1 3]), 1);
        M0map_rcs  = flip(permute(M0map_pcs , [2 1 3]), 1);
    otherwise
end

%--------------------------------------------------------------------------
% Select only the target slice
%--------------------------------------------------------------------------
mask_ground_truth   = mask_rcs(:,:,slice_nr);
T1map_ground_truth  = T1map_rcs(:,:,slice_nr);
T2map_ground_truth  = T2map_rcs(:,:,slice_nr);
T2smap_ground_truth = T2smap_rcs(:,:,slice_nr);
M0map_ground_truth  = M0map_rcs(:,:,slice_nr);

clear mask_rcs T1map_rcs T2map_rcs T2smap_rcs M0map_rcs;
clear mask_pcs T1map_pcs T2map_pcs T2smap_pcs M0map_pcs;

fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

%% Calculate the spatial coordinates [m]
tic; fprintf('Calculating the spatial coordinates... ');

%--------------------------------------------------------------------------
% Calculate a transformation matrix from RCS to GCS [r,c,s] <=> [PE,RO,SL]
%--------------------------------------------------------------------------
rotMatrixRCSToGCS = [0    1    0 ; % [PE]   [0 1 0] * [r]
                     1    0    0 ; % [RO] = [1 0 0] * [c]
                     0    0    1]; % [SL]   [0 0 1] * [s]

%--------------------------------------------------------------------------
% Calculate a rotation matrix from GCS to PCS
%--------------------------------------------------------------------------
[rotMatrixGCSToPCS,PE_sign,RO_sign,main_orientation] = calcMatrixGCSToPCS(dNormalSag, dNormalCor, dNormalTra, dRotAngle);
rotMatrixGCSToPCS(:,1) = PE_sign * rotMatrixGCSToPCS(:,1);
rotMatrixGCSToPCS(:,2) = RO_sign * rotMatrixGCSToPCS(:,2);

%--------------------------------------------------------------------------
% Calculate a rotation matrix from PCS to DCS
%--------------------------------------------------------------------------
rotMatrixPCSToDCS = calcMatrixPCSToDCS('HFS');

%--------------------------------------------------------------------------
% Calculate a rotation matrix from GCS to DCS
%--------------------------------------------------------------------------
rotMatrixGCSToDCS = rotMatrixPCSToDCS * rotMatrixGCSToPCS;

%--------------------------------------------------------------------------
% Calculate a scaling matrix
%--------------------------------------------------------------------------
scaling_matrix = diag(voxel_size); % [m]

%--------------------------------------------------------------------------
% Calculate a slice offset in DCS
%--------------------------------------------------------------------------
DCS_offset = rotMatrixPCSToDCS * PCS_offset; % 3 x 1

%--------------------------------------------------------------------------
% Calculate spatial coordinates in DCS
%--------------------------------------------------------------------------
idx1_range = (-floor(N1/2):ceil(N1/2)-1).';
idx2_range = (-floor(N2/2):ceil(N2/2)-1).';
idx3_range = 0;
[I1,I2,I3] = ndgrid(idx1_range, idx2_range, idx3_range);
r_dcs = (repmat(DCS_offset, [1 N]) + rotMatrixPCSToDCS * rotMatrixGCSToPCS * rotMatrixRCSToGCS * scaling_matrix * cat(2, I1(:), I2(:), I3(:)).').'; % N x 3
x = r_dcs(:,1); % N x 1 [m]
y = r_dcs(:,2); % N x 1 [m]
z = r_dcs(:,3); % N x 1 [m]

%% Calculate a tissue mask
tic; fprintf('Calculating a tissue mask... ');
mask = (mask_ground_truth > 0);
fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

%% Find voxels within a mask
voxel_index = find(mask);
nr_voxels = length(voxel_index);
[sub1,sub2] = ind2sub([N1 N2], voxel_index);

%% Calculate a circular mask (image support constraint)
switch image_ori
    case 'sagittal'
        mask_support = reshape((y.^2 + z.^2 < (FOV1 / 2)^2), [N1 N2]); % [m]
    case 'coronal'
        mask_support = reshape((x.^2 + z.^2 < (FOV1 / 2)^2), [N1 N2]); % [m]
    case 'axial'
        mask_support = reshape((x.^2 + y.^2 < (FOV1 / 2)^2), [N1 N2]); % [m]
    otherwise
end

%% Load coil sensitivity maps
tic; fprintf('Loading coil sensitivity maps... ');
load('smaps_phantom.mat');
Nc = size(smaps,3);

%--------------------------------------------------------------------------
% Resize to fill the empty space
%--------------------------------------------------------------------------
smaps_resized = imresize(smaps, [400 400]);

%--------------------------------------------------------------------------
% Crop and expand to 3D
%--------------------------------------------------------------------------
idx1_range = (-floor(N1/2):ceil(N1/2)-1).' + floor(400/2) + 1;
idx2_range = (-floor(N2/2):ceil(N2/2)-1).' + floor(400/2) + 1;
csm_ground_truth = smaps_resized(idx1_range,idx2_range,:);

%--------------------------------------------------------------------------
% Normalize CSM to remove intensity nonuniformity in reconstructed images
%--------------------------------------------------------------------------
csm_ground_truth = csm_ground_truth ./ sqrt(sum(abs(csm_ground_truth).^2,3));

%--------------------------------------------------------------------------
% Apply a support mask to CSM
%--------------------------------------------------------------------------
%csm_ground_truth = bsxfun(@times, csm_ground_truth, mask_support);

%--------------------------------------------------------------------------
% Clear variables
%--------------------------------------------------------------------------
clear smaps smaps_resized;

fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

%% Calculate a noise level based on the steady-state signal amplitude for a TrueFISP experiment
tic; fprintf('Calculating a noise standard deviation... ');

%--------------------------------------------------------------------------
% Set up imaging parameters for bSSFP
%--------------------------------------------------------------------------
TR    = 15e-3;         % repetition time [sec]
TE    = TR / 2;        % echo time [sec]
theta = 50 * pi / 180; % flip angle [rad]

%--------------------------------------------------------------------------
% Calculate a bSSFP image
%--------------------------------------------------------------------------
im_bssfp = zeros(N1, N2, 'double');
for idx = 1:nr_voxels
    %----------------------------------------------------------------------
    % Get tissue parameters at a given voxel
    %----------------------------------------------------------------------
    idx1 = sub1(idx);
    idx2 = sub2(idx);
    T1 = T1map_ground_truth(idx1,idx2);
    T2 = T2map_ground_truth(idx1,idx2);
    M0 = M0map_ground_truth(idx1,idx2);

    %----------------------------------------------------------------------
    % Calculate a steady-state bSSFP signal without MT effects
    %----------------------------------------------------------------------
    E1 = exp(-TR / T1);
    E2 = exp(-TR / T2);
    kappa = M0 * sin(theta) * exp(-TE / T2);
    im_bssfp(idx1,idx2) = kappa * (1 - E1) / (1 - E1 * E2 - (E1 - E2) * cos(theta)); % only real
end

%--------------------------------------------------------------------------
% SNR = s / sigma, where sigma denotes the standard deviation of the noise
%--------------------------------------------------------------------------
sigma = max(im_bssfp(:)) / SNR;

fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

%% Calculate a ground truth image
tic; fprintf('Calculating a ground truth %s image... ', psdname);

%--------------------------------------------------------------------------
% Set up imaging parameters for FISP
%--------------------------------------------------------------------------
if strcmp(psdname, 'FISP')
    TR    = 60e-3;         % repetition time [sec]
    TE    = 1e-3;          % echo time [sec]
    theta = 30 * pi / 180; % flip angle [rad]
elseif strcmp(psdname, 'FLASH')
    TR    = 60e-3;         % repetition time [sec]
    TE    = 1e-3;          % echo time [sec]
    theta = 40 * pi / 180; % flip angle [rad]    
end

%--------------------------------------------------------------------------
% Calculate a ground truth image
%--------------------------------------------------------------------------
im_ground_truth = zeros(N1, N2, 'double');
for idx = 1:nr_voxels
    %----------------------------------------------------------------------
    % Get tissue parameters at a given voxel
    %----------------------------------------------------------------------
    idx1 = sub1(idx);
    idx2 = sub2(idx);
    T1  = T1map_ground_truth(idx1,idx2);
    T2  = T2map_ground_truth(idx1,idx2);
    T2s = T2smap_ground_truth(idx1,idx2); % T2* relaxation time [sec]
    M0  = M0map_ground_truth(idx1,idx2);

    %----------------------------------------------------------------------
    % Calculate E1 and E2
    %----------------------------------------------------------------------
    E1 = exp(-TR / T1);
    E2 = exp(-TR / T2);

    if strcmp(psdname, 'FISP')
        %------------------------------------------------------------------
        % Calculate a steady-state FISP signal without MT effects
        %------------------------------------------------------------------
        C = E2 * (E1 - 1) * (1 + cos(theta));
        D = 1 - E1 * cos(theta) - (E1 - cos(theta)) * E2^2;
        im_ground_truth(idx1,idx2) = M0 * (1 - E1) / C * ((C + D * E2) / sqrt(D^2 - C^2) - E2) * exp(-TE / T2s); % only real
    elseif strcmp(psdname, 'FLASH')
        %------------------------------------------------------------------
        % Calculate a steady-state FLASH signal without MT effects
        %------------------------------------------------------------------
        kappa = M0 * sin(theta) * exp(-TE / T2s);
        im_ground_truth(idx1,idx2) = kappa * (1 - E1) / (1 - cos(theta) * E1); % only real
    end
end

%--------------------------------------------------------------------------
% Remove voxels with either NaN or Inf values
%--------------------------------------------------------------------------
im_ground_truth(isnan(im_ground_truth(:))) = 0;
im_ground_truth(isinf(im_ground_truth(:))) = 0;

fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

%% Sanity check!
if 0
figure;
surf(reshape(x*1e3, [N1 N2]), reshape(y*1e3, [N1 N2]), reshape(z*1e3, [N1 N2]), im_ground_truth, 'EdgeColor', 'none');
axis equal;
set(gca, 'zdir', 'reverse');
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
end

%% Calculate nominal gradient waveforms in GCS [G/cm]
% duration 9.86 ms from SPI-MRF
%--------------------------------------------------------------------------
% Calculate one spiral interleaf: k in [cycle/cm], g in [G/cm]
%--------------------------------------------------------------------------
[k_spiral_arm,g_spiral_arm,s_spiral_arm,time] = vdsmex(Ni, Fcoeff, res, Gmax, Smax, dt, 10000000);
g_spiral_arm = -g_spiral_arm; % Nk x 2

%--------------------------------------------------------------------------
% Set the dimensions
%--------------------------------------------------------------------------
Nk = length(g_spiral_arm); % number of samples in one spiral interleaf

%--------------------------------------------------------------------------
% Calculate the readout duration
%--------------------------------------------------------------------------
T = Nk * dt; % [sec]

%--------------------------------------------------------------------------
% Rotate the spiral interleaf by 360/Ni degrees every TR
% (Re{g_spiral} + 1j * Im{g_spiral}) * (cos(arg) - 1j * sin(arg))
% RO: real =>  Re{g_spiral} * cos(arg) + Im{g_spiral} * sin(arg)
% PE: imag => -Re{g_spiral} * sin(arg) + Im{g_spiral} * cos(arg)
%--------------------------------------------------------------------------
g_gcs = zeros(Nk, 3, Ni, 'double'); % Nk x 3 x Ni [G/cm] [PE,RO,SL]
for i = 1:Ni
    arg = 2 * pi / Ni * (i - 1); % [rad]
    g_gcs(:,1,i) =  g_spiral_arm(:,1) * cos(arg) + g_spiral_arm(:,2) * sin(arg);
    g_gcs(:,2,i) = -g_spiral_arm(:,1) * sin(arg) + g_spiral_arm(:,2) * cos(arg);
end

%% Calculate spatial coordinates in DCS for reconstruction
%--------------------------------------------------------------------------
% Calculate a scaling matrix
%--------------------------------------------------------------------------
scaling_matrix_recon = diag([res; res; res]) * 1e-3; % [mm] * [1e-3m/mm] => [m]

%--------------------------------------------------------------------------
% Calculate spatial coordinates in DCS
%--------------------------------------------------------------------------
idx1_range = (-floor(N1*osf/2):ceil(N1*osf/2)-1).';
idx2_range = (-floor(N2*osf/2):ceil(N2*osf/2)-1).';
idx3_range = 0;
[I1,I2,I3] = ndgrid(idx1_range, idx2_range, idx3_range);
r_dcs = (repmat(DCS_offset, [1 N*osf^2]) + rotMatrixPCSToDCS * rotMatrixGCSToPCS * rotMatrixRCSToGCS * scaling_matrix_recon * cat(2, I1(:), I2(:), I3(:)).').'; % N x 3
x_recon = r_dcs(:,1); % N * osf^2 x 1 [m]
y_recon = r_dcs(:,2); % N * osf^2 x 1 [m]
z_recon = r_dcs(:,3); % N * osf^2 x 1 [m]
fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

%% Calculate nominal gradient waveforms in DCS [G/cm] [x,y,z]
g_dcs = zeros(Nk, 3, Ni, 'double');
for i = 1:Ni
    g_dcs(:,:,i) = (rotMatrixGCSToDCS * g_gcs(:,:,i).').'; % Nk x 3
end

%% Calculate nominal k-space trajectories in GCS [rad/m] [PE,RO,SL]
%--------------------------------------------------------------------------
% Numerically integrate the coefficients
% [rad/sec/T] * [G/cm] * [T/1e4G] * [1e2cm/m] * [sec] => [rad/m]
%--------------------------------------------------------------------------
k_gcs = cumsum(gamma * g_gcs * 1e-2 * double(dt)); % [rad/m]

%% Calculate nominal k-space trajectories in DCS [rad/m] [x,y,z]
k_dcs = zeros(Nk, 3, Ni, 'double');
for i = 1:Ni
    k_dcs(:,:,i) = (rotMatrixGCSToDCS * k_gcs(:,:,i).').'; % Nk x 3
end

%% Calculate nominal k-space trajectories in RCS [rad/m] [R,C,S]
k_rcs = zeros(Nk, 3, Ni, 'double');
for i = 1:Ni
    k_rcs(:,:,i) = (rotMatrixRCSToGCS.' * k_gcs(:,:,i).').'; % Nk x 3
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

%% Calculate NUFFT structures for lowrank MaxGIRF reconstruction
%--------------------------------------------------------------------------
% Prepare an NUFFT structure using k-space trajectories in RCS
% Note: RCS <=> GCS (PRS)
% ku:[PE]   [0 1 0] * [r]
% kv:[RO] = [1 0 0] * [c]
% kw:[SL]   [0 0 1] * [s]
%--------------------------------------------------------------------------
st = cell(Ni,1);
for i = 1:Ni
    tic; fprintf('(%2d/%2d): Calculating NUFFT structure... ', i, Ni);
    % scaled to [-0.5,0.5] and then [-pi,pi]
    % [rad/m] / ([cycle/cm] * [2pi rad/cycle] * [1e2cm/m]) => [rad/m] / [rad/m] = [unitless] ([-0.5,0.5]) * 2pi => [-pi,pi]
    % The definition of FFT is opposite in NUFFT
    om = -cat(2, k_rcs(:,1,i), k_rcs(:,2,i)) / (2 * krmax * 1e2); % Nk x 2
    Nd = [N1 N2] * osf; % matrix size * oversampling factor
    Jd = [6 6];         % kernel size
    Kd = Nd * 2;        % oversampled matrix size
    st{i} = nufft_init(om, Nd, Jd, Kd, Nd/2, 'minmax:kb');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
end

%% Calculate NUFFT structure for NUFFT reconstruction
om = -cat(2, reshape(k_rcs(:,1,:), [Nk*Ni 1]), reshape(k_rcs(:,2,:), [Nk*Ni 1])) / (2 * krmax * 1e2); % Nk*Ni x 2
Nd = [N1 N2] * osf; % matrix size * oversampling factor
Jd = [6 6];         % kernel size
Kd = Nd * 2;        % oversampled matrix size
nufft_st = nufft_init(om, Nd, Jd, Kd, Nd/2, 'minmax:kb');

%% Calculate the time courses of phase coefficients (Nk x Nl x Ni) [rad/m], [rad/m^2], [rad/m^3]
tic; fprintf('Calculating the time courses of phase coefficients... ');
k = calculate_concomitant_field_coefficients(reshape(g_dcs(:,1,:), [Nk Ni]), reshape(g_dcs(:,2,:), [Nk Ni]), reshape(g_dcs(:,3,:), [Nk Ni]), Nl, B0, gamma, dt);
k = double(k);
fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

%% Calculate concomitant field basis functions (N x Nl) [m], [m^2], [m^3]
tic; fprintf('Calculating concomitant field basis functions... ');
p = calculate_concomitant_field_basis(x, y, z, Nl);
p = double(p);

p_recon = calculate_concomitant_field_basis(x_recon, y_recon, z_recon, Nl);
p_recon = double(p_recon);
fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

%% Calculate non-Cartesian k-space data with concomitant fields
kspace = complex(zeros(Nk, Ni, Nc, 'double')); % Nk x Ni x Nc

for i = 1:Ni % interleaves
    %----------------------------------------------------------------------
    % Calculate the i-th encoding matrix Ei (Nk x N)
    %----------------------------------------------------------------------
    tic; fprintf('Calculating the ith encoding matrix E (i=%d/%d)... ', i, Ni);
    phi = k(:,:,i) * p.'; % (Nk x Nl) * (N x Nl).' => Nk x N
    E = exp(1j * phi);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

    %----------------------------------------------------------------------
    % Add complex Gaussian noise to the k-space data
    % SNR = s / sigma, where sigma denotes the standard deviation of the noise
    % variance = sigma squared
    % Note: var(a * X) = a^2 var(X)
    % When noise is multiplied by E, the Gaussian random variables form a
    % linear combination with complex-valued weighting coefficients.
    % For simplicity, we assume the coefficients are 1. Then the variance of
    % a linear combination becomes
    % var(sum_i^N X_i) = sum_i^N var(X_i) = N * var(X).
    %----------------------------------------------------------------------
    noise = 1 / sqrt(2) * sigma * complex(randn(Nk, Nc, 'double'), randn(Nk, Nc, 'double'));
    kspace_noise = sqrt(N) * noise;

    tic; fprintf('Calculating non-Cartesian k-space data (i=%d/%d)... ', i, Ni);
    %----------------------------------------------------------------------
    % Calculate non-Cartesian k-space data (all channels)
    %----------------------------------------------------------------------
    kspace_interleaf = E * reshape(bsxfun(@times, im_ground_truth, csm_ground_truth), [N Nc]); % Nk x Nc

    %----------------------------------------------------------------------
    % Add complex Gaussian noise to the k-space data
    %----------------------------------------------------------------------
    kspace(:,i,:) = kspace_interleaf + kspace_noise;
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
end
clear phi E;

%% Demodulate k-space data (Nk x Ni x Nc)
for i = 1:Ni
    kspace(:,i,:) = bsxfun(@times, kspace(:,i,:), exp(-1j * k_dcs(:,:,i) * DCS_offset)); % (Nk x 3) * (3 x 1) => Nk x 1
end

%% Perform NUFFT reconstruction
imc_nufft = complex(zeros(N1*osf, N2*osf, Nc, 'double'));
for c = 1:Nc
    tic; fprintf('NUFFT reconstruction (c=%2d/%2d)... ', c, Nc);
    %----------------------------------------------------------------------
    % Apply the adjoint of 2D NUFFT
    %----------------------------------------------------------------------
    imc_nufft(:,:,c) = nufft_adj(reshape(kspace(:,:,c) .* w, [Nk*Ni 1]), nufft_st) / sqrt(prod(Nd));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
end

%% Calculate coil sensitivity maps
tstart = tic; fprintf('Calculating coil sensitivity maps with Walsh method... ');
%--------------------------------------------------------------------------
% IFFT to k-space (k-space <=> image-space)
%--------------------------------------------------------------------------
kspace_gridded = ifft2c(imc_nufft);

%--------------------------------------------------------------------------
% Calculate the calibration region of k-space
%--------------------------------------------------------------------------
cal_shape = [32 32];
cal_data = crop(kspace_gridded, [cal_shape Nc]);
cal_data = bsxfun(@times, cal_data, hamming(cal_shape(1)) * hamming(cal_shape(2)).');

%--------------------------------------------------------------------------
% Calculate coil sensitivity maps
%--------------------------------------------------------------------------
cal_im = fft2c(zpad(cal_data, [N1*osf N2*osf Nc]));
csm = ismrm_estimate_csm_walsh(cal_im);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Perform optimal coil combination for NUFFT reconstruction
im_nufft = sum(imc_nufft .* conj(csm_ground_truth), 3);

%% Sanity check!
if 1
figure;
surf(reshape(x_recon*1e3, [N1 N2]*osf), reshape(y_recon*1e3, [N1 N2]*osf), reshape(z_recon*1e3, [N1 N2]*osf), abs(im_nufft), 'EdgeColor', 'none');
axis equal; colormap(gray(256)); colorbar;
set(gca, 'zdir', 'reverse');
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
end

%
% %%
% c1 = floor(N1/2) + 1 - 10;
% c2 = floor(N2/2) + 1;
% 
% d = im_ground_truth / max(im_ground_truth(c1,c2)) - im_nufft  / max(im_nufft(c1,c2));
% figure, imagesc(abs(flip(d,1))); axis image; colormap(gray(256));
% 

%% Perform King's method for concomitant field correction
tic; fprintf('Performing King''s method...\n');
[im_king, im_fs, tc, fc, fc_XYZ] = perform_deblurring_king_method(kspace, nufft_st, w, csm_ground_truth, reshape(g_dcs(:,1,:), [Nk Ni]), reshape(g_dcs(:,2,:), [Nk Ni]), reshape(g_dcs(:,3,:), [Nk Ni]), x_recon, y_recon, z_recon, rotMatrixRCSToGCS, rotMatrixGCSToPCS, rotMatrixPCSToDCS, field_of_view_mm, DCS_offset, gamma, B0, dt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

%% Calculate the SVD of the higher-order encoding matrix (Nk x N)
start_time_svd = tic;
os = 5; % oversampling parameter for randomized SVD
u_tilde = zeros(Nk, Lmax, Ni, 'double');
v_tilde = zeros(N*osf^2, Lmax, Ni, 'double');
singular_values = zeros(Lmax, Lmax, Ni, 'double');
for i = 1:Ni
    tstart = tic; fprintf('Calculating the randomized SVD (i=%2d/%2d)... \n', i, Ni);
    [U,S,V] = calculate_rsvd_higher_order_encoding_matrix(k(:,4:end,i), p_recon(:,4:end), Lmax, os, vec(zeros(N1*osf,N2*osf, 'double')), zeros(Nk,1), 0);
    u_tilde(:,:,i) = U(:,1:Lmax); % Nk x Lmax
    v_tilde(:,:,i) = V(:,1:Lmax) * S(1:Lmax,1:Lmax)'; % N x Lmax
    singular_values(:,:,i) = S(1:Lmax,1:Lmax);
    fprintf('Calculating the randomized SVD (i=%2d/%2d)... done! (%6.4f/%6.4f sec)\n', i, Ni, toc(tstart), toc(start_time));
end
computation_time_svd = toc(start_time_svd);

%% Perform conjugate phase reconstruction using the MaxGIRF encoding model
b = zeros(N*osf^2, Lmax, 'double');
for i = 1:Ni
    Nd = st{i}.Nd;
    tic; fprintf('Calculating conjugate phase reconstruction (i=%d/%d)... ', i, Ni);

    for c = 1:Nc
        %-----------------------------------------------------------------
        % Caclulate d_{i,c}
        %-----------------------------------------------------------------
        d = double(kspace(:,i,c)); % kspace: Nk x Ni x Nc

        %-----------------------------------------------------------------
        % Calculate sum_{ell=1}^L ...
        % diag(v_tilde(ell,i)) * Fc^H * diag(conj(u_tilde(ell,i)))
        %-----------------------------------------------------------------
        AHd  = zeros(N*osf^2, Lmax, 'double');
        AHd_ = zeros(N*osf^2, 1, 'double');
        for ell = 1:Lmax
            % Preconditioning with density compensation
            FHDuHd = nufft_adj((conj(u_tilde(:,ell,i)) .* d) .* w(:,i), st{i}) / sqrt(prod(Nd));
            AHd_ = AHd_ + v_tilde(:,ell,i) .* reshape(FHDuHd, [N*osf^2 1]);
            AHd(:,ell) = AHd_;
        end

        %-----------------------------------------------------------------
        % Calculate Sc^H * Ei^H * d_{i,c}
        %-----------------------------------------------------------------
        for ell = 1:Lmax
            AHd(:,ell) = reshape(conj(double(csm_ground_truth(:,:,c))), [N*osf^2 1]) .* AHd(:,ell);
        end

        %-----------------------------------------------------------------
        % Calculate b (N x 1)
        %-----------------------------------------------------------------
        for ell = 1:Lmax
            b(:,ell) = b(:,ell) + AHd(:,ell);
        end
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
end
im_cpr = reshape(b, [N1*osf N2*osf Lmax]);

%% Perform lowrank MaxGIRF reconstruction
tstart = tic; fprintf('Performing lowrank MaxGIRF reconstruction...\n');
max_iterations = 20;
limit = 1e-5;
E = @(x,tr) encoding_lowrank_MaxGIRF(x, csm_ground_truth, u_tilde(:,1:L,:), v_tilde(:,1:L,:), w, st, tr);
[m_maxgirf, flag, relres, iter, resvec, lsvec] = lsqr(E, b(:,L), limit, max_iterations, [], [], []); % NL x 1
im_maxgirf = reshape(m_maxgirf, [N1 N2]*osf);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Create output directory
output_directory = fullfile(pwd, sprintf('%s_%s_s%d_osf%d_res%4.2f_SNR%d_%gT_G%g_S%g_T%g_Ni%d_offset%d_L%d_Lmax%d', psdname, image_ori, slice_nr, osf, res, SNR, B0, Gmax*1e1, Smax*1e-2, T*1e3, Ni, offset, L, Lmax));
mkdir(output_directory);

%% Save k and p
output_fullpath = fullfile(output_directory, 'KP');
tstart = tic; fprintf('Saving results: %s... ', output_fullpath);
save(output_fullpath, 'x_recon', 'y_recon', 'z_recon', 'k', 'p_recon', 'tc', 'fc', 'fc_XYZ', '-v7.3');
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Save ground truth
output_fullpath = fullfile(output_directory, 'ground_truth');
tic; fprintf('Saving results: %s... ', output_fullpath);
save(output_fullpath, 'im_ground_truth', '-v7.3');
fprintf('done! (%6.4f sec)\n', toc);

%% Save NUFFT reconstruction
output_fullpath = fullfile(output_directory, 'nufft');
tic; fprintf('Saving results: %s... ', output_fullpath);
save(output_fullpath, 'im_nufft', '-v7.3');
fprintf('done! (%6.4f sec)\n', toc);

%% Save King's method
output_fullpath = fullfile(output_directory, 'king');
tic; fprintf('Saving results: %s... ', output_fullpath);
save(output_fullpath, 'im_king', '-v7.3');
fprintf('done! (%6.4f sec)\n', toc);

%% Save conjugate phase reconstruction (MaxGIRF)
output_fullpath = fullfile(output_directory, 'maxgirf_cpr');
tic; fprintf('Saving results: %s... ', output_fullpath);
save(output_fullpath, 'im_cpr', 'singular_values', 'computation_time_svd', '-v7.3');
fprintf('done! (%6.4f sec)\n', toc);

%% Save iterative MaxGIRF reconstruction with a lowrank approximation
output_fullpath = fullfile(output_directory, 'maxgirf_lowrank');
tic; fprintf('Saving results: %s... ', output_fullpath);
save(output_fullpath, 'im_maxgirf', '-v7.3');
fprintf('done! (%6.4f sec)\n', toc);
