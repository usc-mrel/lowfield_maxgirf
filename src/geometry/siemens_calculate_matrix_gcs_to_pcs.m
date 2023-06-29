function [R_gcs2pcs,phase_sign,read_sign,main_orientation] = siemens_calculate_matrix_gcs_to_pcs(dNormalSag, dNormalCor, dNormalTra, dRotAngle)
% Calculate a rotation matrix from GCS to PCS
% [Sag]               [PE]
% [Cor] = R_gcs2pcs * [RO]
% [Tra]               [SS]
%
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/16/2022, Last modified: 01/16/2022

%% Define constants
SAGITTAL    = 0; % Patient axis perpendicular to the sagittal plane
CORONAL     = 1; % Patient axis perpendicular to the coroal plane
TRANSVERSAL = 2; % Patient axis perpendicular to the transvers plane

%% Calculate the signs of gradient waveforms due to ("swap Fre/Pha")
main_orientation = fGSLClassOri(dNormalSag, dNormalCor, dNormalTra);

switch main_orientation
    case SAGITTAL % =0
        if dRotAngle >= 0
            if ((dRotAngle < pi/4) && ~fGSLAlmEqual(dRotAngle, pi/4)) % 0 <= angle < 45
                read_sign = -1;
                phase_sign =  1;
            elseif (((dRotAngle > pi/4) || fGSLAlmEqual(dRotAngle, pi/4)) && (dRotAngle < 3*pi/4)) % 45 <= angle < 135
                read_sign =  1;
                phase_sign =  1;
            elseif ((dRotAngle > 3*pi/4) || fGSLAlmEqual(dRotAngle, 3*pi/4)) % 135 <= angle <= 180
                read_sign =  1;
                phase_sign = -1;
            end
        elseif dRotAngle < 0
            if ((dRotAngle > -pi/4) && ~fGSLAlmEqual(dRotAngle, -pi/4)) % -45 < angle <= 0
                read_sign = -1;
                phase_sign =  1;
            elseif ((dRotAngle > -3*pi/4) && ((dRotAngle < -pi/4) || fGSLAlmEqual(dRotAngle, -pi/4))) % -135 <= angle <= -45
                read_sign = -1;
                phase_sign = -1;
            elseif (dRotAngle < -3*pi/4) % angle < -135
                read_sign =  1;
                phase_sign = -1;
            end
        end
    case CORONAL % =1
            if ((dRotAngle >= -pi/4) && (dRotAngle <= pi/4)) % -45 <= angle <= 45
                read_sign =  1;
                phase_sign =  1;
            else % -180 <= angle < -45 or 45 < angle <= 180
                read_sign = -1;
                phase_sign = -1;
            end
    case TRANSVERSAL % =2
        if dRotAngle >= 0
            if (dRotAngle <= pi/4) % 0 <= angle <= 45
                read_sign = -1;
                phase_sign =  1;
            elseif ((dRotAngle > pi/4) && (dRotAngle < 3*pi/4) && ~fGSLAlmEqual(dRotAngle, 3*pi/4)) % 45 < angle < 135
                read_sign =  1;
                phase_sign =  1;
            elseif ((dRotAngle > 3*pi/4) || fGSLAlmEqual(dRotAngle, 3*pi/4)) % 135 <= angle <= 180
                read_sign =  1;
                phase_sign = -1;
            end
        elseif dRotAngle < 0
            if ((dRotAngle > -pi/4) || fGSLAlmEqual(dRotAngle, -pi/4)) % -45 <= angle <= 0
                read_sign = -1;
                phase_sign =  1;
            elseif ((dRotAngle > -3*pi/4) && ~fGSLAlmEqual(dRotAngle, -3*pi/4) && (dRotAngle < -pi/4)) % -135 < angle < -45
                read_sign = -1;
                phase_sign = -1;
            else % -180 <= angle < -135
                read_sign =  1;
                phase_sign = -1;
            end
        end
    otherwise
end

%% Calculate a rotation matrix from GCS to PCS
R_gcs2pcs = zeros(3,3, 'double');

%--------------------------------------------------------------------------
% PE direction
%--------------------------------------------------------------------------
dGp = 1; dGr = 0; dGs = 0;
[dSag, dCor, dTra] = transformGCSToPCS(dGp, dGr, dGs, dNormalSag, dNormalCor, dNormalTra, dRotAngle);
R_gcs2pcs(1,1) = dSag;
R_gcs2pcs(2,1) = dCor;
R_gcs2pcs(3,1) = dTra;

%--------------------------------------------------------------------------
% RO direction
%--------------------------------------------------------------------------
dGp = 0; dGr = 1; dGs = 0;
[dSag, dCor, dTra] = transformGCSToPCS(dGp, dGr, dGs, dNormalSag, dNormalCor, dNormalTra, dRotAngle);
R_gcs2pcs(1,2) = dSag;
R_gcs2pcs(2,2) = dCor;
R_gcs2pcs(3,2) = dTra;

%--------------------------------------------------------------------------
% SS direction
%--------------------------------------------------------------------------
dGp = 0; dGr = 0; dGs = 1;
[dSag, dCor, dTra] = transformGCSToPCS(dGp, dGr, dGs, dNormalSag, dNormalCor, dNormalTra, dRotAngle);
R_gcs2pcs(1,3) = dSag;
R_gcs2pcs(2,3) = dCor;
R_gcs2pcs(3,3) = dTra;

%--------------------------------------------------------------------------
% Change the sign of gradients
%--------------------------------------------------------------------------
R_gcs2pcs(:,1) = phase_sign * R_gcs2pcs(:,1);
R_gcs2pcs(:,2) = read_sign  * R_gcs2pcs(:,2);

end
