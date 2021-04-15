function [rotMatrixGCSToPCS,PE_sign,RO_sign,main_orientation] = calcMatrixGCSToPCS(dNormalSag, dNormalCor, dNormalTra, dRotAngle)
% Calculate a rotation matrix from GCS to PCS
% [Sag]                       [PE]
% [Cor] = rotMatrixGCSToPCS * [RO]
% [Tra]                       [SS]

%% Define constants
SAGITTAL   = 0; % Patient axis perpendicular to the sagittal plane
CORONAL    = 1; % Patient axis perpendicular to the coroal plane
TRANSVERSE = 2; % Patient axis perpendicular to the transvers plane

%% Calculate the signs of gradient waveforms due to ("swap Fre/Pha")
main_orientation = fGSLClassOri(dNormalSag, dNormalCor, dNormalTra);

switch main_orientation
    case SAGITTAL % =0
        if dRotAngle >= 0
            if ((dRotAngle < pi/4) && ~fGSLAlmEqual(dRotAngle, pi/4)) % 0 <= angle < 45
                RO_sign = -1;
                PE_sign =  1;
            elseif (((dRotAngle > pi/4) || fGSLAlmEqual(dRotAngle, pi/4)) && (dRotAngle < 3*pi/4)) % 45 <= angle < 135
                RO_sign =  1;
                PE_sign =  1;
            elseif ((dRotAngle > 3*pi/4) || fGSLAlmEqual(dRotAngle, 3*pi/4)) % 135 <= angle <= 180
                RO_sign =  1;
                PE_sign = -1;
            end
        elseif dRotAngle < 0
            if ((dRotAngle > -pi/4) && ~fGSLAlmEqual(dRotAngle, -pi/4)) % -45 < angle <= 0
                RO_sign = -1;
                PE_sign =  1;
            elseif ((dRotAngle > -3*pi/4) && ((dRotAngle < -pi/4) || fGSLAlmEqual(dRotAngle, -pi/4))) % -135 <= angle <= -45
                RO_sign = -1;
                PE_sign = -1;
            elseif (dRotAngle < -3*pi/4) % angle < -135
                RO_sign =  1;
                PE_sign = -1;
            end
        end
    case CORONAL % =1
            if ((dRotAngle >= -pi/4) && (dRotAngle <= pi/4)) % -45 <= angle <= 45
                RO_sign =  1;
                PE_sign =  1;
            else % -180 <= angle < -45 or 45 < angle <= 180
                RO_sign = -1;
                PE_sign = -1;
            end
    case TRANSVERSE % =2
        if dRotAngle >= 0
            if (dRotAngle <= pi/4) % 0 <= angle <= 45
                RO_sign = -1;
                PE_sign =  1;
            elseif ((dRotAngle > pi/4) && (dRotAngle < 3*pi/4) && ~fGSLAlmEqual(dRotAngle, 3*pi/4)) % 45 < angle < 135
                RO_sign =  1;
                PE_sign =  1;
            elseif ((dRotAngle > 3*pi/4) || fGSLAlmEqual(dRotAngle, 3*pi/4)) % 135 <= angle <= 180
                RO_sign =  1;
                PE_sign = -1;
            end
        elseif dRotAngle < 0
            if ((dRotAngle > -pi/4) || fGSLAlmEqual(dRotAngle, -pi/4)) % -45 <= angle <= 0
                RO_sign = -1;
                PE_sign =  1;
            elseif ((dRotAngle > -3*pi/4) && ~fGSLAlmEqual(dRotAngle, -3*pi/4) && (dRotAngle < -pi/4)) % -135 < angle < -45
                RO_sign = -1;
                PE_sign = -1;
            else % -180 <= angle < -135
                RO_sign =  1;
                PE_sign = -1;
            end
        end
    otherwise
end

%% Calculate a rotation matrix from GCS to PCS
rotMatrixGCSToPCS = zeros(3,3, 'double');

%--------------------------------------------------------------------------
% PE direction
%--------------------------------------------------------------------------
dGp = 1; dGr = 0; dGs = 0;
[dSag, dCor, dTra] = transformGCSToPCS(dGp, dGr, dGs, dNormalSag, dNormalCor, dNormalTra, dRotAngle);
rotMatrixGCSToPCS(1,1) = dSag;
rotMatrixGCSToPCS(2,1) = dCor;
rotMatrixGCSToPCS(3,1) = dTra;

%--------------------------------------------------------------------------
% RO direction
%--------------------------------------------------------------------------
dGp = 0; dGr = 1; dGs = 0;
[dSag, dCor, dTra] = transformGCSToPCS(dGp, dGr, dGs, dNormalSag, dNormalCor, dNormalTra, dRotAngle);
rotMatrixGCSToPCS(1,2) = dSag;
rotMatrixGCSToPCS(2,2) = dCor;
rotMatrixGCSToPCS(3,2) = dTra;

%--------------------------------------------------------------------------
% SS direction
%--------------------------------------------------------------------------
dGp = 0; dGr = 0; dGs = 1;
[dSag, dCor, dTra] = transformGCSToPCS(dGp, dGr, dGs, dNormalSag, dNormalCor, dNormalTra, dRotAngle);
rotMatrixGCSToPCS(1,3) = dSag;
rotMatrixGCSToPCS(2,3) = dCor;
rotMatrixGCSToPCS(3,3) = dTra;

%--------------------------------------------------------------------------
% Change the sign of gradients
%--------------------------------------------------------------------------
rotMatrixGCSToPCS(:,1) = PE_sign * rotMatrixGCSToPCS(:,1);
rotMatrixGCSToPCS(:,2) = RO_sign * rotMatrixGCSToPCS(:,2);

end
