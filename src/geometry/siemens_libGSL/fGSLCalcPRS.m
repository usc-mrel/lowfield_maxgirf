function [dGp, dGr] = fGSLCalcPRS(dGs, dPhi)
%--------------------------------------------------------------------------
%
% Name        : fGSLCalcLin
%
% Description : calculates the two vectors of phase encoding and and
%               readout direction.
%               The calculation depends on the slice orientation (the slice
%               normal vector) and on the angle of rotation around the s axis.
%		            Every vector (Gp, Gr and Gs) has the three components sag, cor
%               and tra, the description is patient oriented. All three
%               vectors have a length of 1.0. The biggest component of Gs must
%               have a positive sign.
%
%		Formulas for the rotation around the vector s:
%		    (a = cos (dPhi), b = sin (dPhi))
%
%		    new             old               rotation     base
%		    vector  =   coordinate system   *  matrix    * vector
%
%		    (P_sag)   (P_sag  R_sag  S_sag)   ( a  b  0)   (1)
%		    (P_cor) = (P_cor  R_cor  S_cor) * (-b  a  0) * (0)
%		    (P_tra)   (P_tra  R_tra  S_tra)   ( 0  0  1)   (0)
%
%		    (R_sag)   (P_sag  R_sag  S_sag)   ( a  b  0)   (0)
%		    (R_cor) = (P_cor  R_cor  S_cor) * (-b  a  0) * (1)
%		    (R_tra)   (P_tra  R_tra  S_tra)   ( 0  0  1)   (0)
%
%		    (S_sag)   (P_sag  R_sag  S_sag)   ( a  b  0)   (0)
%		    (S_cor) = (P_cor  R_cor  S_cor) * (-b  a  0) * (0)
%		    (S_tra)   (P_tra  R_tra  S_tra)   ( 0  0  1)   (1)
%
%		    This multiplied:
%
%		    (P_sag)   (a * P_sag - b * R_sag)
%		    (P_cor) = (a * P_cor - b * R_cor)
%		    (P_tra)   (a * P_tra - b * R_tra)
%
%		    (R_sag)   (b * P_sag + a * R_sag)
%		    (R_cor) = (b * P_cor + a * R_cor)
%		    (R_tra)   (b * P_tra + a * R_tra)
%
%		    (S_sag)   (S_sag)
%		    (S_cor) = (S_cor)	well!
%		    (S_tra)   (P_tra)
%
% Return      : NLS error status
%
%--------------------------------------------------------------------------
% dGp[3]    EXP: The GP vector
% dGr[3]    EXP: The GR vector
% dGs[3]    IMP: The GS vector (= slice normal vector)
% dPhi      IMP: The rotational angle around GS
%--------------------------------------------------------------------------
% Source code from MrProtMath.h (MrProt_fGSLCalcPRS)


SAGITTAL   = 0;
CORONAL    = 1;
TRANSVERSE = 2;

dGp = zeros(3, 1, 'double');
dGr = zeros(3, 1, 'double');

lOrientation = fGSLClassOri(dGs(SAGITTAL+1), dGs(CORONAL+1), dGs(TRANSVERSE+1));

switch (lOrientation)
    case TRANSVERSE
        dGp(1) = 0.;
        dGp(2) =  dGs(3) * sqrt (1. / (dGs(2) * dGs(2) + dGs(3) * dGs(3)));
        dGp(3) = -dGs(2) * sqrt (1. / (dGs(2) * dGs(2) + dGs(3) * dGs(3)));

    case CORONAL
        dGp(1) =  dGs(2) * sqrt (1. / (dGs(1) * dGs(1) + dGs(2) * dGs(2)));
        dGp(2) = -dGs(1) * sqrt (1. / (dGs(1) * dGs(1) + dGs(2) * dGs(2)));
        dGp(3) = 0.;

    case SAGITTAL
        dGp(1) = -dGs(2) * sqrt (1. / (dGs(1) * dGs(1) + dGs(2) * dGs(2)));
        dGp(2) =  dGs(1) * sqrt (1. / (dGs(1) * dGs(1) + dGs(2) * dGs(2)));
        dGp(3) = 0.;
    otherwise
end

% PRS => PE x RO x SS
%--------------------------------------------------------------------------
%  Calculate GR = GS x GP
%--------------------------------------------------------------------------
dGr(1) = dGs(2) * dGp(3) - dGs(3) * dGp(2);
dGr(2) = dGs(3) * dGp(1) - dGs(1) * dGp(3);
dGr(3) = dGs(1) * dGp(2) - dGs(2) * dGp(1); 

if (dPhi ~= 0.)
    %----------------------------------------------------------------------
    % Rotate around the S axis
    %----------------------------------------------------------------------
    dGp(1) = cos(dPhi) * dGp(1) - sin(dPhi) * dGr(1);
    dGp(2) = cos(dPhi) * dGp(2) - sin(dPhi) * dGr(2);
    dGp(3) = cos(dPhi) * dGp(3) - sin(dPhi) * dGr(3);

    %----------------------------------------------------------------------
    % Calculate new GR = GS x GP
    %----------------------------------------------------------------------
    dGr(1) = dGs(2) * dGp(3) - dGs(3) * dGp(2);
    dGr(2) = dGs(3) * dGp(1) - dGs(1) * dGp(3);
    dGr(3) = dGs(1) * dGp(2) - dGs(2) * dGp(1);
end 

end