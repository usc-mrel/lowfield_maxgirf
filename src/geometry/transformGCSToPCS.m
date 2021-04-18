function [dSag, dCor, dTra] = transformGCSToPCS(dGp, dGr, dGs, dNormalSag, dNormalCor, dNormalTra, dRotAngle)
%--------------------------------------------------------------------------
%
% Name        : transformGCSToPCS
%
% Description : Transforms a vector given in gradient coordinate system
%               into a vector in patient coordinate system.
%               For this the slice orientation is needed.
%
% Return      : NLS status code
%
%--------------------------------------------------------------------------
% dGp         Import: Phase encoding component
% dGr         Import: Readout component
% dGs         Import: Slice selection component
% dSag        Export: Sagittal component
% dCor        Export: Coronal component
% dTra        Export: Transverse component
% dNormalSag  Import: Sagittal component of slice normal vector
% dNormalCor  Import: Coronal component of slice normal vector
% dNormalTra  Import: Transverse component of slice normal vector
% dRotAngle   Import: Slice rotation angle ("swap Fre/Pha")
%--------------------------------------------------------------------------

dVectorGP = zeros(3,1, 'double'); % The GP vector (in patient coordinates)
dVectorGR = zeros(3,1, 'double'); % The GR vector (in patient coordinates)
dVectorGS = zeros(3,1, 'double'); % The GS vector (in patient coordinates)

dVectorGS(1) = dNormalSag;
dVectorGS(2) = dNormalCor;
dVectorGS(3) = dNormalTra;

[dVectorGP, dVectorGR] = fGSLCalcPRS(dVectorGS, dRotAngle);

dSag = dVectorGP(1) * dGp + dVectorGR(1) * dGr + dVectorGS(1) * dGs;
dCor = dVectorGP(2) * dGp + dVectorGR(2) * dGr + dVectorGS(2) * dGs;
dTra = dVectorGP(3) * dGp + dVectorGR(3) * dGr + dVectorGS(3) * dGs;

end
