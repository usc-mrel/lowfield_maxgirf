function plCase = fGSLClassOri(dSagComp, dCorComp, dTraComp)
%
% Description : This function determines whether a supplied normal vector
%                 describes a sagittal, coronal or transverse slice.
%       Result:
%         CASE = 0: Sagittal
%         CASE = 1: Coronal
%         CASE = 2: Transverse
%
%         This was created with th following table:
%           (Note: ~means 'about equal')
%
%         |sag|~|cor| and |sag|~|tra|  -->  tra
%
%         |sag|~|cor| and |sag|<|tra|  -->  tra
%         |sag|~|cor| and |sag|>|tra|  -->  cor
%
%         |sag|~|tra| and |sag|<|cor|  -->  cor
%         |sag|~|tra| and |sag|>|cor|  -->  tra
%
%         |cor|~|tra| and |cor|<|sag|  -->  sag
%         |cor|~|tra| and |cor|>|sag|  -->  tra
%
%         |sag|>|cor| and |sag|>|tra|  -->  sag
%         |sag|>|cor| and |sag|<|tra|  -->  tra
%         |sag|<|cor| and |cor|>|tra|  -->  cor
%         |sag|<|cor| and |cor|<|tra|  -->  tra
%
%         |sag|>|tra| and |sag|<|cor|  -->  cor
%         |sag|>|tra| and |sag|>|cor|  -->  sag
%         |sag|<|tra| and |tra|<|cor|  -->  cor
%         |sag|<|tra| and |tra|>|cor|  -->  tra
%
%         |cor|>|tra| and |cor|<|sag|  -->  sag
%         |cor|>|tra| and |cor|>|sag|  -->  cor
%         |cor|<|tra| and |tra|<|sag|  -->  sag
%         |cor|<|tra| and |tra|>|sag|  -->  tra
%
% Return      : NLS error status
%
%--------------------------------------------------------------------------
% dSagComp    IMP: Sagittal component of normal vector
% dCorComp    IMP: Coronal component of normal vector
% dTraComp    IMP: Transverse component of normal vector
% plCase      EXP: Case (0=Sagittal, 1=Coronal or 2=Transverse)
%--------------------------------------------------------------------------

SAGITTAL   = 0;
CORONAL    = 1;
TRANSVERSE = 2;

%--------------------------------------------------------------------------
% For performance reasons, calculate some tmp values
%--------------------------------------------------------------------------
dAbsSagComp     = abs(dSagComp);
dAbsCorComp     = abs(dCorComp);
dAbsTraComp     = abs(dTraComp);
bAlmEqualSagCor = fGSLAlmEqual(dAbsSagComp, dAbsCorComp);
bAlmEqualSagTra = fGSLAlmEqual(dAbsSagComp, dAbsTraComp);
bAlmEqualCorTra = fGSLAlmEqual(dAbsCorComp, dAbsTraComp);

%--------------------------------------------------------------------------
% Check all values to determine the slice orientation (sag, cor, tra)
%--------------------------------------------------------------------------
if ((bAlmEqualSagCor              &&  bAlmEqualSagTra)             || ...
    (bAlmEqualSagCor              &&  (dAbsSagComp < dAbsTraComp)) || ...
    (bAlmEqualSagTra              &&  (dAbsSagComp > dAbsCorComp)) || ...
    (bAlmEqualCorTra              &&  (dAbsCorComp > dAbsSagComp)) || ...
    ((dAbsSagComp > dAbsCorComp)  &&  (dAbsSagComp < dAbsTraComp)) || ...
    ((dAbsSagComp < dAbsCorComp)  &&  (dAbsCorComp < dAbsTraComp)) || ...
    ((dAbsSagComp < dAbsTraComp)  &&  (dAbsTraComp > dAbsCorComp)) || ...
    ((dAbsCorComp < dAbsTraComp)  &&  (dAbsTraComp > dAbsSagComp)))

    %----------------------------------------------------------------------
    % Mainly transverse...
    %----------------------------------------------------------------------
    plCase = TRANSVERSE;

elseif ((bAlmEqualSagCor              &&  (dAbsSagComp > dAbsTraComp)) || ...
        (bAlmEqualSagTra              &&  (dAbsSagComp < dAbsCorComp)) || ...
        ((dAbsSagComp < dAbsCorComp)  &&  (dAbsCorComp > dAbsTraComp)) || ...
        ((dAbsSagComp > dAbsTraComp)  &&  (dAbsSagComp < dAbsCorComp)) || ...
        ((dAbsSagComp < dAbsTraComp)  &&  (dAbsTraComp < dAbsCorComp)))

    %----------------------------------------------------------------------
    % Mainly coronal...
    %----------------------------------------------------------------------
    plCase = CORONAL;

elseif ((bAlmEqualCorTra              &&  (dAbsCorComp < dAbsSagComp)) || ...
        ((dAbsSagComp > dAbsCorComp)  &&  (dAbsSagComp > dAbsTraComp)) || ...
        ((dAbsCorComp > dAbsTraComp)  &&  (dAbsCorComp < dAbsSagComp)) || ...
        ((dAbsCorComp < dAbsTraComp)  &&  (dAbsTraComp < dAbsSagComp)))

    %----------------------------------------------------------------------
    % Mainly coronal...
    %----------------------------------------------------------------------
    plCase = SAGITTAL;
end

end