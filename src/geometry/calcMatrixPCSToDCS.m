function rotMatrixPCSToDCS = calcMatrixPCSToDCS(patient_position)
% [X]                       [Sag]
% [Y] = rotMatrixPCSToDCS * [Cor]
% [z]                       [Tra]


switch patient_position
    case 'HFP' % head first / prone
        rotMatrixPCSToDCS = [-1    0    0 ;
                              0    1    0 ;
                              0    0   -1];
    case 'HFS' % head first / supine
        rotMatrixPCSToDCS = [ 1    0    0 ;
                              0   -1    0 ;
                              0    0   -1];
    case 'HFDR' % head first / decubitus right
        rotMatrixPCSToDCS = [ 0    1    0 ;
                              1    0    0 ;
                              0    0   -1];
    case 'HFDL' % head first / decubitus left
        rotMatrixPCSToDCS = [ 0   -1    0 ;
                             -1    0    0 ;
                              0    0   -1];
    case 'FFP' % feet first / prone
        rotMatrixPCSToDCS = [ 1    0    0 ;
                              0    1    0 ;
                              0    0    1];
    case 'FFS' % feet first / supine
        rotMatrixPCSToDCS = [-1    0    0 ;
                              0   -1    0 ;
                              0    0    1];
    case 'FFDR' % feet first / decubitus right
        rotMatrixPCSToDCS = [ 0   -1    0 ;
                              1    0    0 ;
                              0    0    1];
    case 'FFDL' % feet first / decubitus left
        rotMatrixPCSToDCS = [ 0    1    0 ;
                             -1    0    0 ;
                              0    0    1];
    otherwise
end

end
