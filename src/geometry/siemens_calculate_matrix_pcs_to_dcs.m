function R_pcs2dcs = siemens_calculate_matrix_pcs_to_dcs(patient_position)
% [X]               [Sag]
% [Y] = R_pcs2dcs * [Cor]
% [z]               [Tra]
%
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/16/2022, Last modified: 01/16/2022

switch patient_position
    case 'HFP' % head first / prone
        R_pcs2dcs = [-1    0    0 ;
                      0    1    0 ;
                      0    0   -1];
    case 'HFS' % head first / supine
        R_pcs2dcs = [ 1    0    0 ;
                      0   -1    0 ;
                      0    0   -1];
    case 'HFDR' % head first / decubitus right
        R_pcs2dcs = [ 0    1    0 ;
                      1    0    0 ;
                      0    0   -1];
    case 'HFDL' % head first / decubitus left
        R_pcs2dcs = [ 0   -1    0 ;
                     -1    0    0 ;
                      0    0   -1];
    case 'FFP' % feet first / prone
        R_pcs2dcs = [ 1    0    0 ;
                      0    1    0 ;
                      0    0    1];
    case 'FFS' % feet first / supine
        R_pcs2dcs = [-1    0    0 ;
                      0   -1    0 ;
                      0    0    1];
    case 'FFDR' % feet first / decubitus right
        R_pcs2dcs = [ 0   -1    0 ;
                      1    0    0 ;
                      0    0    1];
    case 'FFDL' % feet first / decubitus left
        R_pcs2dcs = [ 0    1    0 ;
                     -1    0    0 ;
                      0    0    1];
    otherwise
end

end