%--------------------------------------------------------------------------
% Compile iterative_image_reconstruction_precompute.c
%--------------------------------------------------------------------------
debug_mode = 0;

computer_type = computer;

if strcmp(computer_type, 'PCWIN64')
    blas_lib = 'C:\Program Files\MATLAB\R2020b\extern\lib\win64\microsoft\libmwblas.lib';
    if debug_mode
        mex('-v', '-R2018a', 'COMPFLAGS="$COMPFLAGS -DDEBUG"', 'iterative_image_reconstruction_precompute.c', blas_lib);
    else
        mex('-v', '-R2018a', 'COMPFLAGS="$COMPFLAGS /arch:AVX2"', 'iterative_image_reconstruction_precompute.c', blas_lib);
    end
elseif strcmp(computer_type, 'MACI64')
    mex -v -R2018a iterative_image_reconstruction_precompute.c -lmwblas
elseif strcmp(computer_type, 'GLNXA64')
    mex -v -R2018a iterative_image_reconstruction_precompute.c -lmwblas
end
