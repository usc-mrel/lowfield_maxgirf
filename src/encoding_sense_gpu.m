function outp = encoding_sense_gpu(inp, csm_device, w_device, nufft_st_device, transpose_indicator)
% Written by Nam Gyun Lee
% Email: nmgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/16/2022, Last modified: 01/16/2022

%% Declare a persistent variable
persistent cg_iter;
if isempty(cg_iter)
    cg_iter = 0;
end

%% Define dimensions
[N1,N2,Nc] = size(csm_device);
N = N1 * N2; % number of spatial locations in a 2D grid
Nd = nufft_st_device.Nd;

%% Determine the operator type
if strcmp(transpose_indicator, 'transp')
    operator_type = 'adjoint';
elseif strcmp(transpose_indicator, 'notransp')
    operator_type = 'forward';
    cg_iter = cg_iter + 1;
end

%% Calculate the forward or adjoint operator
outp = complex(zeros(N, 1, 'double', 'gpuArray'));

start_time = tic;
for c = 1:Nc
    tic; fprintf('(CG=%2d): Calculating the %s operator (c=%2d/%2d)... ', cg_iter, operator_type, c, Nc);
    %----------------------------------------------------------------------
    % Calculate Sc * m (N x 1)
    %----------------------------------------------------------------------
    Sm = reshape(csm_device(:,:,c), [N 1]) .* inp;

    %----------------------------------------------------------------------
    % Calculate F * (Sc * m) (Nk*Ni x 1)
    %----------------------------------------------------------------------
    scale_factor = gpuArray(1 / sqrt(prod(Nd)));
    FSm = nufft_gpu(reshape(Sm, [N1 N2]), nufft_st_device) * scale_factor; % Nk*Ni x 1

    %----------------------------------------------------------------------
    % Calculate F^H * (F * Sc * m) (N x 1)
    %----------------------------------------------------------------------
    % Preconditioning with density compensation
    AHd = reshape(nufft_adj_gpu(FSm .* w_device(:), nufft_st_device) * scale_factor, [N 1]); % N1 x N2 => N x 1

    %----------------------------------------------------------------------
    % Calculate Sc^H * (F^H * F * Sc * m)
    %----------------------------------------------------------------------
    AHd = reshape(conj(csm_device(:,:,c)), [N 1]) .* AHd; % N x 1

    %----------------------------------------------------------------------
    % Calculate output (N x 1)
    %----------------------------------------------------------------------
    outp = outp + AHd;
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
end

end