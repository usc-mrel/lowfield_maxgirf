function outp = encoding_sense(inp, csm, w, nufft_st, transpose_indicator)
% Written by Nam Gyun Lee
% Email: nmgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/16/2022, Last modified: 01/16/2022

%% Declare a persistent variable
persistent cg_iter;
if isempty(cg_iter)
    cg_iter = 0;
end

%% Define dimensions
[N1,N2,Nc] = size(csm);
N = N1 * N2; % number of spatial locations in a 2D grid
Nd = nufft_st.Nd;

%% Determine the operator type
if strcmp(transpose_indicator, 'transp')
    operator_type = 'adjoint';
elseif strcmp(transpose_indicator, 'notransp')
    operator_type = 'forward';
    cg_iter = cg_iter + 1;
end

%% Calculate the forward or adjoint operator
outp = complex(zeros(N, 1, 'double'));

start_time = tic;
for c = 1:Nc
    tic; fprintf('(CG=%2d): Calculating the %s operator (c=%2d/%2d)... ', cg_iter, operator_type, c, Nc);
    %----------------------------------------------------------------------
    % Calculate Sc * m (N x 1)
    %----------------------------------------------------------------------
    Sm = reshape(csm(:,:,c), [N 1]) .* inp;

    %----------------------------------------------------------------------
    % Calculate F * (Sc * m) (Nk*Ni x 1)
    %----------------------------------------------------------------------
    FSm = nufft(reshape(Sm, [N1 N2]), nufft_st) / sqrt(prod(Nd)); % Nk*Ni x 1

    %----------------------------------------------------------------------
    % Calculate F^H * (F * Sc * m) (N x 1)
    %----------------------------------------------------------------------
    % Preconditioning with density compensation
    AHd = reshape(nufft_adj(FSm .* w(:), nufft_st) / sqrt(prod(Nd)), [N 1]); % N1 x N2 => N x 1

    %----------------------------------------------------------------------
    % Calculate Sc^H * (F^H * F * Sc * m)
    %----------------------------------------------------------------------
    AHd = reshape(conj(csm(:,:,c)), [N 1]) .* AHd; % N x 1

    %----------------------------------------------------------------------
    % Calculate output (N x 1)
    %----------------------------------------------------------------------
    outp = outp + AHd;
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
end

end