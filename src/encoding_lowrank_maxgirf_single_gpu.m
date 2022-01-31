function outp = encoding_lowrank_maxgirf_single_gpu(inp, csm_device, U_device, V_device, w_device, st_device, transpose_indicator)
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/18/2022, Last modified: 01/18/2022

%% Declare a persistent variable
persistent cg_iter;
if isempty(cg_iter)
    cg_iter = 0;
end

%% Define dimensions
[N1,N2,Nc] = size(csm_device);
N = N1 * N2; % number of spatial locations in a 2D grid
[Nk,L,Ni] = size(U_device);

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
for i = 1:Ni
    Nd = st_device{i}.Nd;
    tic; fprintf('(CG=%2d): Calculating the %s operator (i=%2d/%2d)... ', cg_iter, operator_type, i, Ni);
    for c = 1:Nc
        %------------------------------------------------------------------
        % Calculate Sc * m (N x 1)
        %------------------------------------------------------------------
        Sm = reshape(csm_device(:,:,c), [N 1]) .* inp;

        %------------------------------------------------------------------
        % Calculate sum_{ell=1}^L ...
        % diag(U(:,ell,i)) * Fi * diag(conj(V(:,ell,i))) * Sc * m
        %------------------------------------------------------------------
        scale_factor = gpuArray(1 / sqrt(prod(Nd)));
        ESm = complex(zeros(Nk, 1, 'double', 'gpuArray'));
        for ell = 1:L
            FDvHSm = nufft_gpu(reshape(conj(V_device(:,ell,i)) .* Sm, [N1 N2]), st_device{i}) * scale_factor; % Nk x 1
            ESm = ESm + U_device(:,ell,i) .* FDvHSm;
        end

        %------------------------------------------------------------------
        % Calculate sum_{ell=1}^L ...
        % diag(V(:,ell,i)) * Fi^H * diag(conj(U(:,ell,i))) * (Ei * Sc * m)
        %------------------------------------------------------------------
        AHd = complex(zeros(N, 1, 'double', 'gpuArray'));
        for ell = 1:L
            % Preconditioning with density compensation
            FHDuHd = nufft_adj_gpu((conj(U_device(:,ell,i)) .* ESm) .* w_device(:,i), st_device{i}) * scale_factor;
            AHd = AHd + V_device(:,ell,i) .* reshape(FHDuHd, [N 1]);
        end

        %------------------------------------------------------------------
        % Calculate Sc^H * (Ei^H * Ei * Sc * m)
        %------------------------------------------------------------------
        AHd = reshape(conj(csm_device(:,:,c)), [N 1]) .* AHd;

        %------------------------------------------------------------------
        % Calculate output (N x 1)
        %------------------------------------------------------------------
        outp = outp + AHd;
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
end

end