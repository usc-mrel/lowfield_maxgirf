function outp = encoding_lowrank_maxgirf(inp, csm, U, V, w, st, transpose_indicator)
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/18/2022, Last modified: 01/18/2022

%% Declare a persistent variable
persistent cg_iter;
if isempty(cg_iter)
    cg_iter = 0;
end

%% Define dimensions
[N1,N2,Nc] = size(csm);
N = N1 * N2; % number of spatial locations in a 2D grid
[Nk,L,Ni] = size(U);

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
for i = 1:Ni
    Nd = st{i}.Nd;
    tic; fprintf('(CG=%2d): Calculating the %s operator (i=%2d/%2d)... ', cg_iter, operator_type, i, Ni);
    for c = 1:Nc
        %------------------------------------------------------------------
        % Calculate Sc * m (N x 1)
        %------------------------------------------------------------------
        Sm = reshape(csm(:,:,c), [N 1]) .* inp;

        %------------------------------------------------------------------
        % Calculate sum_{ell=1}^L ...
        % diag(U(:,ell,i)) * Fi * diag(conj(V(:,ell,i))) * Sc * m
        %------------------------------------------------------------------
        ESm = complex(zeros(Nk, 1, 'double'));
        for ell = 1:L
            FDvHSm = nufft(reshape(conj(V(:,ell,i)) .* Sm, [N1 N2]), st{i}) / sqrt(prod(Nd)); % Nk x 1
            ESm = ESm + U(:,ell,i) .* FDvHSm;
        end

        %------------------------------------------------------------------
        % Calculate sum_{ell=1}^L ...
        % diag(V(:,ell,i)) * Fi^H * diag(conj(U(:,ell,i))) * (Ei * Sc * m)
        %------------------------------------------------------------------
        AHd = complex(zeros(N, 1, 'double'));
        for ell = 1:L
            % Preconditioning with density compensation
            FHDuHd = nufft_adj((conj(U(:,ell,i)) .* ESm) .* w(:,i), st{i}) / sqrt(prod(Nd));
            AHd = AHd + V(:,ell,i) .* reshape(FHDuHd, [N 1]);
        end

        %------------------------------------------------------------------
        % Calculate Sc^H * (Ei^H * Ei * Sc * m)
        %------------------------------------------------------------------
        AHd = reshape(conj(csm(:,:,c)), [N 1]) .* AHd;

        %------------------------------------------------------------------
        % Calculate output (N x 1)
        %------------------------------------------------------------------
        outp = outp + AHd;
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
end

end