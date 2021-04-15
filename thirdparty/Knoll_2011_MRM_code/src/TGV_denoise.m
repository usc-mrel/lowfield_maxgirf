function u = TGV_denoise(f, lambda, maxiter, verbose)
% TGV_denoise -- Primal-dual method for TGV denoising
% (c) Nam Gyun Lee 2019
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Created: 05/17/2019, Last modified: 05/17/2019


[M,N] = size(f);

alpha0 = 2;
alpha1 = 1;
h      = 1; % mesh-size
sigma  = 1 / sqrt(12);
tau    = 1 / sqrt(12);

%--------------------------------------------------------------------------
% Initialization
%--------------------------------------------------------------------------
u     = f; % M x N
u_bar = f; % M x N

v     = complex(zeros(M,N,2, 'double')); % M x N x 2
v_bar = complex(zeros(M,N,2, 'double')); % M x N x 2

p = complex(zeros(M,N,2, 'double')); % M x N x 2
q = complex(zeros(M,N,3, 'double')); % M x N x 3

%--------------------------------------------------------------------------
% main loop
%--------------------------------------------------------------------------
for k = 1:maxiter

    start_time = tic;

    % Algorithm 1 (primal-dual method for TGV denoising)
    p = projP(p + sigma * (grad(u_bar,h) - v_bar), alpha1);
    q = projQ(q + sigma * symmetrized_derivative(v_bar,h), alpha0);
    u_old = u;
    u = prox1(u + tau * div1(p,h), f, tau, lambda);
    u_bar = 2 * u - u_old;
    v_old = v;
    v = v + tau * (p + div2(q,h));
    v_bar = 2 * v - v_old;

    elapsed_time = toc(start_time); 

    %----------------------------------------------------------------------
    % Display progress
    %----------------------------------------------------------------------
    if (verbose)
        fprintf('iteration = %4d, elapsed time: %7.4f sec\n', k, elapsed_time);
    end
end

end