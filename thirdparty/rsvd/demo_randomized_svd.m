% demo_randomized_svd.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/01/2021, Last modified: 01/01/2021

%% Clean slate
close all; clear; clc;

%% Set parameters
m = 1000;  % number of time frames
n = 1000;  % number of tissue-type parameters
q = 1;     % order of power iteration
k = 40;    % rank of approximation
p = 5;     % oversampling parameter

%% Calculate a dictionary
D = phantom(m) + 1j * phantom(m);

%% Initialization
B = complex(zeros(k+p, n  , 'double'));
W = complex(zeros(m  , k+p, 'double'));
Z = complex(zeros(m  , m  , 'double')); % D * D' => (m x n) * (m x n)' => m x m

%% Calculate W and Z
for idx = 1:n
    %----------------------------------------------------------------------
    % Calculate a dictionary atom
    %----------------------------------------------------------------------
    d = D(:,idx); % m x 1

    %----------------------------------------------------------------------
    % Draw a 1 x k Gaussian random vector w
    %----------------------------------------------------------------------
    w = randn(1, k+p, 'double'); % 1 x k+p

    %----------------------------------------------------------------------
    % Calculate W
    %----------------------------------------------------------------------
    W = W + d * w; % m x k+p

    %----------------------------------------------------------------------
    % Calculate Z
    %----------------------------------------------------------------------
    if q > 0
        Z = Z + d * d';
    end
end

%% Calculate Y (m x k+p)
Y = Z^q * W; % (m x m) * (m x k+p) => m x k+p

%% Calculate the QR decomposition of Y
[Q,R] = qr(Y,0); % Q:m x k+p, R:k+p x k+p

%% Calculate B (k+p x n)
for idx = 1:n
    d = D(:,idx);
    B(:,idx) = Q' * d;
end

%% Calculate the SVD of B
[U_hat,S,V] = svd(B,0); % U_hat:k x k, S:k x k, V:k x n

%% Calculate U
U = Q * U_hat; % (m x k+p) * (k+p x k+p) => m x k+p

%% Calculate the deterministic SVD of D
[U2,S2,V2] = svd(D,0);

%% Calculate the approximations
D_rsvd = U(:,1:k) * S(1:k,1:k) * V(:,1:k)';
D_svd  = U2(:,1:k) * S2(1:k,1:k) * V2(:,1:k)';

%% Display results
figure('Color', 'w');
imagesc(real(cat(2, D, D_rsvd, D_svd, (D_rsvd - D_svd)*1e1))); 
caxis([0 1]); axis image;
impixelinfo;

figure('Color', 'w');
imagesc(imag(cat(2, D, D_rsvd, D_svd, (D_rsvd - D_svd)*1e1)));
caxis([0 1]); axis image;
impixelinfo;
