function v = grad(u,h)
% Gradient operator: U -> V (C^MN -> C^2MN)

[M,N] = size(u);
v = complex(zeros(M,N,2, 'double')); % M x N x 2

v(:,:,1) = dxp(u,h);
v(:,:,2) = dyp(u,h);

end