function dxnu = dxn(u,h)
% Calculate the backward difference along the first dimension
% u in U = C^(M x N)

% 0 < i < M - 1
dxnu = (u - u([1 1:end-1],:,:)) / h;
% i = 0
dxnu(1,:,:) = u(1,:,:) / h;
% i = M - 1
dxnu(end,:,:) = -u(end-1,:,:) / h;

end