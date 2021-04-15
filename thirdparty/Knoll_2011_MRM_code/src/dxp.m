function dxpu = dxp(u,h)
% Calculate the forward difference along the first dimension
% u in U = C^(M x N)

dxpu = (u([2:end end],:,:) - u) / h;

end