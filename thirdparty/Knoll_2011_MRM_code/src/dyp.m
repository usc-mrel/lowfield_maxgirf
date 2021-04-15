function dypu = dyp(u,h)
% Calculate the forward difference along the second dimension
% u in U = C^(M x N)

dypu = (u(:,[2:end end],:) - u) / h;

end