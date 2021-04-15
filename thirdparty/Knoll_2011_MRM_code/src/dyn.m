function dynu = dyn(u,h)
% Calculate the backward difference along the second dimension
% u in U = C^(M x N)

% 0 < j < N - 1
dynu = (u - u(:,[1 1:end-1],:)) / h;
% j = 0
dynu(:,1,:) = u(:,1,:) / h;
% j = N - 1
dynu(:,end,:) = -u(:,end-1,:) / h;

end