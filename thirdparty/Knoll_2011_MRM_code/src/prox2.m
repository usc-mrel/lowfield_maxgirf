function prox2_r = prox2(r,sigma,lambda)
% the proximal map

prox2_r = (1 / (1 + sigma * lambda)) * r;

end