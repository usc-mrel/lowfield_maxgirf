function prox1_u = prox1(u,f,tau,lambda)
% the proximal map

prox1_u = (lambda * u + tau * f) ./ (lambda + tau);

end