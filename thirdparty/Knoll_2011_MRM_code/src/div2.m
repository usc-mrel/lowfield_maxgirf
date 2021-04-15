function v = div2(w,h)
% the adjoint of symmetrized second-order derivative operator
% div2: W -> V (C^3MN -> C^2MN)

v = cat(3, dxp(w(:,:,1),h) + dyp(w(:,:,3),h), ...
           dxp(w(:,:,3),h) + dyp(w(:,:,2),h));
end