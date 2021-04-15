function w = symmetrized_derivative(v,h)
% symmetrized second-order derivative operator: V -> W (C^2MN -> C^3MN)

[M,N,~] = size(v);

w = complex(zeros(M,N,3, 'double'));

w(:,:,1) = dxn(v(:,:,1),h); 
w(:,:,2) = dyn(v(:,:,2),h); 
w(:,:,3) = (dyn(v(:,:,1),h) + dxn(v(:,:,2),h)) / 2;

end