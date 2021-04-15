function p = projP(p,alpha1)
% Euclidean projector onto the convex set P, where P = {p in C^(2MN): ||p||_inf <= alpha1}

p1 = p(:,:,1);
p2 = p(:,:,2);
abs_p = sqrt(abs(p1).^2 + abs(p2).^2);

p = p ./ max(1, abs_p / alpha1);

end