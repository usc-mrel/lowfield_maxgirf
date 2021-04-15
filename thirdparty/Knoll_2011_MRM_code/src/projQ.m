function q = projQ(q,alpha0)
% Euclidean projector onto the convex set Q, where Q = {q in C^(3MN): ||q||_inf <= alpha0}

q1 = q(:,:,1);
q2 = q(:,:,2);
q3 = q(:,:,3);
abs_q = sqrt(abs(q1).^2 + abs(q2).^2 + 2 * abs(q3).^2);

q = q ./ max(1, abs_q / alpha0);

end