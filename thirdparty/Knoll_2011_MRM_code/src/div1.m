function u = div1(v,h)
% divergence operator: V -> U (C^2MN -> C^MN)

u = dxn(v(:,:,1),h) + dyn(v(:,:,2),h);

end