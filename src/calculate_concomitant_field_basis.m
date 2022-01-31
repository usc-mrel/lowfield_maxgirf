function p = calculate_concomitant_field_basis(x, y, z, Nl)
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/17/2022, Last modified: 01/17/2022

N = length(x);
p = zeros(N, Nl, 'double');
if Nl >= 1 , p(:,1)  = x;           end
if Nl >= 2 , p(:,2)  = y;           end
if Nl >= 3 , p(:,3)  = z;           end
if Nl >= 4 , p(:,4)  = x.^2;        end
if Nl >= 5 , p(:,5)  = y.^2;        end
if Nl >= 6 , p(:,6)  = z.^2;        end
if Nl >= 7 , p(:,7)  = x .* y;      end
if Nl >= 8 , p(:,8)  = y .* z;      end
if Nl >= 9 , p(:,9)  = x .* z;      end
if Nl >= 10, p(:,10)  = x.^3;       end
if Nl >= 11, p(:,11)  = y.^3;       end
if Nl >= 12, p(:,12)  = z.^3;       end
if Nl >= 13, p(:,13) = x.^2 .* y;   end
if Nl >= 14, p(:,14) = x.^2 .* z;   end
if Nl >= 15, p(:,15) = x.* y.^2;    end
if Nl >= 16, p(:,16) = y.^2 .* z;   end
if Nl >= 17, p(:,17) = x.* z.^2;    end
if Nl >= 18, p(:,18) = y.* z.^2;    end
if Nl >= 19, p(:,19) = x .* y .* z; end

end
