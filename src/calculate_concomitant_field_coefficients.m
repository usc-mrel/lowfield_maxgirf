function k = calculate_concomitant_field_coefficients(gx, gy, gz, Nl, B0, gamma, dt)
% gx, gy, gz in [G/cm]
% gamma in [rad/sec/T]
% dt in [sec]
%
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/17/2022, Last modified: 01/17/2022

%% Determine parameters
[Nk,Ni] = size(gx);
c = zeros(Nk, Nl, Ni, 'double');

%% Calculate the coefficients of gradients [T/m]
%--------------------------------------------------------------------------
% [G/cm] * [T/1e4G] * [1e2cm/m] => [T/m]
%--------------------------------------------------------------------------
if Nl >= 1, c(:,1,:) = gx * 1e-2; end % x
if Nl >= 2, c(:,2,:) = gy * 1e-2; end % y
if Nl >= 3, c(:,3,:) = gz * 1e-2; end % z

%% Calculate the coefficients of concomitant fields [T/m^2]
%--------------------------------------------------------------------------
% 1/B0 order concomitant gradient terms (x2, y2, z2, xy, yz, xz)
% 1 / [T] * ( ([G/cm] * [T/1e4G] * [1e2cm/m])^2 ) => [T/m^2]
%--------------------------------------------------------------------------
if Nl >= 4, c(:,4,:) = ( 1 / (8 * B0) * (gz * 1e-2).^2);                    end % x2
if Nl >= 5, c(:,5,:) = ( 1 / (8 * B0) * (gz * 1e-2).^2);                    end % y2
if Nl >= 6, c(:,6,:) = ( 1 / (2 * B0) * ((gx * 1e-2).^2 + (gy * 1e-2).^2)); end % z2
if Nl >= 7, c(:,7,:) = 0;                                                   end % xy
if Nl >= 8, c(:,8,:) = (-1 / (2 * B0) * (gy * 1e-2) .* (gz * 1e-2));        end % yz
if Nl >= 9, c(:,9,:) = (-1 / (2 * B0) * (gx * 1e-2) .* (gz * 1e-2));        end % xz

%% Calculate the coefficients of concomitant fields [T/m^3]
%--------------------------------------------------------------------------
% 1/B0^2 order concomitant gradient terms (10 terms)
% (x3, y3, z3, x2y, x2z, xy2, y2z, xz2, yz2, xyz)
%--------------------------------------------------------------------------
if Nl >= 10, c(:,10,:) = (-1 / (8 * B0^2) * (gx * 1e-2) .* (gz * 1e-2).^2);                                                      end % x3
if Nl >= 11, c(:,11,:) = (-1 / (8 * B0^2) * (gy * 1e-2) .* (gz * 1e-2).^2);                                                      end % y3
if Nl >= 12, c(:,12,:) = (-1 / (2 * B0^2) * (gz * 1e-2) .* ((gx * 1e-2).^2 + (gy * 1e-2).^2));                                   end % z3
if Nl >= 13, c(:,13,:) = (-1 / (8 * B0^2) * (gy * 1e-2) .* (gz * 1e-2).^2);                                                      end % x2y
if Nl >= 14, c(:,14,:) = (-1 / (2 * B0^2) * (1 / 4 * (gz * 1e-2).^3 - (gx * 1e-2).^2 .* (gz * 1e-2)));                           end % x2z
if Nl >= 15, c(:,15,:) = (-1 / (8 * B0^2) * (gx * 1e-2) .* (gz * 1e-2).^2);                                                      end % xy2
if Nl >= 16, c(:,16,:) = (-1 / (2 * B0^2) * (1 / 4 * (gz * 1e-2).^3 - (gy * 1e-2).^2 .* (gz * 1e-2)));                           end % y2z
if Nl >= 17, c(:,17,:) = (-1 / (2 * B0^2) * ((gx * 1e-2) .* ((gx * 1e-2).^2 + (gy * 1e-2).^2) - (gx * 1e-2) .* (gz * 1e-2).^2)); end % xz2
if Nl >= 18, c(:,18,:) = (-1 / (2 * B0^2) * ((gy * 1e-2) .* ((gx * 1e-2).^2 + (gy * 1e-2).^2) - (gy * 1e-2) .* (gz * 1e-2).^2)); end % yz2
if Nl >= 19, c(:,19,:) = (1  / (B0^2)     * (gx * 1e-2) .* (gy * 1e-2) .* (gz * 1e-2));                                          end % xyz

%% Numerically integrate the coefficients
%--------------------------------------------------------------------------
% [rad/sec/T] * [T/m] * [sec] => [rad/m]
%--------------------------------------------------------------------------
k = cumsum(gamma * c * dt); % [rad/m], [rad/m^2], [rad/m^3]

end
