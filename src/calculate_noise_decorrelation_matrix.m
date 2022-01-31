function [Psi,inv_L] = calculate_noise_decorrelation_matrix(noise_fullpath)
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/16/2022, Last modified: 01/16/2022

%% Read a noise only ISMRMRD file
tic; fprintf('Reading an ISMRMRD file: %s... ', noise_fullpath);
if exist(noise_fullpath, 'file')
    dset = ismrmrd.Dataset(noise_fullpath, 'dataset');
    fprintf('done! (%6.4f sec)\n', toc);
else
    error('File %s does not exist.  Please generate it.' , noise_fullpath);
end
raw_data = dset.readAcquisition();

%% Sort noise data
is_noise = raw_data.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
meas = raw_data.select(find(is_noise));
nr_repetitions = length(meas.data); % number of repetitions
[nr_samples,nr_channels] = size(meas.data{1});
noise = complex(zeros(nr_channels, nr_samples, nr_repetitions, 'double'));
for idx = 1:nr_repetitions
    noise(:,:,idx) = meas.data{idx}.'; % nr_samples x nr_channels => nr_channels x nr_samples
end

%% Calculate the receiver noise matrix
% Use the definition in Appendix B of Pruessmann et al. (MRM 46:638–651 (2001))
% ns denotes the number of noise samples taken per channel
% eta lists these samples in an nc x ns matrix
tic; fprintf('Calculating the receiver noise matrix... ');
Nc = nr_channels;
Ns = nr_samples * nr_repetitions;
eta = reshape(noise, [Nc Ns]);
Psi = eta * eta' / Ns; % Equation B1
fprintf('done! (%6.4f sec)\n', toc);

%% Calculate the Cholesky decomposition of the receiver noise matrix
tic; fprintf('Calculating the Cholesky decomposition of the receiver noise matrix... ');
L = chol(Psi, 'lower'); % Equation B4
inv_L = inv(L);
fprintf('done! (%6.4f sec)\n', toc);

end
