function [Psi,inv_L] = process_noise_ismrmrd_data(noise_filename)


%% Read an ismrmrd file
tic; fprintf('Reading an ismrmrd file: %s... ', noise_filename);
if exist(noise_filename, 'file')
    dset = ismrmrd.Dataset(noise_filename, 'dataset');
    fprintf('done! (%6.4f sec)\n', toc);
else
    error('File %s does not exist.  Please generate it.' , noise_filename);
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
