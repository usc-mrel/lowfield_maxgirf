function [kPred, GPred] = apply_GIRF(gradients_nominal, dt, R, tRR)
% function [kPred, GPred] = apply_GIRF(gradients_nominal, dt, R, tRR)
%
% % UNITS
% gradients_nominal [ G/cm ]
% dt                [ s ]
%
% Hack to R to handle field strength (backwards compatibility)
% R.R = rotation matrix;
% R.T = field strength {}
%
% tRR is sub-dwell-time offset (-1 to 1) [0]?


% handle "nasty" co-opting of R-variable to include field info.
if isstruct(R)
    field_T = R.T;
    R = R.R;
else
    field_T = 0.55;
end

%% LOAD GIRF (Field/scanner dependent)
% Field selection needs to be handled upstream, or include the headers
% in to this file

if field_T == 1.4940
    % 1.5T Aera (NHLBI 2016)
    girf_file = 'GIRF_20160501.mat';
elseif field_T == 0.55
    % 0.55T Aera (NHLBI 2018)
    girf_file = 'GIRF_20200221_Duyn_method_coil2.mat';
end

% === Load file ===
try
    load(girf_file);
    disp(['Using ' girf_file]);
    %GIRF(:,2) = GIRF(:,1); % Hack by NGL
catch
    disp('Couldn''t find the GIRF file, using ones..');
    GIRF = ones(3800,3);
end

%%
dtGIRF = 10e-6;
dtSim = dt;
l_GIRF = length(GIRF);
[samples, interleaves, gs] = size(gradients_nominal);

% if readout is real long, need to pad the GIRF measurement
if samples*dt > dtGIRF*l_GIRF
    disp('readout length > 38ms, zero-padding GIRF');
    pad_factor = 1.5 * (samples * dt) / (dtGIRF * l_GIRF); % 1.5 factor to ease calculations below
    new_GIRF = zeros(round(l_GIRF * pad_factor), 3);

    for i = 1:3
        fft_GIRF = fftshift(ifft(ifftshift(GIRF(:,i))));
        zeropad = round( abs( (l_GIRF-length(new_GIRF) ) /2 ));
        temp = zeros(length(new_GIRF),1);
        % smoothing of padding:
        H_size = 200;
        H = hanningt(H_size);
        fft_GIRF(1:(H_size/2)) = fft_GIRF(1:(H_size/2)).*reshape(H(1:(H_size/2)),size(fft_GIRF(1:(H_size/2))));
        fft_GIRF(end-(H_size/2 - 1):end) = fft_GIRF(end-(H_size/2 - 1):end).*reshape(H((H_size/2 + 1):H_size),size(fft_GIRF(end-(H_size/2 - 1):end)) );

        temp( (1+zeropad):(zeropad + l_GIRF) )= fft_GIRF;
        % figure, plot(real(temp))
        new_GIRF(:,i) = fftshift(fft(fftshift(temp)));
    end

    old_GIRF = GIRF;
    GIRF = new_GIRF;
    l_GIRF = length(GIRF);

    %     new_GIRF = zeros(round((samples*dt) / (dtGIRF)),3);
    %
    %     for i = 1:3
    %         fft_GIRF = (ifft((GIRF(:,i) ))); % figure, plot(abs(fft_GIRF))
    %
    %         temp = zeros(length(new_GIRF),1);
    % %         H_size = 2000;        H = hanningt(H_size); %figure, plot(H)
    % %         fft_GIRF(end-(H_size/2 - 1):end) = fft_GIRF(end-(H_size/2 - 1):end).*reshape(H((H_size/2 + 1):H_size),size(fft_GIRF(end-(H_size/2 - 1):end)) );
    % %
    %         temp(1:length(fft_GIRF)) = fft_GIRF; % hold on, plot(abs(temp))
    %
    %         new_GIRF(:,i) = (fft((temp)));
    %     end
end
% %%
%     figure,
%     subplot(2,2,1), plot(abs(GIRF)), title('mag G')
%     subplot(2,2,2), plot(abs(new_GIRF)),title('mag G''')
%     subplot(2,2,3), plot(angle(GIRF)), title('angle G')
%     subplot(2,2,4), plot(angle(new_GIRF)),title('angle G''')

%%   GIRF prediction   %
%%%%%%%%%%%%%%%%%%%%%%%
% tRR = -0.7;
% tRR = 0;
% tRR = tRR + .1746;
% tRR = tRR + .3492;
if nargin < 4
    tRR = 0;
end

ADCshift = 0.85e-6 + 0.5 * dt + tRR * dt; %NCO Clock shift

%%
clear G0 GNom GPred kNom kPred

% allocation
Nominal   = zeros(samples, 3, 'double');
Predicted = zeros(samples, 3, 'double');
GNom      = zeros(samples, 3, interleaves, 'double');
GPred     = zeros(samples, 3, interleaves, 'double');
kNom      = zeros(samples, 3, interleaves, 'double');
kPred     = zeros(samples, 3, interleaves, 'double');

% GIRF process
for l = 1:interleaves
    % Select current spiral arm
    G0(:,1) = gradients_nominal(:,l,1);
    G0(:,2) = gradients_nominal(:,l,2);
    G0(:,3) = zeros(length(gradients_nominal),1);
    %Rotate into physical coordinates
    G0 = (R * G0.').';

    %--Loop through x,y,z gradient trajectories--%
    for ax = 1:3

        % Zeropad in time domain to match frequency resolution of GIRF (match readout length)
        L = round(dtGIRF * l_GIRF / dtSim); % when waveform not at GRT
        G = zeros(L,1, 'double');

        % REQUIRES SOME INTELLIGENCE TO DETERMINE WHICH SPIRAL FORM

        %------------------------------------------------------------------
        % Modified by NGL -- start
        %------------------------------------------------------------------
        % Original code:
        % SPIRAL OUT
        %G(1:length(G0)) = G0(:,ax); % figure, plot(G)
        % Make waveform periodic by returning to zero
        %          H = G(end)*hanningt((200)*2); % figure, plot(H)
        %          G((length(G0)+1):(length(G0)+length(H))) = H;
        %H = G(length(G0))*hanningt((200)*2); % figure, plot(H)
        %G((length(G0)+1):(length(G0)+0.5*length(H))) = H(length(H)*0.5 + 1:end); % figure, plot(G)
        %
        % RR shift to middle to avoid offset?
        %G = circshift(G,round(L/2 - samples/2)); % figure, plot(G)
        %------------------------------------------------------------------
        % SPIRAL OUT
        N = length(G0);
        index_range = (-floor(N/2):ceil(N/2)-1).' + floor(L/2) + 1;
        G(index_range) = G0(:,ax);

        % Make a waveform periodic by returning to zero
        H = G(index_range(N)) * hanningt(400);
        G(index_range(N)+(1:length(H)*0.5)) = H(length(H)*0.5+1:end);
        %------------------------------------------------------------------
        % Modified by NGL -- end
        %------------------------------------------------------------------

        % SPIRAL IN-OUT
        %          G((1:length(G0))) = G0(:,ax); % pad first point to 0

        %FFT nominal gradient
        dw = 1 / (L * dtSim); % frequency resolution [Hz]
        %w = 0:dw:dw*(length(G)-1); % Modified by NGL
        w = (-floor(L/2):ceil(L/2)-1).' * dw; % [Hz]
        %I = fftshift(fft(fftshift(G))); % Modified by NGL
        I = fftshift(fft(ifftshift(G)));

        %Zeropad GIRF and I to bandwidth of sampling (when waveform not at GRT)
        GIRF1 = zeros(L,1, 'double');

        if dt > dtGIRF
            % RR crop
            GIRF1 = GIRF(round(l_GIRF/2 - L/2 + 1):round(l_GIRF/2 + L/2),ax);
            % RR .. padding
            temp = hanningt(10);
            GIRF1(1) = 0; GIRF1(end) = 0;
            GIRF1(2:round(length(temp)/2) + 1) = GIRF1(2:round(length(temp)/2)+1).*reshape(temp(1:round(length(temp)/2)),size(GIRF1(2:round(length(temp)/2)+1)));
            GIRF1(end-round(length(temp)/2):end-1) = GIRF1(end-round(length(temp)/2):end-1).*reshape(temp((round(length(temp)/2) + 1):length(temp)),size(GIRF1(end-round(length(temp)/2):end-1)));
            %             figure, plot(real(GIRF1))
        else
            %--------------------------------------------------------------
            % Modified by NGL
            % Usual operation.. (ACW)
            %zeropad = round(abs((l_GIRF-L)/2)); %amount of zeropadding
            %GIRF1((1+zeropad):(zeropad+l_GIRF))= GIRF(:,ax);
            %--------------------------------------------------------------
            index_range = (-floor(l_GIRF/2):ceil(l_GIRF/2)-1).' + floor(L/2) + 1;
            GIRF1(index_range) = GIRF(:,ax);
        end

        %Predict Gradient and apply clock shift
        %------------------------------------------------------------------
        % Modified by NGL
        %------------------------------------------------------------------
        %P = I.*GIRF1.*exp(-1i.*ADCshift*2*pi*w)'; % P = I.*GIRF1.*exp(1i.*ADCshift*2*pi*w)';
        P = I .* GIRF1 .* exp(1j * ADCshift * 2 * pi * w);

        %zeropad to required bandwidth (when waveform is at GRT)
        BW = 1 / dt;
        L = round(BW/dw);    %number of points required

        PredGrad = zeros(L,1, 'double');
        NomGrad = zeros(L,1, 'double');
        zeropad = round(abs((length(G)-L)/2)); %amount of zeropadding

        PredGrad((1+zeropad):(zeropad+length(P)))= P;
        NomGrad((1+zeropad):(zeropad+length(I)))= I;

        %------------------------------------------------------------------
        % Modified by NGL
        % PredGrad = ifftshift(ifft(ifftshift(PredGrad)));
        % NomGrad  = ifftshift(ifft(ifftshift(NomGrad)));
        %------------------------------------------------------------------
        %FFT back to time domain
        PredGrad = fftshift(ifft(ifftshift(PredGrad))) * L / length(G);
        NomGrad  = fftshift(ifft(ifftshift(NomGrad)))  * L / length(G);

        %Correct polarity of gradients
        multiplier = zeros(length(PredGrad),1);
        for i = 1:length(PredGrad)
            if real(PredGrad(i))>0; multiplier(i) = 1;
            else; multiplier(i) = -1;
            end
        end
        PredGrad = abs(PredGrad).*multiplier;

        multiplier = zeros(length(NomGrad),1);
        for i = 1:length(NomGrad)
            if real(NomGrad(i))>0; multiplier(i) = 1;
            else; multiplier(i) = -1;
            end
        end
        NomGrad = abs(NomGrad).*multiplier;

        %------------------------------------------------------------------
        % Modified by NGL
        % RR - circle back
        %NomGrad = circshift(NomGrad, -round(L/2 - samples/2) -1);
        %PredGrad = circshift(PredGrad, -round(L/2 - samples/2) -1);
        %
        %Only take the samples relevant to the readout
        %Nominal(:,ax) = NomGrad(1:samples);
        %Predicted(:,ax) = PredGrad(1:samples);
        %------------------------------------------------------------------
        index_range = (-floor(samples/2):ceil(samples/2)-1).' + floor(L/2) + 1;
        Nominal(:,ax) = NomGrad(index_range);
        Predicted(:,ax) = PredGrad(index_range);
    end

    %rotate back to logical coordinates
    GNom(:,:,l)= (R.'*Nominal.').';
    GPred(:,:,l) = (R.'*Predicted.').';

    %Integrate to get k-space trajectory from gradient
    kNom(:,:,l)  = cumsum(GNom(:,:,l));
    kPred(:,:,l) = cumsum(GPred(:,:,l));
end

% Scale k-space in units of 1/cm
gamma = 2.67522212e+8; % gyromagnetic ratio for 1H [rad/sec/T]
% Modified by NGL: 2.675e8 => gamma
kPred = 0.01*(gamma/(2*pi))*(kPred*0.01)*dt; % (kPred*0.01):: assuming gradients in are in G/cm!!!
kNom  = 0.01*(gamma/(2*pi))*(kNom*0.01)*dt; % (kPred*0.01):: assuming gradients in are in G/cm!!!

%Permute
kPred = permute(kPred,[1 3 2]);
kNom  = permute(kNom,[1 3 2]);
GPred = permute(GPred,[1 3 2]);
GNom  = permute(GNom,[1 3 2]);

% figure,
% subplot(2,2,1); plot(kNom(:,:,1),kPred(:,:,1),'-')
% subplot(2,2,2); plot(kNom(:,:,2),kPred(:,:,2),'-')
% subplot(2,2,3); plot(GNom(:,:,1),GPred(:,:,1),'-')
% subplot(2,2,4); plot(GNom(:,:,2),GPred(:,:,2),'-')

%%%%%%%%%%%%%%%%%
%%  calculate gridder  %
%%%%%%%%%%%%%%%%%
%
% omega = kPred(:,:,1:2)*(pi/max(max(max(kPred))));
% omega = reshape(omega,interleaves*samples,2);
%
% gradients_nominal = reshape(GPred(:,:,1:2),interleaves*samples,2);
% grad = complex(gradients_nominal(:,1),gradients_nominal(:,2));
% kk = complex(omega(:,1),omega(:,2));
% weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
%
% st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);


end

function wind = hanningt(windowLength)
% If license toolbox being an idiot..
% ripped from: https://www.mathworks.com/matlabcentral/fileexchange/48925-hann-window

if license('checkout','Signal_Toolbox')
    wind = hanning(windowLength);
else
    N = windowLength - 1;
    num = linspace(0,N,windowLength);
    wind =  0.5*(1 - cos(2*pi*num/N));
end

% matlab:
% function w = sym_hanning(n)
% %SYM_HANNING   Symmetric Hanning window.
% %   SYM_HANNING Returns an exactly symmetric N point window by evaluating
% %   the first half and then flipping the same samples over the other half.
%
% if ~rem(n,2)
%    % Even length window
%    half = n/2;
%    w = calc_hanning(half,n);
%    w = [w; w(end:-1:1)];
% else
%    % Odd length window
%    half = (n+1)/2;
%    w = calc_hanning(half,n);
%    w = [w; w(end-1:-1:1)];
% end
%
% %---------------------------------------------------------------------
% function w = calc_hanning(m,n)
% %CALC_HANNING   Calculates Hanning window samples.
% %   CALC_HANNING Calculates and returns the first M points of an N point
% %   Hanning window.
%
% w = .5*(1 - cos(2*pi*(1:m)'/(n+1)));

end