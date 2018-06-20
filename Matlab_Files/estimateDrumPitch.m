function [  ] = estimateDrumPitch(path)

% read the audio file
[signal, fs] = audioread(path);

% compute length and create window
n = length(signal);
hannWindow = hann(n,'periodic');

% --- PRE PROCESSING ---

% down-mixing
if (size(signal,2)> 1)
    signal = mean(signal,2);
end
% normalization
signal = signal/max(abs(signal));

% plot signal
subplot(3,1,1);
plot(signal);
title('Input Signal');
xlabel('Sample');
ylabel('Normalized Amplitude');

% -----------------------

% --- COMPUTE MAGNITUDE SPECTRUM ---

% create spectrum using dft instead of spectrogram func
windowedSignal = signal;%.* hannWindow;

% Because the signal is real-valued, you only need power estimates for the positive or negative frequencies. 
% In order to conserve the total power, multiply all frequencies that occur in both sets by a factor of 2.
magnitudeSpectrum = abs(fft(windowedSignal))*2/n; % /n because of DFT
magnitudeSpectrum = magnitudeSpectrum(1:n/2+1); % cut negative half as we dont need it

% plot spectrum
subplot(3,1,2);
plot (10*log10(magnitudeSpectrum));

set(gca,'xscale','log');
xlim([10 21000]);

title('Magnitude Spectrum');
xlabel('Frequency[Hz]');
ylabel('Magnitude[dB]');

% -----------------------

% --- COMPUTE ACF AND FREQUENCY ---

% initialize acf variables
f_min = 100; % minimum frequency to be detected (threshold)

eta_min = round(f_min/fs * n); % Convert min frequency to index in acf

% allocate (not needed yet because we just have a single block)
f = zeros(1, size(magnitudeSpectrum,2));

% use spectral symmetry for robustness
% what this does is to add a mirrored version of the signal to make it
% symmetric
magnitudeSpectrum = [flipud(magnitudeSpectrum); magnitudeSpectrum];

% compute the ACF (for each block) --> size(X,2) = 1
for (i = 1: size(magnitudeSpectrum,2))
    afCorr          = xcorr(magnitudeSpectrum(:,i),'coeff'); % Normalizes the sequence so that the autocorrelations at zero lag equal 1
    afCorr          = afCorr((ceil((length(afCorr)/2))+1):end); % Cut left half
    [fDummy, f(i)]  = max(afCorr(1+eta_min:end)); % Store row (index) of maximum value in f
end


% f is now our index with the highest peak
f = f + eta_min; % compensate lag introduced in max function

% plot acf
subplot(3,1,3)
plot (10*log10(afCorr));

xlim([0 8000])

title('ACF');
xlabel('Lag [Samples]');
ylabel('Correlation');

% convert max index to Hz
f = f / n * fs;
f % print result

% -----------------------


% ---- old stuff from testing ------

% compute spectrogram over whole file (block length = file length)
%[X, f, t] = spectrogram(signal, hannWindow, 0, n, fs);
% magnitude spectrum
%X = abs(X)*2/n; % why do we do that?
        
% compute frequency
%f           = hPitchFunc(X, f_s);


%[freq, time] = ComputePitch('SpectralAcf', signal, fs, [], length(signal), length(signal));

% pre-processing -> mono audio channels
%lv = ~sum(signal == 0, 2);                          % Rows Without Zeros (Logical Vector)
%mono = sum(signal, 2);                              % Sum Across Columns
%mono(lv) = mono(lv)/2;                              % Divide Rows Without Zeros By ‘2’ To Get Mean

%spectrogram = psd(spectrum.periodogram,signal,'Fs',fs,'NFFT',length(signal));



%plot(spectrogram);

%mag = abs(fft(signal, fs));
%mag = mag(1:21000);
%plot(mag2db(mag));

%set(gca, 'XScale', 'log');
%axis([10 21000 0 inf])
%xlabel('Frequency [Hz]');
%ylabel('Magnitude [dB]');
%title('Magnitude Spectrum of a drum hit');



%Fs = fs;


%x = signal;

%N = length(x);
%xdft = fft(x);
%xdft = xdft(1:N/2+1);
%psdx = (1/(Fs*N)) * abs(xdft).^2;
%psdx(2:end-1) = 2*psdx(2:end-1);
%freq = 0:Fs/length(x):Fs/2;

%plot(freq,10*log10(psdx))
%grid on
%set(gca, 'XScale', 'log');
%xlim([10 21000])
%xlabel('Frequency [Hz]');
%ylabel('Magnitude [dB]');
%title('Magnitude Spectrum of a drum hit');

%corr = xcorr(psdx,psdx);
%plot(10*log10(corr));



end

