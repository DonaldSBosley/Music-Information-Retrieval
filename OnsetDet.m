%beatdet
%
%wavin = the input wave file to be read
%winLength = window length in sample
%winType: window type: should pass one of the following strings
%           - ‘rect’
%           - ‘hamming’
%           - ‘hann’
%           - ‘blackman’
%           - ‘bartlett’

function DonBosleybeatdet(wavin, winLength, winType)
%% READ IN WAVE FILE - CHECK VALIDITY

[isitawave] = wavfinfo (wavin);         %get file information

switch isitawave                        %check for a valid wave
    case 'Sound (WAV) file'             %valid wave
        disp (['Successfully opened ', wavin '.wav'])
    case ''                             %invalid wave + error message
        wavin = input('Please input a valid wavefile:'); 
end

[wavin, fs, nbits] = wavread (wavin);   %Read in wavfile
 
%% NUMBER OF CHANNELS FOR SIGNAL

[~ , wavecolumns] = size(wavin);        %get wave matrix dimensions
if wavecolumns == 1                     %it's a mono file,  
elseif wavecolumns == 2                 %the file is stereo
        wavin = mean (wavin, 2);        %sum down to mono
elseif wavecolumns >= 3                 %format isn't valid
        error ('You appear to have a wave file that is not mono or stereo.')
end

%% WINDOW TYPE 
switch winType 
    case 'rectwin'
        window = rectwin(winLength); %rectangular window
    case 'hamming'
        window = hamming(winLength);% hamming window
    case 'hann'
        window = hann(winLength);% hann window
    case 'blackman'
        window = blackman(winLength);%blackman window
    case 'bartlett'
        window = bartlett(winLength); %bartlett window
end

%% USE SPECTRAL FLUX
%% Global Variables for Frequency Domain
olapLen = winLength/4;
hop = winLength - olapLen;
%This aligns the peaks properly
zerowav = [zeros(winLength/2, 1); wavin];       %Add 1/2 window zeropad
specout = spectrogram(zerowav, window, olapLen);%Take spectrogram
[N, M] = size(specout);                         %Size
timef = (1:M-1)* hop/fs;                        %Time Vector for plotting

%% Calculate Spectral Flux Half-Wave
specdifabsx = diff(abs(specout), 1, 2);         %1:Nyquist differences
specdifx = diff(specout, 1, 2);                 %Non-abs differences
specHx = (specdifabsx + specdifx) / 2;          %Final calc
specmeans = mean(specHx(:, 1:M-1));             %Mean energy across all freqs
% Bring down average (psuedo normalize throught averaging)
buffcompare = abs(mean(buffer(specmeans, 10,9)));%Find local means
newhalfwav = specmeans - buffcompare;           %Subtract local means
compmat = zeros(1, length(specmeans));          %Vector of zeros
halfwavspec = gt (newhalfwav, compmat);         %Greater than zero?
avghalfwav = abs(newhalfwav .* halfwavspec);    %Multiply < 0 by 0
avghalfwav = avghalfwav ./ max(avghalfwav);     %Normalize

%% Remove Local Mean + Half Wav Rec
buffcompare = mean(buffer(avghalfwav, 10,9));
NOV = avghalfwav - buffcompare;
compmat = zeros(1, length(NOV));                %Vector of zeros
comp = gt (NOV, compmat);         %Greater than zero?
Nov = abs(NOV.* comp);    %Multiply < 0 by 0
%plot(Nov)

%% Create Initial Novelty Matrix
%buffNov = buffer(Nov, winLength/2, winLength/4);  %Place novelty into buf
buffNov = buffer(Nov, winLength, winLength/4);
[L, R] = size(buffNov);                             %size
pad = zeros(L,R);                                   %create zero pad
buffNovpad = [buffNov; pad];                        %concatenate pad
[Q, P] = size(buffNovpad);                          %find new size

%% TAKE ACF OF NOVELTY FUNCTION

novACF1 = real( ifft( (abs (fft(buffNovpad) ) .^2 ) ));%take autocorrelation of novelty
novACF = (novACF1((1:L), :));                 %elminate above Nyquist
lagsvec = (L: -1 : 1)';            %THIS MIGHT BE BACKWARDS
lagsmat = repmat(lagsvec,1, P);
newNovMat = novACF ./ lagsmat ;                  %Unbiased novelty matrix

%% CREATE COMB FILTER MATRIX 
%{
L = length(newNovMat);

comb = zeros(L, 4*L); % figure out what horizontal size is - must correspond to ACF size

b = 50;
for lag = 1:L % maybe nos start in lag 1, calculate min BPM
   R_w = (lag/(b^2)) * exp(-lag^2/(2*b^2)); % R*** weighting function
   comb(lag, lag:lag:lag*nFB) = R_w; % comb filter bank - again with a specific hor. size..
end
%}

%% COMB FILTER MATRIX (ala TAE MIN)
nFB = 4;                    %Number of pulses per filter
Lmax = ceil(fs/hop);    %find min tempo/max lags
Lmin = ceil(fs/4/hop);
b = 43;                     %part of distribution function; responds to 120 bpm?
comb = zeros(Lmax - Lmin+1, L);

for k = 1:Lmax-Lmin+1
    lag = k + Lmin - 1;
    weight = (lag/(b^2)) * exp (-lag^2 / (2*b^2));
    comb(k, lag:lag:lag*nFB) = weight;
end

%imagesc(comb)
%% ORIGINAL FILTER BANK
%{
bpmforACF = zeros(length(newNovMat));           %preallocate tempo matrix

for k = (1: Q/2)
    pulseholder = zeros((length(newNovMat)),1);
    pulseholder(1:k:end) = 1;
    pulseholder = pulseholder ./ sum(pulseholder);
    bpmforACF(:,k) = pulseholder;
end
%imagesc(bpmforACF)
%}

%% NOVELTY MATRIX
novMat = comb * newNovMat; 
[~, index] = max(novMat);
lags = index + Lmin - 1;

%% Global average tempo - to check output of function up to this point
avgidx = mean(lags);                   %average bin for max lag values
charles = length(newNovMat);
avgbps = (avgidx * charles)/fs;             %convert to beats per second
avgbpm = 60/avgbps;                     %convert to beats per minute

%% Local tempi in each window
localbps = lags .* hop ./fs;
localbpm = 60 ./ localbps;

%% PLOTS:

subplot (3,1,1)
plot(Nov)                                       %novelty function
subplot(3,1,2:3)
imagesc(novMat)                                 %tempo gram
hold on
plot(localbpm)
hold off


end

