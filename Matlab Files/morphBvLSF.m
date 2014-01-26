function [y,yh,ys,yu] = morphBvLSF(bv,lv,fs,w,N,t,nH,minf0,maxf0,vuvError,...
    maxhd,stocf,p,fadeLen,expon,htintp,rintp)

% Author: Matthew Roddy, August 2013
% Builds upon code published in DAFX: Digital Audio Effects, 2nd ed.
% by J. Bonada, X. Serra, X. Amatriain, A. Loscos
%
% morph between two sounds using LSFs and the harmonic plus stochastic model
% 
% bv: backing vocal, lv: lead vocal 
% fs: sampling rate, w: analysis window (odd size),
% N: FFT size (minimum 512), t: threshold in negative dB,
% nH: maximum number of harmonics, minf0: minimum f0 frequency in Hz,
% maxf0: maximim f0 frequency in Hz,
% vuvError: error threshold in the f0 detection (ex: 0.2),
% maxhd: max. relative deviation in harmonic detection (ex: 5)
% stocf: decimation factor of mag spectrum for stochastic analysis
% p: linear prediction order (ex: 60)
% fadeLen: length of crossfade, in windows. Must be a whole number (ex: 3)
% expon: exponential factor for crossfade (ex: 2)
% htintp: harmonic timbre interpolation factor,
% rintp: residual interpolation factor,f0et
% yh: harmonic component, ys: stochastic component, yu: unvoiced component 
% y: output sound 
 
if length(htintp)==1
    htintp =[ 0      length(bv)/fs;    % input time
        htintp htintp       ];  % control value
end
if length(rintp)==1
    rintp =[ 0     length(bv)/fs;      % input time
        rintp rintp         ];   % control value
end
M = length(w);                        % analysis window size
Ns = 1024;                            % FFT size for synthesis
H = 256;                              % hop size for analysis and synthesis
N2 = N/2+1;                           % half-size of spectrum
soundlength = length(bv);              % length of input sound array
hNs = Ns/2;                           % half synthesis window size
hM = (M-1)/2;                         % half analysis window size
pin = max(H+1,1+hM);                  % initialize sound pointer to middle of analysis window
pend = soundlength-max(hM,H);         % last sample to start a frame
yh = zeros(soundlength+Ns/2,1);       % output sine component
ys = zeros(soundlength+Ns/2,1);       % output residual component
yu = zeros(soundlength+Ns/2,1);       % output residual component
w = w/sum(w);                         % normalize analysis window
sw = zeros(Ns,1);
ow = triang(2*H-1);                   % overlapping window
ovidx = Ns/2+1-H+1:Ns/2+H;            % overlap indexes
sw(ovidx) = ow(1:2*H-1);
bh = blackmanharris(Ns);              % synthesis window
bh = bh ./ sum(bh);                   % normalize synthesis window
wr = bh;                              % window for residual
sw(ovidx) = sw(ovidx) ./ bh(ovidx);
sws = H*hanning(Ns)/2;                % synthesis window for stochastic
minpin2 = max(H+1,1+hM);              % minimum sample value for x2
maxpin2 = min(length(lv)-hM-1);       % maximum sample value for x2
inFadePos = Ns;


while pin<pend

    %%%% First sound analysis %%%%
    [f0Bv,hlocBv,mXsenvBv,lsfBv,gainBv,hiBv,hphaseBv,bvVuV] = analysis(bv,fs,w,wr,pin,M,hM,N,N2,Ns,hNs,...
        nH,t,vuvError,minf0,maxf0,maxhd,stocf,p);
    %%%% Second sound analysis %%%%
    pin2 = round(pin/length(bv)*length(lv)); % linear time mapping between inputs
    pin2 = min(maxpin2,max(minpin2,pin2));
    [f0Lv,hlocLv,mXsenvLv,lsfLv,gainLv,hiLv,hphaseLv,lvVuV] = analysis(lv,fs,w,wr,pin2,M,hM,N,N2,Ns,hNs,...
        nH,t,vuvError,minf0,maxf0,maxhd,stocf,p);   

    
    %%%% Morph %%%%
    chtintp = interp1(htintp(1,:),htintp(2,:),pin/fs);   % get control value
    crintp = interp1(rintp(1,:),rintp(2,:),pin/fs);      % get control value
    mYsenv = mXsenvBv*(1-crintp) + mXsenvLv*crintp; % morph residual envelope
    
    %%%% LSFs %%%%
    ylsf = lsfBv*(1-chtintp) + lsfLv*chtintp; % interpolate the lsfs
    ya = lsf2poly(ylsf);  %Linear prediction coefficents
    A=-20*log10(abs(fft(ya,N2)));%  Magnitude spectrum from LP
    envelope=A+gainBv; % spectral envelope
    bins = 1:N2; 
    hMag = interp1(bins,envelope,hlocBv(1:hiBv-1)); % magnitudes of harmonics
    
    %%%%% Re-synthesis %%%%
    Yh = genspecsines(hlocBv(1:hiBv-1),hMag,hphaseBv,Ns);             % generate sines
    mYs = interp(mYsenv,stocf);                % interpolate to original size
    roffset = ceil(stocf/2)-1;                 % interpolated array offset
    mYs = [ mYs(1)*ones(roffset,1); mYs(1:Ns/2+1-roffset) ];
    mYs = 10.^(mYs/20);                           % dB to linear magnitude
    pYs = 2*pi*rand(Ns/2+1,1);                 % generate phase random values
    mYs1 = [mYs(1:Ns/2+1); mYs(Ns/2:-1:2)];    % create magnitude spectrum
    pYs1 = [pYs(1:Ns/2+1); -1*pYs(Ns/2:-1:2)]; % create phase spectrum
    Ys = mYs1.*cos(pYs1)+1i*mYs1.*sin(pYs1);   % compute complex spectrum
    yhw = fftshift(real(ifft(Yh)));            % sines in time domain using IFFT
    ysw = fftshift(real(ifft(Ys)));            % stoc. in time domain using IFFT
    ro = pin-hNs;                         % output sound pointer for overlap
    yh(ro:ro+Ns-1) = yh(ro:ro+Ns-1)+yhw(1:Ns).*sw;  % overlap-add for sines
    ys(ro:ro+Ns-1) = ys(ro:ro+Ns-1)+ysw(1:Ns).*sws; % overlap-add for stochastic
   
    [ synthFade,uvFade,outFadePos ] = bvFade( inFadePos,Ns,H,bvVuV,lvVuV,fadeLen,expon ); % crossfade algorithm
    inFadePos = outFadePos;    

    
    yh(pin-hM:pin-hM+H-1) = yh(pin-hM:pin-hM+H-1).*synthFade(:); %apply fade to hop size amount of samples
    ys(pin-hM:pin-hM+H-1) = ys(pin-hM:pin-hM+H-1).*synthFade(:);
    yu(pin-hM:pin-hM+H-1) = bv(pin-hM:pin-hM+H-1).*uvFade(:);
    pin = pin+H;                                    % advance the sound pointer

end
 y = yh+ys+yu;                                          % sum sines and stochastic

end

function [f0,hloc,mXsenv,lsf,gain,hi,hphase,yVuV] = analysis(x,fs,w,wr,pin,M,hM,N,N2,...
    Ns,hNs,nH,t,vuvError,minf0,maxf0,maxhd,stocf,p)

xw = x(pin-hM:pin+hM).*w(1:M);              % window the input sound
fftbuffer = zeros(N,1);                     % initialize buffer for FFT
fftbuffer(1:(M+1)/2) = xw((M+1)/2:M);       % zero-phase window in fftbuffer
fftbuffer(N-(M-1)/2+1:N) = xw(1:(M-1)/2);
X = fft(fftbuffer);                    % compute the FFT
mX = 20*log10(abs(X(1:N2)));           % magnitude spectrum
pX = unwrap(angle(X(1:N/2+1)));        % unwrapped phase spectrum
ploc = 1 + find((mX(2:N2-1)>t) .* (mX(2:N2-1)>mX(3:N2)) ...
    .* (mX(2:N2-1)>mX(1:N2-2)));          % find peaks
[ploc,pmag,pphase] = peakinterp(mX,pX,ploc);          % refine peak values
[f0, error] = TWM(ploc,pmag,N,fs,minf0,maxf0);    % find f0

if(error>vuvError)
    yVuV= 0;
else
    yVuV = 1;
end

hloc = zeros(nH,1);                    % initialize harmonic locations
hmag = zeros(nH,1)-100;                % initialize harmonic magnitudes
hphase = zeros(nH,1);                  % initialize harmonic phases
hf = (f0>0).*(f0.*(1:nH));             % initialize harmonic frequencies
hi = 1;                                % initialize harmonic index
npeaks = length(ploc);                 % number of peaks found

while (f0>0 && hi<=nH && hf(hi)<fs/2)  % find harmonic peaks
    [dev,pei] = min(abs((ploc(1:npeaks)-1)/N*fs-hf(hi)));    % closest peak
    if ((hi==1 || ~any(hloc(1:hi-1)==ploc(pei))) && dev<maxhd*hf(hi))
        hloc(hi) = ploc(pei);              % harmonic locations
        hmag(hi) = pmag(pei);              % harmonic magnitudes
        hphase(hi) = pphase(pei);          % harmonic phases
    end
    hi = hi+1;                           % increase harmonic index
end
hloc(1:hi-1) = (hloc(1:hi-1)~=0).*((hloc(1:hi-1)-1)*Ns/N);  % synth. locs

%%%%%%%%%%
% Calculate Power Spectral Density from harmonic Peaks %%%%%%%%%%
% First get rid of zeros in peak data and change harmonic index (hi) value
% Then add endpoints at zero and Fs/2, followed by interpolation.
I = find(hloc); %Get indices of nonzero elements
hi2 = length(I)+2; 
hloc2 = zeros(hi2,1);hloc2(1)=0;hloc2(2:hi2-1)=hloc(I);
hloc2(hi2) = N2;
hmag2 = zeros(hi2,1); hmag2(2:hi2-1) = hmag(I);
hmag2(1) = -150; hmag2(hi2) = -150; % Arbitrarily low magnitude for 0 Hz and fs/2 Hz
hPSD = power(10,hmag2(1:hi2).*((1/20)*2)); % First half of power spectrum
hFreq = hloc2(1:hi2).*(fs/N2);
xi = 1:(0.5*fs)/N2:0.5*fs;
yi = interp1(hFreq,hPSD,xi, 'linear'); % yi are squared absolute value
yPSD = [yi(1:N2-1) fliplr(yi(1:N2-1))]; % power spectral density


%%%%%%%%%%
% Linear Predicion Coefficients from Power Spectral Density
rxx = real(ifft(yPSD)); % Autocorrelation is IFFT of PSD
R = rxx(1:p);
if norm(R)~=0
    a=levinson(R,p);    % Levinson-Durbin recursion
else
    a=[1, zeros(1,p)];
end
lsf = poly2lsf(a);  % Calculate LSFs from LPC coefficients
R=R(:)'; a=a(:)';   % row vectors
g=sqrt(sum(a.*R));  % gain factor
gain=20*log10(g);   % Linear prediction gain in dB


%%%%%%%%%%
% Calculate Residual envelope by subtracting complex harmonic spectrum
% from original signal.

ri= pin-hNs;                     % input sound pointer for residual analysis
xr = x(ri:ri+Ns-1).*wr(1:Ns);          % window the input sound
Xr = fft(fftshift(xr));                % compute FFT for residual analysis
Xh = genspecsines(hloc,hmag,hphase,Ns);             % generate sine
Xr = Xr-Xh;                              % get the residual complex spectrum
mXr = 20*log10(abs(Xr(1:Ns/2+1)));       % magnitude spectrum of residual
mXsenv = decimate(max(-200,mXr),stocf);  % decimate the magnitude spectrum
end