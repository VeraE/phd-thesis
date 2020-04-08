%Mean HpTF of the FABIAN HATS with headphones type AKG K601, compensation filter
%and compensated result
%
%To run this script, you need to download the HpIRs from:
%http://doi.org/10.14279/depositonce-5718.3
%
%The calculation of the headphone compensation filter is based on the code
%provided at: http://doi.org/10.5281/zenodo.401042

%*****************************************************************************
% Copyright (c) 2020      Vera Erbes                                         *
%                         Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock, Germany  *
%                                                                            *
% This file is part of the supplementary material for Vera Erbes' PhD        *
% thesis.                                                                    *
%                                                                            *
% You can redistribute the material and/or modify it under the terms of the  *
% GNU General Public License as published by the Free Software Foundation,   *
% either version 3 of the License, or (at your option) any later version.    *
%                                                                            *
% This Material is distributed in the hope that it will be useful, but       *
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
% or FITNESS FOR A PARTICULAR PURPOSE.                                       *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy of the GNU General Public License along   *
% with this program. If not, see <http://www.gnu.org/licenses/>.             *
%                                                                            *
%                                                  vera.erbes@uni-rostock.de *
%*****************************************************************************

close all; clear; clc

SFS_start
SOFAstart

%% Choose parameters
clen = 2048; %target length of compensation filter
wlen = 1024; %length of window for raw IR measurements
%regularisation
fc = 6000; %upper corner frequency for regularisation high-shelf
beta = .4; %regularisation effort
%1: same filter for left/right ear averaged, 0: different filter for left/right
one_filter = 1;

%% Load headphone measurements
N = 12; %number of measurements per channel

%select measurements of left and right channel (outliers removed)
idxl = [1 2 3 4 5 6 7 8   10 11   ];
idxr = [1 2 3 4 5 6   8 9 10 11 12];

%load all
sofa = SOFAload('HpIRs.sofa');
hl_raw = squeeze(sofa.Data.IR(:,1,:))';
hr_raw = squeeze(sofa.Data.IR(:,2,:))';
fs = sofa.Data.SamplingRate; %sampling frequency

%% Window and select IRs
%normalise all IRs
for n = 0:N/12-1
    idx = (1:12)+n*12;
    maxvalue(n+1) = max(max(abs([hl_raw(:,idx); hr_raw(:,idx)])));
    hl_raw(:,idx) = hl_raw(:,idx)./maxvalue(n+1);
    hr_raw(:,idx) = hr_raw(:,idx)./maxvalue(n+1);
end

%window IRs and truncate to target length
win = blackmanharris(wlen);
win(1:wlen/2) = 1;
win = [win; zeros(clen-wlen,1)];
hl_win = zeros(clen,N); %preallocation
hr_win = zeros(clen,N); %preallocation
for n = 1:N
    hl_win(:,n) = hl_raw(1:clen,n).*win;
    hr_win(:,n) = hr_raw(1:clen,n).*win;
end

%pick only selected windowed IRs
hl = hl_win(:,idxl);
hr = hr_win(:,idxr);

%% Calculate complex mean of selected HpTFs
%if only one filter for both ears: concatenate IRs of left and right channel
if one_filter
    hl = [hl hr];
    hr = hl;
end

%normalise selected IRs
maxvalue(end+1) = max(abs([hl(:); hr(:)]));
hl = hl./maxvalue(end);
hr = hr./maxvalue(end);

%FFT of selected IRs
Hl = fft(hl,[],1);
Hr = fft(hr,[],1);

%complex mean
Hml = sum(Hl,2)/size(Hl,2);
Hmr = sum(Hr,2)/size(Hr,2);

%% Design target bandpass
fl = 20; %lower corner frequency
fh = 20000; %upper corner frequency
stopatt = 60; %stopband attenuation

%calculate beta for kaiser window to satisfy desired stopband attenuation
beta_kaiser = .1102*(stopatt-8.7);
w = kaiser(clen+1,beta_kaiser);
b = fir1(clen,[fl/(fs/2) fh/(fs/2)],w);
b = b(2:end);
D = fft(b);
D = D.';

%% Design regularisation filter (high-shelf)
freg = [0 2000 fc fs/2]/(fs/2); %specified frequencies
G = [-20 -20 0 0]; %gain in dB at specified frequencies
g = 10.^(G/20); %linear gain at specified frequencies
b = fir2(50,freg,g);
b = [b'; zeros(clen-length(b),1)];
B = fft(b);

%% Calculate inverse filter in frequency domain
Hcl = D.*conj(Hml)./(Hml.*conj(Hml)+beta*B.*conj(B));
Hcr = D.*conj(Hmr)./(Hmr.*conj(Hmr)+beta*B.*conj(B));

hcl = ifft(Hcl);
hcr = ifft(Hcr);

%% Calculate compensation results
%compensation result for mean HpTF
hml_comp = conv(hcl,ifft(Hml));
hmr_comp = conv(hcr,ifft(Hmr));

%% Plot for diss
%colour for plot
colours = moreland(4);

%freq. vector for shortened IRs
f = (0:clen-1)/clen*fs;
%freq. vector for comp. results
fcomp = (0:length(hml_comp)-1)/length(hml_comp)*fs;

figure
    semilogx(f,db(abs(Hml)),'Color',colours(1,:)), hold on
    semilogx(f,db(abs(Hcl)),'Color',colours(3,:))
    semilogx(fcomp,db(abs(fft(hml_comp))),'Color',colours(4,:))
    grid
    xlabel('frequency / Hz'), ylabel('magnitude / dB')
    axis([20 fs/2 -30 20])
    set(gca,'XTickLabel',{'100','1k','10k'})
    legend('averaged HpTF','compensation filter','compensated result',...
        'Location','SW')

savefig('fig3_02')