%Illustration of a comb filter

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

%% Create signal
N = 2^16;
signal = zeros(1,N);
signal(1) = 1;
signal(10) = 1;

%calculate spectrum
Signal = fft(signal);
%frequency vector
fs = 44100; %sampling rate
f = (0:N-1)/N*fs;

%% Plot
colours = moreland(2);

figure
    plot(f,abs(Signal),'Color',colours(1,:))
    grid
    xlabel('frequency / Hz'), ylabel('|X(f)|')
    axis([0 fs/2 0 2])
    set(gca,'XTickLabel',{'0','5000','10000','15000','20000'})

savefig('fig2_17')