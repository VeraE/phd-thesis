%Magnitude of the wave number spectrum of the 3D free-field Green's function

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

%% General definitions
f = 0:5:3000; %frequency / Hz
w = 2*pi*f;
kx = -80:.5:80;
c = 343; %speed of sound / m/s

yref = 1; %distance reference line / m
z = 0; %evaluate only z-plane

[kxkx,ww] = meshgrid(kx,w);

%% Calculate wave number spectrum
%initiate with all entries infinty
Gtilde = ones(length(w),length(kx))*inf;
%condition to discriminate between propagating and evenscent part
condition = abs(kxkx)-abs(ww)/c;
%propagating part
indH = find(condition<0);
Gtilde(indH) = -1i/4 ...
    *(  besselj(0,sqrt((ww(indH)/c).^2-kxkx(indH).^2)*yref)...
    -1i*bessely(0,sqrt((ww(indH)/c).^2-kxkx(indH).^2)*yref));
%evanescent part
indK = find(condition>0);
Gtilde(indK) = 1/(2*pi)...
    *besselk(0,sqrt(kxkx(indK).^2-(ww(indK)/c).^2)*yref);

%% Plot
figure
    imagesc(kx,f,20*log10(abs(Gtilde)))
    hold on
    %lines for |kx|=|w/c|
    plot(w/c,f,'k');
    plot(-w/c,f,'k')
    text(16,1300,'k_x = \omega /c','rotation',72)
    %colorbar
	colormap yellowred(256)
    set(gca,'CLim',[-100 0])
    c = colorbar;
    c.Ticks = -100:20:0;
    c.TickLabels = {'-100','-80','-60','-40','-20','0 dB'};
    %axes properties
    axis square
    axis([kx(1) kx(end) 0 f(end)])
    set(gca,'YDir','normal')
    set(gca,'XTick',-80:40:80)
    xlabel('k_x / rad/m'), ylabel('frequency / Hz')

savefig('fig2_06')