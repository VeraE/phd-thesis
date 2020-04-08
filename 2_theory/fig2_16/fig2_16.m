%Impulse responses of a real and a virtual point source in 2.5D WFS in free field
%and inside a listening room

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

%% Load impulse responses from localisation experiment
load(['Real_source_free_field.mat'])
load(['Real_source_room.mat'])
load(['WFS_free_field.mat'])
load(['WFS_room.mat'])

%zeropad shorter irs for plots
max_N = max([length(Real_IR) length(Real_RIR) length(WFS_IR) length(WFS_RIR)]);
Real_IR = [Real_IR; zeros(max_N-length(Real_IR),1)];
Real_RIR = [Real_RIR; zeros(max_N-length(Real_RIR),1)];
WFS_IR = [WFS_IR; zeros(max_N-length(WFS_IR),1)];
WFS_RIR = [WFS_RIR; zeros(max_N-length(WFS_RIR),1)];

%% Plots
colours = moreland(2);

%time vector
fs = 44100;
t = (0:max_N-1)/fs*1000;

figure
    plot(t,Real_IR,'Color',colours(1,:))
    grid
    axis([0 80 -.01 .03])
    xlabel('time / ms'), ylabel('amplitude')
    set(gca,'YTick',-.01:.01:.03)

savefig('fig2_16a')

figure
    plot(t,WFS_IR,'Color',colours(2,:))
    grid
    axis([0 80 -.01 .03])
    xlabel('time / ms'), ylabel('amplitude')
    set(gca,'YTick',-.01:.01:.03)

savefig('fig2_16c')

figure
    plot(t,Real_RIR,'Color',colours(1,:))
    grid
    axis([0 80 -.01 .03])
    xlabel('time / ms'), ylabel('amplitude')
    set(gca,'YTick',-.01:.01:.03)

savefig('fig2_16b')

figure
    plot(t,WFS_RIR,'Color',colours(2,:))
    grid
    axis([0 80 -.01 .03])
    xlabel('time / ms'), ylabel('amplitude')
    set(gca,'YTick',-.01:.01:.03)

savefig('fig2_16d')