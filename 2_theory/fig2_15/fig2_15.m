%Illustration of a raised-cosine shaped tapering window

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

%% Parameters for SFS Toolbox
conf = struct();

%sound field synthesis
conf.usetapwin = true;
conf.tapwinlen = 0.3;

%secondary sources
conf.secondary_sources.number = 401;
conf.secondary_sources.size = 4;
conf.secondary_sources.center = [0 0 0];
conf.secondary_sources.geometry = 'line';
conf.secondary_sources.logspread = 1.0;
x0 = secondary_source_positions(conf);

%% Generate and plot tapering window
win = tapering_window(x0,conf);

colours = moreland(2);

figure
    plot(x0(:,1),win,'Color',colours(1,:))
    grid
    xlabel('x0 / m'), ylabel('tapering window weight')
    ylim([0 1.2])

savefig('fig2_15')