%Sound field of a virtual plane wave in 2.5D WFS without and with tapering window
%for the linear secondary source distribution

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

%misc
conf.tmpdir = '';
conf.showprogress = false;

%audio
conf.c = 343;

%sound field synthesis
conf.dimension = '2.5D';
conf.driving_functions = 'default';
conf.xref = [0 -1 0];
conf.tapwinlen = 0.3;

%sound field simulations
conf.resolution = 600;
conf.phase = 0;

%secondary sources
conf.secondary_sources.number = 401;
conf.secondary_sources.size = 4;
conf.secondary_sources.center = [0 0 0];
conf.secondary_sources.geometry = 'line';
conf.secondary_sources.logspread = 1.0;

%plotting
conf.plot.useplot = true;
conf.plot.usenormalisation = true;
conf.plot.normalisation = 'center';
conf.plot.mode = 'monitor';
conf.plot.usedb = false;
conf.plot.caxis = [];
conf.plot.colormap = 'default';
conf.plot.loudspeakers = false;
conf.plot.realloudspeakers = false;
conf.plot.lssize = 0.16;
conf.plot.size_unit = 'px';
conf.plot.size = [540 404];
conf.plot.cmd = '';
conf.plot.usefile = false;
conf.plot.file = '';

%% Plots
%define plane for plot
X = [-3 3];
Y = [-5 1];
Z = 0;
%virtual source
xs = [0 -1 0];
src = 'pw';
%frequency / Hz
f = 1000;

%create points for secondary source distribution
x = -2:.01:2;
y = zeros(size(x));

%WITHOUT TAPERING WINDOW
conf.usetapwin = false;
%calculate and plot sound field
sound_field_mono_wfs(X,Y,Z,xs,src,f,conf)
%add secondary source distribution
hold on
plot(x,y,'k','LineWidth',2)
%colorbar
colorbar

savefig('fig2_14a')

%WITH TAPERING WINDOW
conf.usetapwin = true;
%calculate and plot sound field
sound_field_mono_wfs(X,Y,Z,xs,src,f,conf)
%add secondary source distribution
hold on
plot(x,y,'k','LineWidth',2)
%colorbar
colorbar

savefig('fig2_14b')