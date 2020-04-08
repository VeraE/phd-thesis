%Demonstration of outer field when using the single layer potential

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
conf.dimension = '3D';
conf.xref = [0 0 0];
conf.usetapwin = false;
conf.tapwinlen = 0.3;

%sound field simulations
conf.resolution = 600;
conf.phase = 0;

%secondary sources
conf.secondary_sources.number = 2*50^2; %must be 2*n^2 for Gauss grid
conf.secondary_sources.size = 3;
conf.secondary_sources.center = [0 0 0];
conf.secondary_sources.geometry = 'sphere';
conf.secondary_sources.grid = 'gauss';

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
X = [-2 2];
Y = [-2 2];
Z = 0;
%frequency / Hz
f = 1000;

%create points for secondary source distribution
phi = 0:.01:2*pi;
r = conf.secondary_sources.size/2;
x = r*cos(phi);
y = r*sin(phi);

%PLANE WAVE AS VIRTUAL SOURCE
conf.driving_functions = 'default';
xs = [0 -1 0];
src = 'pw';
%calculate and plot sound field
sound_field_mono_wfs(X,Y,Z,xs,src,f,conf)
%plot settings
set(gca,'YTick',-2:2)
%add secondary source distribution
hold on
plot(x,y,'k','LineWidth',2)
%colorbar
colorbar

savefig('fig2_02a')

%POINT SOURCE AS VIRTUAL SOURCE
conf.driving_functions = 'point_source';
xs = [0 2 0]; %position
src = 'ps'; %type
%calculate and plot sound field
sound_field_mono_wfs(X,Y,Z,xs,src,f,conf)
%plot settings
set(gca,'YTick',-2:2)
%add secondary source distribution
hold on
plot(x,y,'k','LineWidth',2)
%colorbar
colorbar

savefig('fig2_02b')