%Sound pressure of a virtual point source in 2.D WFS at different frequencies
%(demonstrating spatial aliasing)

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
conf.usetapwin = true;
conf.tapwinlen = 0.3;

%sound field simulations
conf.resolution = 600;
conf.phase = 0;

%secondary sources
conf.secondary_sources.number = 16;
conf.secondary_sources.size = 3;
conf.secondary_sources.center = [0 0 0];
conf.secondary_sources.geometry = 'line';
conf.secondary_sources.logspread = 1.0;

%plotting
conf.plot.useplot = false;
conf.plot.usenormalisation = false;
conf.plot.mode = 'monitor';
conf.plot.usedb = false;
conf.plot.caxis = [];
conf.plot.colormap = 'default';
conf.plot.loudspeakers = true;
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
Y = [-4 0];
Z = 0;
%virtual source
xs = [0 1 0];
src = 'ps';

%NO ALIASING
%frequency / Hz
f = 600;
%calculate and plot sound field
[P,x,y] = sound_field_mono_wfs(X,Y,Z,xs,src,f,conf);
%normalise to reference position and plot
[~,xnorm] = min(abs(x(1,:)-conf.xref(1)));
[~,ynorm] = min(abs(y(:,1)-conf.xref(2)));
P = P/max(abs(P(xnorm(1),ynorm(1))));
plot_sound_field(P,X,Y,Z,secondary_source_positions(conf),conf)
%plot settings
set(gca,'YTick',-4:0)
ph = get(gca,'Children');
set(ph(1),'MarkerSize',3)
%colorbar
colorbar

savefig('fig2_07a')

%SOME ALIASING
%frequency / Hz
f = 1200;
%calculate and plot sound field
[P,x,y] = sound_field_mono_wfs(X,Y,Z,xs,src,f,conf);
%normalise to reference position and plot
[~,xnorm] = min(abs(x(1,:)-conf.xref(1)));
[~,ynorm] = min(abs(y(:,1)-conf.xref(2)));
P = P/max(abs(P(xnorm(1),ynorm(1))));
plot_sound_field(P,X,Y,Z,secondary_source_positions(conf),conf)
%plot settings
set(gca,'YTick',-4:0)
ph = get(gca,'Children');
set(ph(1),'MarkerSize',3)
%colorbar
colorbar

savefig('fig2_07b')

%HEAVY ALIASING
%frequency / Hz
f = 2400;
%calculate and plot sound field
[P,x,y] = sound_field_mono_wfs(X,Y,Z,xs,src,f,conf);
%normalise to reference position and plot
[~,xnorm] = min(abs(x(1,:)-conf.xref(1)));
[~,ynorm] = min(abs(y(:,1)-conf.xref(2)));
P = P/max(abs(P(xnorm(1),ynorm(1))));
plot_sound_field(P,X,Y,Z,secondary_source_positions(conf),conf)
%plot settings
set(gca,'YTick',-4:0)
ph = get(gca,'Children');
set(ph(1),'MarkerSize',3)
%colorbar
colorbar

savefig('fig2_07c')