%Demonstration of amplitude deviations in 2.5D WFS

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

%% Definitions for SFS Toolbox
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
conf.usetapwin = false;
conf.tapwinlen = 0.3;

%sound field simulations
conf.resolution = 600;
conf.phase = 0;

%secondary sources
conf.secondary_sources.number = 201;
conf.secondary_sources.size = 20;
conf.secondary_sources.center = [0 0 0];
conf.secondary_sources.geometry = 'line';
conf.secondary_sources.logspread = 1.0;

%plotting
conf.plot.useplot = true;
conf.plot.usenormalisation = false;
conf.plot.mode = 'monitor';
conf.plot.realloudspeakers = false;
conf.plot.lssize = 0.16;
conf.plot.size_unit = 'px';
conf.plot.size = [540 404];
conf.plot.cmd = '';
conf.plot.usefile = false;
conf.plot.file = '';

%% Plots sound fields
%define plane for plot
X = [-2 2];
Y = [-3 1];
Z = 0;
%point source/virtual source
xs = [0 1 0];
%frequency / Hz
f = 1000;

%create points for secondary source distribution
x = -2:.01:2;
y = zeros(size(x));

%POINT SOURCE FIELD
conf.plot.usedb = false;
conf.plot.caxis = [-1 1]*.05;
conf.plot.colormap = 'default';
conf.plot.loudspeakers = true;
%calculate and plot sound field
sound_field_mono_point_source(X,Y,Z,xs,f,conf)
%plot settings
set(gca,'YTick',-3:1)
ph = get(gca,'Children');
set(ph(1),'MarkerSize',3)
%change colorbar
c = colorbar;
c.Ticks = -.05:.05:.05;

savefig('fig2_03a')

%POINT SOURCE LEVEL
conf.plot.usedb = true;
conf.plot.caxis = [-42 0];
conf.plot.colormap = yellowred(14);
%calculate and plot sound field
sound_field_mono_point_source(X,Y,Z,xs,f,conf)
%plot settings
set(gca,'YTick',-3:1)
ph = get(gca,'Children');
set(ph(1),'MarkerSize',3)
%change colorbar
c = colorbar;
c.Ticks = -36:12:0;
c.TickLabels = {'-36','-24','-12','0 dB'};

savefig('fig2_03b')

%VIRTUAL SOURCE FIELD
conf.plot.usedb = false;
conf.plot.caxis = [-1 1]*.05;
conf.plot.colormap = 'default';
conf.plot.loudspeakers = false;
src = 'ps';
%calculate and plot sound field
sound_field_mono_wfs(X,Y,Z,xs,src,f,conf)
%add secondary source distribution
hold on
plot(x,y,'k','LineWidth',2)
%plot settings
set(gca,'xLim',X)
set(gca,'YTick',-3:1)
ph = get(gca,'Children');
set(ph(1),'MarkerSize',3)
%change colorbar
c = colorbar;
c.Ticks = -.05:.05:.05;

savefig('fig2_03c')

%VIRTUAL SOURCE LEVEL
conf.plot.usedb = true;
conf.plot.caxis = [-42 0];
conf.plot.colormap = yellowred(14);
%calculate and plot sound field
sound_field_mono_wfs(X,Y,Z,xs,src,f,conf)
%add secondary source distribution
hold on
plot(x,y,'k','LineWidth',2)
%plot settings
set(gca,'XLim',X)
set(gca,'YTick',-3:1)
ph = get(gca,'Children');
set(ph(1),'MarkerSize',3)
%change colorbar
c = colorbar;
c.Ticks = -36:12:0;
c.TickLabels = {'-36','-24','-12','0 dB'};

savefig('fig2_03d')

%% Plot amplitude over distance
%define examined field points
X = 0;
Y = 0:-.01:-3;
Z = 0;

%calculate amplitude over distance
conf.plot.useplot = false;
P_real = zeros(size(Y)); %preallocation
P_virtual = zeros(size(Y)); %preallocation
for n = 1:length(Y);
    P_real(n) = sound_field_mono_point_source(X,Y(n),Z,xs,f,conf);
    P_virtual(n) = sound_field_mono_wfs(X,Y(n),Z,xs,src,f,conf);
end

%plot
colours = moreland(2);
figure
    plot(Y,db(abs(P_real)),'Color',colours(1,:)), hold on
    plot(Y,db(abs(P_virtual)),'Color',colours(2,:))
    grid
    xlabel('y / m'), ylabel('amplitude / dB')
    set(gca,'XDir','reverse')
    set(gca,'YLim',[-42 0])
    set(gca,'YTick',-36:12:0)
    lh = legend('real point source','virtual point source','Location','NE');
    lh.ItemTokenSize = [10,1];

savefig('fig2_03e')