%Sound field of a virtual point source in 2.5D WFS emitting a broadband
%impulse in dB

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
conf.debug = 0;
conf.showprogress = false;

%audio
conf.fs = 44100;
conf.c = 343;

%delayline
conf.delayline.filter = 'integer';
conf.delayline.resampling = 'matlab';
conf.delayline.resamplingfactor = 100;

%sound field synthesis
conf.dimension = '2.5D';
conf.driving_functions = 'default';
conf.N = 4096;
conf.xref = [0 -1 0];
conf.usetapwin = true;
conf.tapwinlen = 0.3;
conf.t0 = 'source';
conf.usebandpass = true;
conf.bandpassflow = 10;
conf.bandpassfhigh = 20000;

%sound field simulations
conf.resolution = 600;

%secondary sources
conf.secondary_sources.size = 10;
conf.secondary_sources.center = [0 0 0];
conf.secondary_sources.geometry = 'line';
conf.secondary_sources.logspread = 1.0;

%Wave Field Synthesis
conf.wfs.usehpre = true;
conf.wfs.hpretype = 'FIR';
conf.wfs.hpreflow = 50;
conf.wfs.hpreFIRorder = 512;

%plotting
conf.plot.useplot = true;
conf.plot.usenormalisation = true;
conf.plot.normalisation = 'max';
conf.plot.mode = 'monitor';
conf.plot.usedb = true;
conf.plot.caxis = [-45 0];
conf.plot.colormap = 'default';
conf.plot.loudspeakers = true;
conf.plot.realloudspeakers = false;
conf.plot.lssize = 0.16;
conf.plot.size_unit = 'px';
conf.plot.size = [540 404];
conf.plot.cmd = '';
conf.plot.usefile = false;
conf.plot.file = '';

%% Plot
%define plane for plot
X = [-2 2];
Y = [-4 0];
Z = 0;
%virtual source
xs = [0 1 0];
src = 'ps';
%time
t = 3/conf.c;

%NO ALIASING
Delta_x0 = .005;
conf.secondary_sources.number = conf.secondary_sources.size/Delta_x0+1;
conf.wfs.hprefhigh = conf.c/(2*Delta_x0);
%calculate and plot sound field
sound_field_imp_wfs(X,Y,Z,xs,src,t,conf)
%plot settings
set(gca,'XLim',X)
set(gca,'YTick',-4:0)
ph = get(gca,'Children');
set(ph(1),'MarkerSize',3)
%change colorbar
c = colorbar;
c.Ticks = -45:15:0;
c.TickLabels = {'-45','-30','-15','0 dB'};

savefig('fig2_09a')

%SOME ALIASING
Delta_x0 = .2;
conf.secondary_sources.number = conf.secondary_sources.size/Delta_x0+1;
conf.wfs.hprefhigh = conf.c/(2*Delta_x0);
%calculate and plot sound field
sound_field_imp_wfs(X,Y,Z,xs,src,t,conf)
%plot settings
set(gca,'XLim',X)
set(gca,'YTick',-4:0)
ph = get(gca,'Children');
set(ph(1),'MarkerSize',3)
%change colorbar
c = colorbar;
c.Ticks = -45:15:0
c.TickLabels = {'-45','-30','-15','0 dB'};

savefig('fig2_09b')

%HEAVY ALIASING
Delta_x0 = .5;
conf.secondary_sources.number = conf.secondary_sources.size/Delta_x0+1;
conf.wfs.hprefhigh = conf.c/(2*Delta_x0);
%calculate and plot sound field
sound_field_imp_wfs(X,Y,Z,xs,src,t,conf)
%plot settings
set(gca,'XLim',X)
set(gca,'YTick',-4:0)
ph = get(gca,'Children');
set(ph(1),'MarkerSize',3)
%change colorbar
c = colorbar;
c.Ticks = -60:20:0;
c.TickLabels = {'-60','-40','-20','0 dB'};

savefig('fig2_09c')