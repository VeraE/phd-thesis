%Impulse responses of real and virtual point sources in 2.5D WFS with two
%different secondary source spacings

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

%% Parameters for SFS Toolbox
conf = struct();

%misc
conf.debug = 0;

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

%binaural reproduction
conf.ir.useinterpolation = false;
conf.ir.hrirpredelay = 100;
conf.ir.usehcomp = false;

%% Calculate impulse responses
%define examined field point
X = [0 -1 0];
%head orientation
head_orientation = [pi/2 0];
%virtual source
xs = [0 1 0];
src = 'ps';
%dummy HRTFs
sofa = dummy_irs(conf);
%insert 100 zeros at beginning of dummy irs to avoid problems with resample,
%this delay can later be substracted by conf.ir.hrirpredelay
sofa.Data.IR = cat(3,zeros(1,2,conf.ir.hrirpredelay),...
    sofa.Data.IR(:,:,1:end-conf.ir.hrirpredelay));

%REAL POINT SOURCE
ir_real = ir_point_source(X,head_orientation,xs,sofa,conf);

%VIRTUAL POINT SOURCE
%secondary source distances
Delta_x0 = [.2 .5];
ir_virtual = zeros(conf.N,length(Delta_x0)); %preallocation
for n = 1:length(Delta_x0)
    conf.secondary_sources.number = conf.secondary_sources.size/Delta_x0(n)+1;
    conf.wfs.hprefhigh = conf.c/(2*Delta_x0(n));
    tmp = ir_wfs(X,head_orientation,xs,src,sofa,conf);
    ir_virtual(:,n) = tmp(:,1);
end

%remove additional delay t0 for focused sources introduced by SFS Toolbox
diameter = secondary_source_diameter(conf);
t0_delay = diameter/conf.c*conf.fs;
ir_virtual = ir_virtual(round(t0_delay)+1:end,:);
%remove pre-filter delay to time-align with real source
ir_virtual = ir_virtual(conf.wfs.hpreFIRorder/2+1:end,:);

%% Plots
colours = moreland(2);

%time vectors
t_real = (0:length(ir_real)-1)/conf.fs*1000; %/ ms
t_virtual = (0:length(ir_virtual)-1)/conf.fs*1000; %/ ms

%real point source
figure
    plot(t_real,ir_real,'Color',colours(1,:))
    grid
    axis([0 40 -.01 .05])
    xlabel('time / ms'), ylabel('amplitude')
    set(gca,'YTick',-.01:.01:.05)

savefig('fig2_10a')
    
%virtual point souce, few LS
figure
    plot(t_virtual,ir_virtual(:,2),'Color',colours(1,:))
    grid
    axis([0 40 -.01 .05])
    xlabel('time / ms'), ylabel('amplitude')
    set(gca,'YTick',-.01:.01:.05)

savefig('fig2_10b')

%virtual point source, many LS
figure
    plot(t_virtual,ir_virtual(:,1),'Color',colours(1,:))
    grid
    axis([0 40 -.01 .05])
    xlabel('time / ms'), ylabel('amplitude')
    set(gca,'YTick',-.01:.01:.05)

savefig('fig2_10c')