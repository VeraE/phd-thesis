%Magnitude responses in dB of a virtual point source in 2.5D WFS with different
%secondary source spacings

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

%% Plot
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

%secondary source distances
Delta_x0 = [.05 .1 .2 .5];

%calculate impulse responses
ir = zeros(conf.N,length(Delta_x0)); %preallocation
for n = 1:length(Delta_x0)
    conf.secondary_sources.number = conf.secondary_sources.size/Delta_x0(n)+1;
    conf.wfs.hprefhigh = conf.c/(2*Delta_x0(n));
    tmp = ir_wfs(X,head_orientation,xs,src,sofa,conf);
    ir(:,n) = tmp(:,1);
end

%remove additional delay t0 for focused sources introduced by SFS Toolbox
diameter = secondary_source_diameter(conf);
t0_delay = diameter/conf.c*conf.fs;
ir = ir(round(t0_delay)+1:end,:);

%plot
N = size(ir,1);
f = (0:N-1)/N*conf.fs;
colours = moreland(4);
figure
    for n = 1:length(Delta_x0)
        semilogx(f,db(abs(fft(ir(:,n))))+40*(n-1)+28,'Color',colours(n,:))
        hold on
        hprefhigh = conf.c/(2*Delta_x0(n));
        line([hprefhigh hprefhigh],[-40 140],'Color',colours(n,:),...
            'LineStyle','--')
    end
    %add labels for secondary source distance
    text(185,127,'\Deltax_0 = 50 cm','Color',colours(4,:),...
        'HorizontalAlignment','right')
    text(185,87,'20 cm','Color',colours(3,:),'HorizontalAlignment','right')
    text(185,47,'10 cm','Color',colours(2,:),'HorizontalAlignment','right')
    text(185,7,'5 cm','Color',colours(1,:),'HorizontalAlignment','right')
    grid
    xlabel('frequency / Hz'), ylabel('amplitude / dB')
    axis([50 20000 -40 140])
    set(gca,'XTick',[100 1000 10000])
    text(conf.c/(2*Delta_x0(4)),-55,'f_{alias,\Deltax_0=50cm}',...
        'HorizontalAlignment','center','Color',colours(4,:))
    set(gca,'XTickLabel',{'100','1k','10k'})

savefig('fig2_08')