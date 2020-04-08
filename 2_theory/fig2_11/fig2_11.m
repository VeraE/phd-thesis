%Sound fields of a real and two virtual point sources in 2.5D WFS, that are
%located close to the secondary source array

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
f = 350; %frequency / Hz
c = 343; %speed of sound / m/s
w = 2*pi*f;
lambda = c/f;

pA = sqrt(2); %amplitude / Pa
p0 = 2*10^-5; %reference sound pressure / Pa

%% Positions of virtual point sources
xps_behind = [0 -.02 0];
xps_between = [.1 -.02 0];

%% xy-plane
res = .01; %/ m
x   = -2:res:2; y = 0:res:4; %/ m
[X,Y] = meshgrid(x,y);

%% Secondary sources
dx0 = .2; %secondary source spacing / m
x0 = -50:dx0:+50; y0 = x0*0; z0 = y0*0;
x0 = [x0;y0;z0]; clear y0; clear z0;
N_Monopole = size(x0,2);
xref = [0; 1; 0]; %reference point / m
yref = xref(2);

%% Calculate field of a real point source
rps = abs(sqrt((X-xps_behind(1)).^2+(Y-xps_behind(2)).^2));
Pps = pA * exp(-1i*w/c*rps) ./ rps;

%% Virtual point source with SDM
%distances r: points in xy-plane to secondary sources
r = zeros(N_Monopole,length(x)*length(y)); %preallocation
for nn = 1:N_Monopole
    r_tmp   = sqrt( (X-x0(1,nn)).^2 + (Y-x0(2,nn)).^2 );
    r_tmp   = r_tmp';
    r(nn,:) = r_tmp(:);
end
clear r_tmp

%Green's function
G_x_x0 = exp(-1i*w/c*r)./(4*pi*r);

D1appr_behind = zeros(N_Monopole,1); %preallocation
%pre-equalisation filter
F1appr_behind =  4*pi*1/2*sqrt(yref/(yref-xps_behind(2)))*1i*w/c*xps_behind(2);
F1appr_between =  4*pi*1/2*sqrt(yref/(yref-xps_between(2)))*1i*w/c*xps_between(2);
%driving functions
for nnn = 1:N_Monopole
    %behind
    r_x02xps_behind = sqrt((x0(1,nnn)-xps_behind(1))^2 ...
        +(x0(2,nnn)-xps_behind(2))^2+(x0(3,nnn)-xps_behind(3))^2);
    D1appr_behind(nnn,:) = pA * F1appr_behind * ...
        (besselj(1,w/c*r_x02xps_behind) - 1i*bessely(1,w/c*r_x02xps_behind))...
        / r_x02xps_behind; 
    %between
    r_x02xps_between = sqrt((x0(1,nnn)-xps_between(1))^2 ...
        +(x0(2,nnn)-xps_between(2))^2+(x0(3,nnn)-xps_between(3))^2);
    D1appr_between(nnn,:) = pA * F1appr_between * ...
        (besselj(1,w/c*r_x02xps_between) - 1i*bessely(1,w/c*r_x02xps_between))...
        / r_x02xps_between; 
end

%calculate sound fields
%behind
D1appr_behind = D1appr_behind.'
P1appr_behind = zeros(1,length(x)*length(y)); %preallocation
P1appr_behind = P1appr_behind + D1appr_behind*G_x_x0;    
P1appr_behind = reshape(P1appr_behind,length(y),length(x)).';
P1appr_behind = P1appr_behind*dx0;
%between
D1appr_between = D1appr_between.';
P1appr_between = zeros(1,length(x)*length(y)); %preallocation
P1appr_between = P1appr_between + D1appr_between*G_x_x0;    
P1appr_between = reshape(P1appr_between,length(y),length(x)).';
P1appr_between = P1appr_between*dx0;

%% Calculate level and phase at reference point
%find index for reference line
indy = find(y>=yref,1,'first');
%find index for y-axis
indx = find(x>=0,1,'first');
%sound pressure, level and phase of real point source at (0,yref,0)
Pps_ref   = Pps(indy,indx);
Pps_refdB = 20*log10(abs(Pps_ref/sqrt(2))/p0)
Pps_phase = angle(Pps_ref);
%sound pressure, level and phase of virtual point sources at (0,yref,0)
%BEHIND
P1appr_ref_behind   = P1appr_behind(indy,indx);
P1appr_refdB_behind = 20*log10(abs(P1appr_ref_behind/sqrt(2))/p0)
P1appr_phase_behind = angle(P1appr_ref_behind);
%BETWEEN
P1appr_ref_between   = P1appr_between(indy,indx);
P1appr_refdB_between = 20*log10(abs(P1appr_ref_between/sqrt(2))/p0)
P1appr_phase_between = angle(P1appr_ref_between);

%% Plots
%FIELDS

%real source
figure  
    imagesc(x,y,real(Pps)), hold on
    %colorbar
    colormap moreland(256)
    set(gca,'CLim',[-2 2])
    c = colorbar;
    c.Ticks = -2:2;
    c.TickLabels = {'-2','-1','0','1','2 Pa'};
    %axes properties
    axis equal
    axis([min(x) max(x) min(y) max(y)])
    set(gca,'XTick',[min(x):max(x)])
    set(gca,'YTick',[min(y):max(y)])
    set(gca,'YTickLabel',{'0','-1','-2','-3','-4'})
    xlabel('x / m'), ylabel('y / m')

savefig('fig2_11a')

%virtual source, behind
figure  
    imagesc(x,y,real(P1appr_behind)), hold on
    %secondary sources
    plot3(x0(1,1:N_Monopole),x0(2,1:N_Monopole),x0(3,1:N_Monopole)+1e6,'ok',...
        'MarkerSize',3,'MarkerFaceColor','k')
    %reference line
    line([min(x) max(x)],[yref yref],'Color','k','LineStyle','--')
    %colorbar
    colormap moreland(256)
    set(gca,'CLim',[-2 2])
    c = colorbar;
    c.Ticks = -2:2;
    c.TickLabels = {'-2','-1','0','1','2 Pa'};
    %axes properties
    axis equal
    axis([min(x) max(x) min(y) max(y)])
    set(gca,'XTick',[min(x):max(x)])
    set(gca,'YTick',[min(y):max(y)])
    set(gca,'YTickLabel',{'0','-1','-2','-3','-4'})
    xlabel('x / m'), ylabel('y / m')

savefig('fig2_11b')

%virtual source, between
figure  
    imagesc(x,y,real(P1appr_between)), hold on
    %secondary sources
    plot3(x0(1,1:N_Monopole),x0(2,1:N_Monopole),x0(3,1:N_Monopole)+1e6,'ok',...
        'MarkerSize',3,'MarkerFaceColor','k')
    %reference line
    line([min(x) max(x)],[yref yref],'Color','k','LineStyle','--')
    %colorbar
    colormap moreland(256)
    set(gca,'CLim',[-2 2])
    c = colorbar;
    c.Ticks = -2:2;
    c.TickLabels = {'-2','-1','0','1','2 Pa'};
    %axes properties
    axis equal
    axis([min(x) max(x) min(y) max(y)])
    set(gca,'XTick',[min(x):max(x)])
    set(gca,'YTick',[min(y):max(y)])
    set(gca,'YTickLabel',{'0','-1','-2','-3','-4'})
    xlabel('x / m'), ylabel('y / m')

savefig('fig2_11c')

%LEVELS

%real source
figure  
    imagesc(x,y,20*log10(abs(Pps/sqrt(2))/p0)), hold on
    %colorbar
    colormap yellowred(20)
    set(gca,'CLim',[64 124])
    c = colorbar;
    c.Ticks = 64:15:124;
    c.TickLabels = {'64','79','94','109','124 dB_{SPL}'};
    %axes properties
    axis equal
    axis([min(x) max(x) min(y) max(y)])
    set(gca,'XTick',[min(x):max(x)])
    set(gca,'YTick',[min(y):max(y)])
    set(gca,'YTickLabel',{'0','-1','-2','-3','-4'})
    xlabel('x / m'), ylabel('y / m')

savefig('fig2_11d')

%virtual source, behind
figure  
    imagesc(x,y,20*log10(abs(P1appr_behind/sqrt(2))/p0)), hold on
    %secondary sources
    plot3(x0(1,1:N_Monopole),x0(2,1:N_Monopole),x0(3,1:N_Monopole)+1e6,'ok',...
        'MarkerSize',3,'MarkerFaceColor','k')
    %reference line
    line([min(x) max(x)],[yref yref],'Color','k','LineStyle','--')
    %colorbar
    colormap yellowred(20)
    set(gca,'CLim',[64 124])
    c = colorbar;
    c.Ticks = 64:15:124;
    c.TickLabels = {'64','79','94','109','124 dB_{SPL}'};
    %axes properties
    axis equal
    axis([min(x) max(x) min(y) max(y)])
    set(gca,'XTick',[min(x):max(x)])
    set(gca,'YTick',[min(y):max(y)])
    set(gca,'YTickLabel',{'0','-1','-2','-3','-4'})
    xlabel('x / m'), ylabel('y / m')

savefig('fig2_11e')

%virtual source, between
figure  
    imagesc(x,y,20*log10(abs(P1appr_between/sqrt(2))/p0)), hold on
    %secondary sources
    plot3(x0(1,1:N_Monopole),x0(2,1:N_Monopole),x0(3,1:N_Monopole)+1e6,'ok',...
        'MarkerSize',3,'MarkerFaceColor','k')
    %reference line
    line([min(x) max(x)],[yref yref],'Color','k','LineStyle','--')
    %colorbar
    colormap yellowred(20)
    set(gca,'CLim',[64 124])
    c = colorbar;
    c.Ticks = 64:15:124;
    c.TickLabels = {'64','79','94','109','124 dB_{SPL}'};
    %axes properties
    axis equal
    axis([min(x) max(x) min(y) max(y)])
    set(gca,'XTick',[min(x):max(x)])
    set(gca,'YTick',[min(y):max(y)])
    set(gca,'YTickLabel',{'0','-1','-2','-3','-4'})
    xlabel('x / m'), ylabel('y / m')

savefig('fig2_11f')