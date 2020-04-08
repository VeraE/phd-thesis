%Sound field components separated in propagating and evanescent parts for virtual
%point sources in 2.5D WFS close to the secondary source array

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
x0 = -20:dx0:+20; y0 = x0*0; z0 = y0*0;
x0 = [x0;y0;z0]; clear y0; clear z0;
N_Monopole = size(x0,2);
xref = [0; 1; 0]; %reference point / m
yref = xref(2);

%% Calculate field of a real point source
rps = abs(sqrt((X-xps_behind(1)).^2+(Y-xps_behind(2)).^2));
Pps = pA * exp(-1i*w/c*rps) ./ rps;

%% Aliasing components
%range of wave numbers
kx = linspace(-50,50,20001);
dkx = kx(2)-kx(1);
[KX,Y] = meshgrid(kx,abs(y));

%evaluated range of eta
maxeta = 10;
eta = -maxeta:maxeta;

%BEHIND
%1. sound field component P_S,pr1
P_Spr1_behind = zeros(length(y),length(kx));
for n = 1:length(eta)
    idx = find(abs(kx-eta(n)*2*pi/dx0) < abs(w/c) & abs(kx) < (w/c));
    P_Spr1_behind(:,idx) = P_Spr1_behind(:,idx)...
        -1i*pi*sqrt(yref/(yref-xps_behind(2)))...
        *(  besselj(0,sqrt((w/c)^2-KX(:,idx).^2).*sqrt(Y(:,idx).^2))...
        -1i*bessely(0,sqrt((w/c)^2-KX(:,idx).^2).*sqrt(Y(:,idx).^2)))...
        .*exp(1i*(KX(:,idx)-2*pi/dx0*eta(n))*xps_behind(1))...
        .*exp(1i*sqrt((w/c)^2-(KX(:,idx)-2*pi/dx0*eta(n)).^2)*xps_behind(2));
end
%2. sound field component P_S,pr2
P_Spr2_behind = zeros(length(y),length(kx));
for n = 1:length(eta)
    idx = find(abs(kx-eta(n)*2*pi/dx0) > abs(w/c) & abs(kx) < (w/c));
    P_Spr2_behind(:,idx) = P_Spr2_behind(:,idx)...
        -1i*pi*sqrt(yref/(yref-xps_behind(2)))...
        *(  besselj(0,sqrt((w/c)^2-KX(:,idx).^2).*sqrt(Y(:,idx).^2))...
        -1i*bessely(0,sqrt((w/c)^2-KX(:,idx).^2).*sqrt(Y(:,idx).^2)))...
        .*exp(1i*(KX(:,idx)-2*pi/dx0*eta(n))*xps_behind(1))...
        .*exp(sqrt((KX(:,idx)-2*pi/dx0*eta(n)).^2-(w/c)^2)*xps_behind(2));
end
%3. sound field component P_S,ev1
P_Sev1_behind = zeros(length(y),length(kx));
for n = 1:length(eta)
    idx = find(abs(kx-eta(n)*2*pi/dx0) < abs(w/c) & abs(kx) > (w/c));
    P_Sev1_behind(:,idx) = P_Sev1_behind(:,idx)...
        +2*sqrt(yref/(yref-xps_behind(2)))...
        *besselk(0,sqrt(KX(:,idx).^2-(w/c)^2).*sqrt(Y(:,idx).^2))...
        .*exp(1i*(KX(:,idx)-2*pi/dx0*eta(n))*xps_behind(1))...
        .*exp(1i*sqrt((w/c)^2-(KX(:,idx)-2*pi/dx0*eta(n)).^2)*xps_behind(2));
end
%4. sound field component P_S,ev2
P_Sev2_behind = zeros(length(y),length(kx));
for n = 1:length(eta)
    idx = find(abs(kx-eta(n)*2*pi/dx0) > abs(w/c) & abs(kx) > (w/c));
    P_Sev2_behind(:,idx) = P_Sev2_behind(:,idx)...
        +2*sqrt(yref/(yref-xps_behind(2)))...
        *besselk(0,sqrt(KX(:,idx).^2-(w/c)^2).*sqrt(Y(:,idx).^2))...
        .*exp(1i*(KX(:,idx)-2*pi/dx0*eta(n))*xps_behind(1))...
        .*exp(sqrt((KX(:,idx)-2*pi/dx0*eta(n)).^2-(w/c)^2)*xps_behind(2));
end

%BETWEEN
%1. sound field component P_S,pr1
P_Spr1_between = zeros(length(y),length(kx));
for n = 1:length(eta)
    idx = find(abs(kx-eta(n)*2*pi/dx0) < abs(w/c) & abs(kx) < (w/c));
    P_Spr1_between(:,idx) = P_Spr1_between(:,idx)...
        -1i*pi*sqrt(yref/(yref-xps_between(2)))...
        *(  besselj(0,sqrt((w/c)^2-KX(:,idx).^2).*sqrt(Y(:,idx).^2))...
        -1i*bessely(0,sqrt((w/c)^2-KX(:,idx).^2).*sqrt(Y(:,idx).^2)))...
        .*exp(1i*(KX(:,idx)-2*pi/dx0*eta(n))*xps_between(1))...
        .*exp(1i*sqrt((w/c)^2-(KX(:,idx)-2*pi/dx0*eta(n)).^2)*xps_between(2));
end
%2. sound field component P_S,pr2
P_Spr2_between = zeros(length(y),length(kx));
for n = 1:length(eta)
    idx = find(abs(kx-eta(n)*2*pi/dx0) > abs(w/c) & abs(kx) < (w/c));
    P_Spr2_between(:,idx) = P_Spr2_between(:,idx)...
        -1i*pi*sqrt(yref/(yref-xps_between(2)))...
        *(  besselj(0,sqrt((w/c)^2-KX(:,idx).^2).*sqrt(Y(:,idx).^2))...
        -1i*bessely(0,sqrt((w/c)^2-KX(:,idx).^2).*sqrt(Y(:,idx).^2)))...
        .*exp(1i*(KX(:,idx)-2*pi/dx0*eta(n))*xps_between(1))...
        .*exp(sqrt((KX(:,idx)-2*pi/dx0*eta(n)).^2-(w/c)^2)*xps_between(2));
end
%3. sound field component P_S,ev1
P_Sev1_between = zeros(length(y),length(kx));
for n = 1:length(eta)
    idx = find(abs(kx-eta(n)*2*pi/dx0) < abs(w/c) & abs(kx) > (w/c));
    P_Sev1_between(:,idx) = P_Sev1_between(:,idx)...
        +2*sqrt(yref/(yref-xps_between(2)))...
        *besselk(0,sqrt(KX(:,idx).^2-(w/c)^2).*sqrt(Y(:,idx).^2))...
        .*exp(1i*(KX(:,idx)-2*pi/dx0*eta(n))*xps_between(1))...
        .*exp(1i*sqrt((w/c)^2-(KX(:,idx)-2*pi/dx0*eta(n)).^2)*xps_between(2));
end
%4. sound field component P_S,ev2
P_Sev2_between = zeros(length(y),length(kx));
for n = 1:length(eta)
    idx = find(abs(kx-eta(n)*2*pi/dx0) > abs(w/c) & abs(kx) > (w/c));
    P_Sev2_between(:,idx) = P_Sev2_between(:,idx)...
        +2*sqrt(yref/(yref-xps_between(2)))...
        *besselk(0,sqrt(KX(:,idx).^2-(w/c)^2).*sqrt(Y(:,idx).^2))...
        .*exp(1i*(KX(:,idx)-2*pi/dx0*eta(n))*xps_between(1))...
        .*exp(sqrt((KX(:,idx)-2*pi/dx0*eta(n)).^2-(w/c)^2)*xps_between(2));
end

%% Inverse spatial Fourier transform
%BEHIND
P1_behind = zeros(length(y),length(x)); %preallocation
P2_behind = zeros(length(y),length(x)); %preallocation
P3_behind = zeros(length(y),length(x)); %preallocation
P4_behind = zeros(length(y),length(x)); %preallocation
for n = 1:length(x)
        P1_behind(:,n) = 1/(2*pi)*sum(P_Spr1_behind.*...
            repmat(exp(-1i*kx*x(n)),length(y),1),2)*dkx;
        P2_behind(:,n) = 1/(2*pi)*sum(P_Spr2_behind.*...
            repmat(exp(-1i*kx*x(n)),length(y),1),2)*dkx;
        P3_behind(:,n) = 1/(2*pi)*sum(P_Sev1_behind.*...
            repmat(exp(-1i*kx*x(n)),length(y),1),2)*dkx;
        P4_behind(:,n) = 1/(2*pi)*sum(P_Sev2_behind.*...
            repmat(exp(-1i*kx*x(n)),length(y),1),2)*dkx;
end
%take signal amplitude into account
P1_behind = pA*P1_behind;
P2_behind = pA*P2_behind;
P3_behind = pA*P3_behind;
P4_behind = pA*P4_behind;
%addition of sound field components
Pges_behind = P1_behind+P2_behind+P3_behind+P4_behind;

%BETWEEN
P1_between = zeros(length(y),length(x)); %preallocation
P2_between = zeros(length(y),length(x)); %preallocation
P3_between = zeros(length(y),length(x)); %preallocation
P4_between = zeros(length(y),length(x)); %preallocation
for n = 1:length(x)
        P1_between(:,n) = 1/(2*pi)*sum(P_Spr1_between.*...
            repmat(exp(-1i*kx*x(n)),length(y),1),2)*dkx;
        P2_between(:,n) = 1/(2*pi)*sum(P_Spr2_between.*...
            repmat(exp(-1i*kx*x(n)),length(y),1),2)*dkx;
        P3_between(:,n) = 1/(2*pi)*sum(P_Sev1_between.*...
            repmat(exp(-1i*kx*x(n)),length(y),1),2)*dkx;
        P4_between(:,n) = 1/(2*pi)*sum(P_Sev2_between.*...
            repmat(exp(-1i*kx*x(n)),length(y),1),2)*dkx;
end
%take signal amplitude into account
P1_between = pA*P1_between;
P2_between = pA*P2_between;
P3_between = pA*P3_between;
P4_between = pA*P4_between;
%addition of sound field components
Pges_between = P1_between+P2_between+P3_between+P4_between;

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
Pges_ref_behind   = Pges_behind(indy,indx);
Pges_refdB_behind = 20*log10(abs(Pges_ref_behind/sqrt(2))/p0)
Pges_phase_behind = angle(Pges_ref_behind);
%BETWEEN
Pges_ref_between   = Pges_between(indy,indx);
Pges_refdB_between = 20*log10(abs(Pges_ref_between/sqrt(2))/p0)
Pges_phase_between = angle(Pges_ref_between);

%% Plots

%BEHIND
%P1
figure  
    imagesc(x,y,real(P1_behind)), hold on
    %secondary sources
    plot3(x0(1,1:N_Monopole),x0(2,1:N_Monopole),x0(3,1:N_Monopole)+1e6,'ok',...
        'MarkerSize',3,'MarkerFaceColor','k')
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

savefig('fig2_12a')

%P2
figure  
    imagesc(x,y,real(P2_behind)), hold on
    %secondary sources
    plot3(x0(1,1:N_Monopole),x0(2,1:N_Monopole),x0(3,1:N_Monopole)+1e6,'ok',...
        'MarkerSize',3,'MarkerFaceColor','k')
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

savefig('fig2_12b')

%P3
figure  
    imagesc(x,y,real(P3_behind)), hold on
    %secondary sources
    plot3(x0(1,1:N_Monopole),x0(2,1:N_Monopole),x0(3,1:N_Monopole)+1e6,'ok',...
        'MarkerSize',3,'MarkerFaceColor','k')
    %colorbar
    colormap moreland(256)
    set(gca,'CLim',[-.5 .5])
    c = colorbar;
    c.Ticks = -.5:.5:.5;
    c.TickLabels = {'-0.5','0','0.5 Pa'};
    %axes properties
    axis equal
    axis([min(x) max(x) min(y) max(y)])
    set(gca,'XTick',[min(x):max(x)])
    set(gca,'YTick',[min(y):max(y)])
    set(gca,'YTickLabel',{'0','-1','-2','-3','-4'})
    xlabel('x / m'), ylabel('y / m')

savefig('fig2_12c')

%P4
figure  
    imagesc(x,y,real(P4_behind)), hold on
    %secondary sources
    plot3(x0(1,1:N_Monopole),x0(2,1:N_Monopole),x0(3,1:N_Monopole)+1e6,'ok',...
        'MarkerSize',3,'MarkerFaceColor','k')
    %colorbar
    colormap moreland(256)
    set(gca,'CLim',[-.5 .5])
    c = colorbar;
    c.Ticks = -.5:.5:.5;
    c.TickLabels = {'-0.5','0','0.5 Pa'};
    %axes properties
    axis equal
    axis([min(x) max(x) min(y) max(y)])
    set(gca,'XTick',[min(x):max(x)])
    set(gca,'YTick',[min(y):max(y)])
    set(gca,'YTickLabel',{'0','-1','-2','-3','-4'})
    xlabel('x / m'), ylabel('y / m')

savefig('fig2_12d')

%BETWEEN
%P1
figure  
    imagesc(x,y,real(P1_between)), hold on
    %secondary sources
    plot3(x0(1,1:N_Monopole),x0(2,1:N_Monopole),x0(3,1:N_Monopole)+1e6,'ok',...
        'MarkerSize',3,'MarkerFaceColor','k')
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

savefig('fig2_13a')

%P2
figure  
    imagesc(x,y,real(P2_between)), hold on
    %secondary sources
    plot3(x0(1,1:N_Monopole),x0(2,1:N_Monopole),x0(3,1:N_Monopole)+1e6,'ok',...
        'MarkerSize',3,'MarkerFaceColor','k')
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

savefig('fig2_13b')

%P3
figure  
    imagesc(x,y,real(P3_between)), hold on
    %secondary sources
    plot3(x0(1,1:N_Monopole),x0(2,1:N_Monopole),x0(3,1:N_Monopole)+1e6,'ok',...
        'MarkerSize',3,'MarkerFaceColor','k')
    %colorbar
    colormap moreland(256)
    set(gca,'CLim',[-.5 .5])
    c = colorbar;
    c.Ticks = -.5:.5:.5;
    c.TickLabels = {'-0.5','0','0.5 Pa'};
    %axes properties
    axis equal
    axis([min(x) max(x) min(y) max(y)])
    set(gca,'XTick',[min(x):max(x)])
    set(gca,'YTick',[min(y):max(y)])
    set(gca,'YTickLabel',{'0','-1','-2','-3','-4'})
    xlabel('x / m'), ylabel('y / m')

savefig('fig2_13c')

%P4
figure  
    imagesc(x,y,real(P4_between)), hold on
    %secondary sources
    plot3(x0(1,1:N_Monopole),x0(2,1:N_Monopole),x0(3,1:N_Monopole)+1e6,'ok',...
        'MarkerSize',3,'MarkerFaceColor','k')
    %colorbar
    colormap moreland(256)
    set(gca,'CLim',[-.5 .5])
    c = colorbar;
    c.Ticks = -.5:.5:.5;
    c.TickLabels = {'-0.5','0','0.5 Pa'};
    %axes properties
    axis equal
    axis([min(x) max(x) min(y) max(y)])
    set(gca,'XTick',[min(x):max(x)])
    set(gca,'YTick',[min(y):max(y)])
    set(gca,'YTickLabel',{'0','-1','-2','-3','-4'})
    xlabel('x / m'), ylabel('y / m')

savefig('fig2_13d')