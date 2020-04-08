%Illustration of test-retest reliability in a scatter plot

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

%% Create datapoints
rng(3); %set seed for random number generator
N = 20;
x = rand(N,1)*10;
y = x+(rand(N,1)-.5)*3;

%% Calculate Pearson's correlation coefficient
rho = corr(x,y);

%% Plot
colours = moreland(2);
figure
    scatter(x,y,10,colours(1,:),'filled'), hold on
    line([0 10],[0 10],'Color','k')
    grid
    xlabel('test score'), ylabel('retest score')
    axis equal
    axis([0 10 0 10])
    set(gca,'XTick',0:10)

savefig('fig3_01')