winstyle = 'docked';
% winstyle = 'normal';

set(0,'DefaultFigureWindowStyle',winstyle)
set(0,'defaultaxesfontsize',18)
set(0,'defaultaxesfontname','Times New Roman')
% set(0,'defaultfigurecolor',[1 1 1])

% clear VARIABLES;
clear
close all
global spatialFactor;
global c_eps_0 c_mu_0 c_c c_eta_0
global simulationStopTimes;
global AsymForcing
global dels
global SurfHxLeft SurfHyLeft SurfEzLeft SurfHxRight SurfHyRight SurfEzRight



dels = 0.75;
spatialFactor = 1;

c_c = 299792458;                  % speed of light
c_eps_0 = 8.8542149e-12;          % vacuum permittivity
c_mu_0 = 1.2566370614e-6;         % vacuum permeability
c_eta_0 = sqrt(c_mu_0/c_eps_0);

% Variables
tSim = 200e-15 % sim time
%f = 230e12; % frequency
f = 250e12;
lambda = c_c/f; % wavelength

% Region setup
xMax{1} = 20e-6; % max x-coordinate for region
nx{1} = 200; % nx grid nodes
ny{1} = 0.75*nx{1}; % ny grid nodes

Reg.n = 1; % number of regions

% Material parameters
mu{1} = ones(nx{1},ny{1})*c_mu_0; % permeability of region

epi{1} = ones(nx{1},ny{1})*c_eps_0; % permittivity of region
% epi{1}(125:150,55:95)= c_eps_0*11.3; % sets up permittivity of inclusion within region

%3a - add more inclusions
epi{1}(125:150,45:70) = c_eps_0*11.3; % new inclusion 1
epi{1}(125:150,80:105) = c_eps_0*11.3; % new inclusion 2

sigma{1} = zeros(nx{1},ny{1}); % conductivity of region
sigmaH{1} = zeros(nx{1},ny{1}); % magnet conductivity probably?

% simulation parameters
dx = xMax{1}/nx{1}; % small spatial step
dt = 0.25*dx/c_c; % small time step
nSteps = round(tSim/dt*2); % number of steps
yMax = ny{1}*dx; % max y coordinate
nsteps_lamda = lambda/dx % number of wavelength steps

% Plotting
movie = 1; % plot as movie
Plot.off = 0; % keep plot on
Plot.pl = 0;
Plot.ori = '13'; % orientation
Plot.N = 100; % number of steps to plot
Plot.MaxEz = 1.1; % max Ez
Plot.MaxH = Plot.MaxEz/c_eta_0; % max H
Plot.pv = [0 0 90]; % plot view
Plot.reglim = [0 xMax{1} 0 yMax]; % plot limits

% Boundary Conditions
bc{1}.NumS = 2; % number of sources
bc{1}.s(1).xpos = nx{1}/(4) + 1; % source x position
bc{1}.s(1).type = 'ss'; % source type
bc{1}.s(1).fct = @PlaneWaveBC; 

% Add second source
bc{1}.s(2).xpos = 140; % source x position -- put in inclusions
bc{1}.s(2).type = 'ss'; % source type
bc{1}.s(2).fct = @PlaneWaveBC; 

% mag = -1/c_eta_0;
mag = 1; % magnitude
phi = 0; % phase
omega = f*2*pi; % radial frequency
betap = 0; % propagation
t0 = 30e-15; % initial time
%st = 15e-15;
st = -0.05;
s = 0;
y0 = yMax/2;
sty = 1.5*lambda;
bc{1}.s(1).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'}; % setup source parameters into a variable
bc{1}.s(2).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'}; % setup source parameters into a variable

Plot.y0 = round(y0/dx);

% Setup types for the boundary conditions
bc{1}.xm.type = 'a';
bc{1}.xp.type = 'a';
bc{1}.ym.type = 'a';
bc{1}.yp.type = 'a';

% absorbing boundary condition parameters
pml{1}.width = 20 * spatialFactor;
pml{1}.m = 3.5;

Reg.n  = 1;
Reg.xoff{1} = 0;
Reg.yoff{1} = 0;

RunYeeReg

%2ci - line 44 and 45 setup the epsilon and the inclusion. 
% Commenting out the inclusion created a constant epsilon that does work in simulation, but now the wave just travels.

%2cii - bc structure is used to define the number of sources, and hold the actual sources, where each source is another struct containing its properties.
% It also seems to hold the type of boundary conditions (xm, xp, ym, yp)
% Ultimately, holds the boundaries and the sources

%2ciii - sets up the source (first source, in this case)

%2civ -- xm/xp/ym/yp setup the boundary type conditions for E and H fields
% for example, in Yee2DeM, there's options for 'a', 's', 'e', 'm', 'p', where each control different sigma or Hx/Hy
% a = absorbing, s = hard source, e = electric field on wall, m = magnetic wall, p = periodic
% xp.type='e' does: Hy{k}(end,:) = Hy{k}(end-1,:); which seems to reflect the wave once it hits the second region/inclusion
% from lecture, 'e' puts an electric field on the wall

%3b - st seems to vary the time step or something similar. Before, it was smooth, after, the simulation is choppy
% Even though simulation was choppy, it was easier to see the effect of the inclusions





