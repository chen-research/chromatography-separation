% function parameter_to_data(nCompIn,parameterIn)
function [t, C, Tsol] = parameter_to_solutions(nCompIn,parameterIn,Injection)

% nCompIn: Number of injected components. Suggestion: nCompIn>=2.
% parameterIn: input matrix of size nCompIn*4. 
% Type: [aI_i, bI_i, aII_i, bII_i], where 
% i=1, ..., nCompIn
% For nCompIn=2: 
% parameterIn=[2, 0.1, 1, 0.05; 4, 0.2, 2, 0.1].
% For nCompIn=3: 
% parameterIn=[2, 0.1, 1, 0.05; 4, 0.2, 2, 0.1; 6, 0.3, 3, 0.15].

% Ye Zhang & Patrik Forssen
% June 2016
% Karlstad University


% PROBLEM DEFINITION
% Problem is defined in a structure "sysOpt"

% Column characteristics
sysOpt.L       = 15;              % Column length [cm]
sysOpt.D       =  0.46;              % Column diameter [cm]
sysOpt.Epsilon =  0.5616;            % Porosity
% Note that often porosity Epsilon is used instead of phase ratio F, here we
% have that F = (1 - Epsilon)/Epsilon

% Linear velocity
sysOpt.Flow    = 0.78;               % Flow [mL/min]
% Note that usually the flow rate is given instead of the linear velocity u, u
% is then calculated from the flow rate and the column characteristics

% Injection (boundary conditions)
sysOpt.Vinj    = 50;              % Injection volume [uL]
sysOpt.Csamp   = Injection;             % Sample concentration [mM]
% sysOpt.Csamp can be a vector with the individual concentrations for the
% components or a scalar when they are equal. Note that here we assume a sqaure
% injection profile.

% Adsorption isotherm
nComp          = nCompIn;                % Number of injected components
% sysOpt.IsoPar  = [...
%   (2:2:nComp*2)', (0.1 :0.1 :0.1 *nComp)', ...
%   (1:1:nComp*1)', (0.05:0.05:0.05*nComp)'];  % [aI, bI, aII, bII]

sysOpt.IsoPar  = parameterIn;
sysOpt.CompAds = 1;                % Set to 0 for non-competativ adsorption
% Here we use an n-site Langmuir adsorption isotherm where the adsorption
% parameters are given in a matrix. In each row the adsorption parameters for an
% injected component is given as a sequence of the a and b parameters for each
% adsorption site.

% Diffusion
sysOpt.Nx      = 9000;             % "Number of theoretical plates"
% sysOpt.Nx can be a vector with the individual number of theoretical plates for
% the components or a scalar when they are equal. Note that usually the "number
% of theoretical plates" is given instead of the diffusion constants Da, here we
% have that Da = L*u/(2*Nx).

% Kinetic constants
%sysOpt.kf      = [...
%  (10000:10000:10000*nComp)', (10:10:10*nComp)'];
sysOpt.kf      = (1:1:1*nComp)';
% sysOpt.kf can be a scalar that is used for all components and adsorption
% sites, it can be a vector with individual kinetic constans for each adsorption
% site that is used for all components or it can be a matrix with the idividual
% kinetic constans for each componnet and each adsorption site.
% NOTE THAT sysOpt.kf SHOULD BE EMPTY FOR THE EQUILIBRIUM-DISPERSIVE MODEL!


% SOLVER OPTIONS
% In the structure "solOpt" (optional) the user can change the default options
% for the solver
solOpt.xPoints     = [];  % Space discretization points (default 200)
solOpt.Tsol        = [];  % Solver time range (default calculated)
% ODE solver options (see "odeset" for more information)
solOpt.NormControl = '';  % 'on' or 'off' (default 'on')
solOpt.RelTol      = [];  % Default 1e-4 for ED model and 1e-3 for TD model
solOpt.AbsTol      = [];  % Default 1e-6


% CALL SOLVER
%sysOpt.kf = [];  % Use ED model (comment this line to use TD model!)
[t, C, Tsol]    = finiteVolumes(sysOpt, solOpt);

% Plot
% figure
% hold on
% box on
% if (isempty(sysOpt.kf))
%   title('EQUILIBRIUM-DISPERSIVE MODEL')
% else
%   title('TRANSPORT-DISPERSIVE MODEL')
% end
% xlabel('\itt\rm [min.]')
% ylabel('\itC\rm [M]')
% plot(t, C, 'linewidth', 1.5)
% yl = ylim;
% ylim([0 yl(2)])

% save dataC.mat t C