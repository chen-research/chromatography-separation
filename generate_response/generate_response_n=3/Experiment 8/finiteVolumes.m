function [tOut, COut, Tsol] = finiteVolumes(sysOpt, solOpt)

% Finite Volumes solver for the Equilbrium-Dispersive model (ED) and
% Transport-Dispersive (TD) model of chromatography

% Ye Zhang & Patrik Forssen
% June 2016
% Karlstad University

% Default input
if (nargin == 1)
  solOpt = [];
end
% Default output
tOut = [];
COut = [];

% Calculate some parameters
colS = sysOpt.D^2*pi/4.0;                         % Column surface [cm^2]
u    = (sysOpt.Flow/60.0)/(sysOpt.Epsilon*colS);  % Linear velocity [cm/s]
T0   = sysOpt.L/u;                                % Dead time [s]
F    = (1 - sysOpt.Epsilon)/sysOpt.Epsilon;       % Phase ratio

% n-Site Langmuir adsorption isotherm, calculate default time range
a      = sysOpt.IsoPar(:, 1:2:end);
nComp  = size(a, 1);    % Number of components injected
nSites = size(a, 2);    % Number of adsorption sites
Tsol   = 1.25*T0*(1 + F*max(sum(a, 2)));     % Default simulation time range [s]

% Calculate diffusion constant
Da   = (sysOpt.L*u)./(2*sysOpt.Nx);          % Diffusion constant
if (numel(Da) == 1)
  % Single value supplied
  Da = Da*ones(1, nComp);
end

% Parameters for square injection profile (boundary conditions)
Tinj  = (sysOpt.Vinj*1e-3)/(sysOpt.Flow/60); % Injection time [s]
% Tinj  = 10*60;
sysOpt.Vinj = Tinj*(sysOpt.Flow/60)*1e3;
dInj  = [Tinj - Tinj/100, Tinj + Tinj/100];
Csamp = 1e-3*sysOpt.Csamp;                   % Sample conc in [M]
if (numel(Csamp) == 1)
  % Single value supplied
  Csamp = Csamp*ones(1, nComp);
end

% Check if kinetic model (TD) should be used
kineticFlag = 0;
if (isfield(sysOpt, 'kf') && ~isempty(sysOpt.kf))
  kineticFlag = 1;
  if (numel(sysOpt.kf) == 1)
    % Single value supplied
    sysOpt.kf = sysOpt.kf*ones(1, nSites);
  end
  if (numel(sysOpt.kf) == nSites)
    % Vector supplied
    kf = repmat(sysOpt.kf(:)', nComp, 1);
  else
    kf = sysOpt.kf;
  end
end


% SOLVER

% Use concentration scaling, this fixes stability problems when very small
% amounts are injected
CScaleFac = 1./min(sysOpt.Vinj*1e-6*Csamp);
Csamp     = CScaleFac.*Csamp;

% Possible user override of default solution time range
if (~isempty(solOpt) && isstruct(solOpt) && isfield(solOpt, 'Tsol') && ...
    ~isempty(solOpt.Tsol) && solOpt.Tsol ~= 0)
  Tsol = solOpt.Tsol;
end

% Cells used
if (~isempty(solOpt) && isstruct(solOpt) && isfield(solOpt, 'xPoints') && ...
    ~isempty(solOpt.xPoints) && solOpt.xPoints ~= 0)
  nPts = solOpt.xPoints;
else
  % Default
  nPts = 200;
end

% Space increment
dz = sysOpt.L/(nPts + 1);

% ODE-solver options
% Use small intial time-step
options = odeset(...
  'OutputFcn'  , '', ...
  'InitialStep', Tinj/100, ...
  'NormControl', 'on');
if (kineticFlag)
  % Need Jacobian sparsity pattern for TD
  JPattern = TDJacPattern;
  options  = odeset(options, 'JPattern', JPattern);
else
  % Increase relative tolerance for ED as default
  options  = odeset(options, 'RelTol', 1e-4);
end
% Possible user override of norm control
if (~isempty(solOpt) && isstruct(solOpt) && ...
    isfield(solOpt, 'NormControl') && ~isempty(solOpt.NormControl))
  options  = odeset(options, 'NormControl', solOpt.NormControl);
end
% Possible user override of relative tolerance
if (~isempty(solOpt) && isstruct(solOpt) && isfield(solOpt, 'RelTol') && ...
    ~isempty(solOpt.RelTol) && solOpt.RelTol ~= 0)
  options  = odeset(options, 'RelTol', solOpt.RelTol);
end
% Possible user override of relative tolerance
if (~isempty(solOpt) && isstruct(solOpt) && isfield(solOpt, 'RelTol') && ...
    ~isempty(solOpt.RelTol) && solOpt.RelTol ~= 0)
  options  = odeset(options, 'RelTol', solOpt.RelTol);
end
% Possible user override of absolute tolerance
if (~isempty(solOpt) && isstruct(solOpt) && isfield(solOpt, 'AbsTol') && ...
    ~isempty(solOpt.AbsTol) && solOpt.AbsTol ~= 0)
  options  = odeset(options, 'AbsTol', solOpt.AbsTol);
end

% Pre-allocate common variables, speeds up calculations
aMat   = zeros(nPts, nComp, nSites);
bMat   = zeros(nPts, nComp, nSites);
aVec   = zeros(nComp, nSites);
bVec   = zeros(nComp, nSites);
for siteNo = 1 : nSites
  aSite              = sysOpt.IsoPar(: , 2*siteNo-1);
  aVec(:, siteNo)    = aSite;
  aMat(:, :, siteNo) = repmat(aSite', nPts, 1);
  bSite              = sysOpt.IsoPar(: , 2*siteNo);
  bVec(:, siteNo)    = bSite;
  bMat(:, :, siteNo) = repmat(bSite', nPts, 1);
end
multiIdMat = repmat(eye(nComp), [1, 1, nPts]);
zerosMat   = zeros(nPts-2, nComp);
thirdMat   = ones( nPts-2, nComp)/3;
twoMat     = 2*ones( nPts-2, nComp);
zerosVec   = zeros(1, nComp);
DadzMat    = repmat(Da/dz, nPts, 1);
if (kineticFlag)
  kfMat          = [];
  kfMat(1, :, :) = kf;
  kfMat          = repmat(kfMat, nPts, 1, 1);
end

% Initial conditions
Cinit = repmat(zeros(1, nComp), nPts, 1);
if (kineticFlag)
  % TD model also needs intial q
  qInit    = LangmuirIsotherm(Cinit, 'q');
  C_q_init = Cinit(:);
  for siteNo = 1 : nSites
    tmp = qInit(:, :, siteNo);
    C_q_init = [C_q_init; tmp(:)*CScaleFac];
  end
else
  % ED model
  Cinit = Cinit(:);
end


% Call ODE-solver
try
  if (kineticFlag)
    % Must use a stiff solver for TD. Note that a custom version of ode23t with
    % column reodering of the non-linear system Jacobian of ode23t is used in
    % the full implementation
    [tOut, COut] = ode23t(@odeDefTD, [0 Tsol], C_q_init, options);
    
  else
    
    % ED
    [tOut, COut] = ode45( @odeDefED, [0 Tsol], Cinit   , options);
  end
catch
  return
end

% Fix output
tOut        = tOut/60;
COut        = COut(:, nPts:nPts:nPts*nComp);
COut        = COut/CScaleFac;



  function dCdt = odeDefED(t, CIn)
    % ODE-defintion for ED. NOTE! CIn i SCALED concentrations!
    
    % Reshape input
    CMat = reshape(CIn, nPts, nComp);
    
    % Injection Profile
    tFac = 1 - (t - dInj(1))/(2*Tinj/100);
    if (t <= dInj(1))
      Cinj = Csamp(:)';
    elseif (t >= dInj(2))
      Cinj = zeros(1, nComp);
    else
      Cinj = Csamp(:)'*tFac;
    end
    
    % Koren flux limiter
    centerFluxMat = u*CMat;
    fluxGradMat   = diff(centerFluxMat);
    r             = (fluxGradMat(2:end, :) + 1e-10)./...
      (fluxGradMat(1:end-1, :) + 1e-10);
    theta         = max(zerosMat, min(2*r, min(thirdMat + 2*r/3, twoMat)));
    bndFluxMat    = centerFluxMat(2:end-1, :) + ...
      0.5*theta.*diff(centerFluxMat(1:end-1, :), 1);
    
    % Left boundary
    CLeft      = (u*dz*Cinj - 2*Da.*CMat(1, :))./(u*dz + 2*Da);
    bndFluxMat = [u*CLeft; centerFluxMat(1, :); bndFluxMat; ...
      centerFluxMat(end, :)];
    
    % dC/dz
    dCdz = [u*(CLeft - Cinj)./Da; diff(CMat)/dz; zerosVec];
    % RHS
    RHS  = -(diff(bndFluxMat, 1)/dz - DadzMat.*diff(dCdz, 1));
   
    % Jacobian
    Jac = LangmuirIsotherm(CMat/CScaleFac, 'dqdC');
    A   = multiIdMat + F*Jac;
             
    % Use sparse solver
    RHS  = RHS';
    RHS  = sparse(RHS(:));
    A    = blocks2sparse(A);
    dCdt = A\RHS;
    dCdt = reshape(dCdt, nComp, nPts);
    dCdt = dCdt'; 
    dCdt = dCdt(:);
  end



  function dCdt_dqdt = odeDefTD(t, C_q_In)
    %  ODE-defintion for TD. NOTE! C_q_In in SCALED concentrations!
    
    % Reshape input
    CMat = reshape(C_q_In(1 : nPts*nComp), nPts, nComp);
    qMat = zeros(nPts, nComp, nSites);
    for i = 1 : nSites
      qMat(:, :, i) = reshape(C_q_In(i*nPts*nComp+1 : (i+1)*nPts*nComp), ...
        nPts, nComp)/CScaleFac;
    end
    
    % Injection Profile
    tFac = 1 - (t - dInj(1))/(2*Tinj/100);
    if (t <= dInj(1))
      Cinj = Csamp(:)';
    elseif (t >= dInj(2))
      Cinj = zeros(1, nComp);
    else
      Cinj = Csamp(:)'*tFac;
    end
        
    % dq/dt
    qstarMat = LangmuirIsotherm(CMat/CScaleFac, 'q');
    dqdt     = zeros(nPts, nComp, nSites);
    dqdtSum  = zeros(nPts, nComp);
    for i = 1  : nSites
      dqdt(:, :, i) = kfMat(:, :, i).* ...
        (qstarMat(:, :, i) - qMat(:, :, i))*CScaleFac;
      dqdtSum = dqdtSum + dqdt(:, :, i); 
    end
        
    % Koren flux limiter
    centerFluxMat = u*CMat;
    fluxGradMat   = diff(centerFluxMat);
    r             = (fluxGradMat(2:end, :) + 1e-10)./ ...
      (fluxGradMat(1:end-1, :) + 1e-10);
    theta         = max(zerosMat, min(2*r, min(thirdMat + 2*r/3, twoMat)));
    bndFluxMat    = centerFluxMat(2:end-1, :) + ...
      0.5*theta.*diff(centerFluxMat(1:end-1, :), 1);
    
    % Left boundary
    CLeft      = (u*dz*Cinj - 2*Da.*CMat(1, :))./(u*dz + 2*Da);
    bndFluxMat = [u*CLeft; centerFluxMat(1, :); bndFluxMat; ...
      centerFluxMat(end, :)];
    
    % dC/dt
    dCdz = [u*(CLeft - Cinj)./Da; diff(CMat)/dz; zerosVec];
    dwdt = -(diff(bndFluxMat, 1)/dz - DadzMat.*diff(dCdz, 1));
    dCdt = dwdt - F*dqdtSum;
    
    % Output
    dCdt_dqdt = [dCdt(:); dqdt(:)];

  end



  function output = LangmuirIsotherm(CMat, mode)
    % Calcaultes stationary phase concentration, q, or the Jacobian of with
    % respect to the mobile phase concentration. Note that this part is coded as
    % a Fortran 90 MEX-file in the full implementation.
    
    switch mode
      case 'q'
        qstarMat = zeros(nPts, nComp, nSites);
        if (sysOpt.CompAds)
          
          % Competitive adosption
          for i = 1 : nSites
            denom = 1 + sum(bMat(:, :, i).*CMat, 2);
            qstarMat(:, :, i) = bsxfun(@rdivide, ...
              (aMat(:, :, i).*CMat), denom);
          end
          
        else
          
          % Non-competativ adsorption
          for i = 1 : nSites
            denom = 1 + bMat(:, :, i).*CMat;
            qstarMat(:, :, i) = aMat(:, :, i).*CMat./denom;
          end
          
        end
        output = qstarMat;
        
      case 'dqdC'
        der = zeros(nPts, nComp);
        Jac = zeros(nComp, nComp, nPts);
        if (sysOpt.CompAds)
         
          % Competitive adosption
          denom = zeros(nPts, nSites);
          for i = 1 : nSites
            denom(:, i) = 1 + sum(bMat(:, :, i).*CMat, 2);
            der = der + aMat(:, :, i).*bsxfun(@rdivide, bsxfun(@minus, ...
              denom(:, i), bMat(:, :, i).*CMat), denom(:, i).^2);
          end
          for i = 1 : nComp
            for j = 1 : nComp
              if (i ~= j)
                offDiagElem = zeros(nPts, 1);
                for k = 1 : nSites
                  offDiagElem = offDiagElem - aVec(i, k)* ...
                    bVec(j, k)*CMat(:, i)./denom(:, k).^2;
                end
                Jac(i, j, :) = offDiagElem;
              else
                Jac(i, j, :) = der(:, i);
              end
            end
          end
          
        else
          
          % Non-competativ adsorption
          for i = 1 : nSites
            denom = 1 + bMat(:, :, i).*CMat;
            der   = der + aMat(:, :, i)./denom - ...
              aMat(:, :, i).*bMat(:, :, i).*CMat./denom.^2;
          end
          for i = 1 : nComp
            Jac(i, i, :) = der(:, i);
          end
          
        end
        output = Jac;
        
    end
    
  end



  function B = blocks2sparse(A)
    % Converts a 3D-matrix to a block-diagonal sparse matrix
    
    blockSize = size(A, 1);
    nBlocks   = size(A, 3);
    
    indVec = 1 : blockSize : blockSize*nBlocks;
    rowInd = zeros(blockSize*blockSize, length(indVec));
    colInd = zeros(blockSize*blockSize, length(indVec));
    valVec = zeros(blockSize*blockSize, length(indVec));
    count  = 0;
    for j = 0 : blockSize-1
      for i = 0 : blockSize-1
        count  = count + 1;
        rowInd(count, :) = indVec+i;
        colInd(count, :) = indVec+j;
        valVec(count, :) = A(i+1, j+1, :);
      end
    end
    rowInd = rowInd(:);
    colInd = colInd(:);
    valVec = valVec(:);
    B      = sparse(rowInd, colInd, valVec, ...
      blockSize*nBlocks, blockSize*nBlocks);

  end



  function jPattern = TDJacPattern
    % Calculate the Jacobian sparsity pattern that is used
    % for the Transport-Dispersive Model of Chromatography
    
    % Default output
    jPattern = zeros(nPts*nComp*(nSites+1), nPts*nComp*(nSites+1));
    
    % Base blocks
    diag4Block = ...
      diag(ones(nPts-1, 1),  1) + diag(ones(nPts, 1)  ,  0) + ...
      diag(ones(nPts-1, 1), -1) + diag(ones(nPts-2, 1), -2);
    diag1Block = diag(ones(nPts, 1), 0);
    zeroBlock  = zeros(nPts);
    
    % Iterate to build matrix
    tmpMat = [];
    for siteBlockRow = 1 : nSites+1
      tmpMatRow = [];
      for siteBlockCol = 1 : nSites+1
        tmpCompMat = [];
        for compBlockRow = 1 : nComp
          tmpCompRow = [];
          for compBlockCol = 1 : nComp
            currBlock = zeroBlock;
            if (sysOpt.CompAds)
              
              % Competative adsorption
              if (siteBlockRow == 1 && siteBlockCol == 1 && ...
                  compBlockRow == compBlockCol)
                currBlock = diag4Block;
              elseif (siteBlockRow == 1 && ...
                  compBlockRow == compBlockCol)
                currBlock = diag1Block;
              elseif (siteBlockCol == 1)
                currBlock = diag1Block;
              elseif (siteBlockRow == siteBlockCol && ...
                  compBlockRow == compBlockCol)
                currBlock = diag1Block;
              end
              
            else
              
               % Non-competative adsorption
              if (siteBlockRow == 1 && siteBlockCol == 1 && ...
                  compBlockRow == compBlockCol)
                currBlock = diag4Block;
              elseif (siteBlockRow == 1 && ...
                  compBlockRow == compBlockCol)
                currBlock = diag1Block;
              elseif (siteBlockCol == 1 && ...
                  compBlockRow == compBlockCol)
                currBlock = diag1Block;
              elseif (siteBlockRow == siteBlockCol && ...
                  compBlockRow == compBlockCol)
                currBlock = diag1Block;  
              end
              
            end
            % Add component columnn
            tmpCompRow = [tmpCompRow, currBlock];
          end
          % Add component row
          tmpCompMat = [tmpCompMat; tmpCompRow];
        end
        % Add site columnn
        tmpMatRow = [tmpMatRow, tmpCompMat];
      end
      % Add site row
      tmpMat = [tmpMat; tmpMatRow];
    end
    
    % Set sparsity pattern
    jPattern(1:size(tmpMat, 1), 1:size(tmpMat, 2)) = tmpMat;
    jPattern = sparse(jPattern);
  end

end