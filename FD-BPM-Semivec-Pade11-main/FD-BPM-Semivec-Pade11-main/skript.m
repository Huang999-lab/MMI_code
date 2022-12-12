%% Create the concentration profile of an open MMI-Element
% Exemplary BPM calculation for a rectangular step-index MMI-element

% Parameter definition
clear
close all
clc

PROGRESS = 'cl';

lambda              = 1330e-9;       % Wavelength
alpha               = 0.5;           % Crank-Nicolson-Constant
delta               = 1e-6;
alpha               = alpha + delta;
solver_tolerance    = 1e-12;         % Tolerance of BiCGM solver

POLARIZATION = 'TE';
FIELDCOMPONENT = 'Ex';
BC = 'ABC';
ABSORBER = 0;                     % Small value needed to prevent diverging field, 0 to deactivate.
    
dx_fine   = .5e-6;                   % Step size for fine step
dy_fine   = .5e-6;                   % Step size for fine step
dx_coarse = 1e-6;                    % Step coarse for fine step               
dy_coarse = 1e-6;                    % Step coarse for fine step
dz        = 1e-6;                    % Propagationstep

Excitation.type = 'gauss';           % Type of Excitation
Excitation.visualize = 1;            % Plot Excitation field
Excitation.sigma_x = 5e-6;           % Only needed for fieldtype 'gauss'
Excitation.sigma_y = 3e-6;           % Only needed for fieldtype 'gauss'
Excitation.threshold = 1-1/exp(1);   % Only needed for fieldtype 'full'
Excitation.numberInterpolation = 0;  % Only needed for fieldtype 'modal'
Excitation.mode = 'k_bar';           % Only needed for fieldtype 'modal'

% Some possible manipulations of index profile
append = 0;
append_dist = 0e-6;     % um
prepend = 0;
prepend_dist = 0e-6;    % um
upper_cladding = 1;

%% Generation of the refractive index profile
% Loading index profile and base vectors
load('Index_GI_MM_Diffused_Waveguide.mat')
n_max = max(max(n));
n_min = min(min(n));
nc = repmat(n,[1 1 101]); % 重复三维块数组副本
beta_0 = 2*pi/lambda;       % Wave number
k_bar  = (n_max - n_min)/2; % Mean wave number

% Grid
z = [0:1e-6:100e-6];                                                                % Base vector
[xg,yg,zg] = meshgrid(x,y,z);                                                       % Base grid 对于 136x355x101 个点都有对应的 xg,yg,zg 坐标值

xgb = squeeze(xg(:,:,1));  % Transversal grid of first slice
ygb = squeeze(yg(:,:,1));  % Transversal grid of first slice
%% BPM flags and parameters
dim_y   = size(nc,1); % Global dimension
dim_x   = size(nc,2); 

dim_yl   = dim_y - 2; % Local dimensions (without boundary values)
dim_xl   = dim_x - 2;

grid.localIndex = zeros(size(nc,1),size(nc,2));   % Grid variables
grid.globalIndex = zeros(size(nc,1),size(nc,2));

grid.boundaryFlag = zeros(size(nc,1),size(nc,2)); % Boundary flag
grid.boundaryFlag([1 end],[1:end 1:end]) = 1;
grid.boundaryFlag([1:end 1:end],[1 end]) = 1;

grid.localIndex(2:end-1,2:end-1) = reshape(1:1:dim_xl*dim_yl',dim_yl,dim_xl);% 节点添加标签
grid.globalIndex(1:end) = 1:1:length(grid.globalIndex(1:end)); % 节点添加标签
%% Execute BPM
out = [POLARIZATION '-BPM with ' BC ' boundary condition started for graded-index multimode diffused waveguide at 1330nm wavelength.'];
disp(out);

%%%% Definition of variables and constants
format long
beta_0 = 2*pi/lambda;           % wave number
beta_z = beta_0 * k_bar;         % propagation constant as defined by neff
n_max_b = max(max(nc(:,:,1)));     % maximum value of n
n_min_b = min(min(nc(:,:,1)));     % minimum value of n
delta_n_b = n_max_b - n_min_b;        % maximum refractive index contrast

n_max = squeeze(max(max(nc(:,:,:))));
n_min = squeeze(min(min(nc(:,:,:))));
delta_n = n_max - n_min;

%%%% Check field components input: FX --> Ex

%%%% BPM specific parameters
globalIndexSlgs   = grid.globalIndex(2:end-1,2:end-1); % Vector for global adresses of all lokal elements
globalAdrSlgs   = reshape(globalIndexSlgs,size(globalIndexSlgs,1)*size(globalIndexSlgs,2),1); % 地址列向量

%%%% Definition of initial fields depending on chosen excitation: gaussian
sigma_x = Excitation.sigma_x; % sigma of gaussian distribution as measure for the extention of the gaussian beam in x-direction in [um]
sigma_y = Excitation.sigma_y; % sigma of gaussian distribution as measure for the extention of the gaussian beam in y-direction in [um]

xg1 = squeeze(xg(:,:,1));
yg1 = squeeze(yg(:,:,1));

[r_max,c_max] = find(squeeze(nc(:,:,1)) == max(max(nc(:,:,1)))); % Those parameters assure, that the gaussian beam hits the profile at its maximum preventing mismatch

xg1 = xg1 - xg1(r_max(1),c_max(1)); % x-grid for the definition of the gaussian beam
yg1 = yg1 - yg1(r_max(1),c_max(1)); % y-grid for the definition of the gaussian beam

phiInput = 1*exp(-xg1.^2/(2*(sigma_x))^2 -yg1.^2/(2*(sigma_y))^2);  % Definition of gaussian beam
phi(:,:,1) = phiInput;
% figure
% surface(xg1,yg1,phiInput)

%%%% Defining gamma for handling different boundary conditions: ABC
gammaBoundaryCondition = 0;

%%%% Generate multistep paramters
UX = zeros(1,1);
VX = zeros(1,1);

R=2*beta_z;

C=1/R^2-1j*dz*(1-alpha)*1/R;
D=1/R^2+1j*dz*( alpha )*1/R;

UX(1)= C;
VX(1)= D;
tic
c = 1; % Global progress counter (waitbar)

%% Propagation in z direction
    %%%% Address space initialization
    globalIndexSlgs   = grid.globalIndex(2:end-1,2:end-1); %(without boundary values)
    globalAdrSlgs   = reshape(globalIndexSlgs,size(globalIndexSlgs,1)*size(globalIndexSlgs,2),1);% 地址列向量

    globalIndexN = grid.globalIndex(3:end-1,2:end-1); 
    globalAdrN = reshape(globalIndexN,size(globalIndexN,1)*size(globalIndexN,2),1);
    
    globalIndexS = grid.globalIndex(2:end-2,2:end-1); 
    globalAdrS = reshape(globalIndexS,size(globalIndexS,1)*size(globalIndexS,2),1);

    globalIndexE = grid.globalIndex(2:end-1,2:end-2);
    globalAdrE = reshape(globalIndexE,size(globalIndexE,1)*size(globalIndexE,2),1);
    
    globalIndexW = grid.globalIndex(2:end-1,3:end-1);
    globalAdrW = reshape(globalIndexW,size(globalIndexW,1)*size(globalIndexW,2),1);
    
    localIndexSlgs    = grid.localIndex(2:end-1,2:end-1);
    localAdrSlgs    = reshape(localIndexSlgs,size(localIndexSlgs,1)*size(localIndexSlgs,2),1);
    
    localIndexN  = grid.localIndex(3:end-1,2:end-1); 
    localAdrN  = reshape(localIndexN,size(localIndexN,1)*size(localIndexN,2),1);
    
    localIndexS  = grid.localIndex(2:end-2,2:end-1); 
    localAdrS  = reshape(localIndexS,size(localIndexS,1)*size(localIndexS,2),1);
    
    localIndexE = grid.localIndex(2:end-1,2:end-2);
    localAdrE = reshape(localIndexE,size(localIndexE,1)*size(localIndexE,2),1);
    
    localIndexW = grid.localIndex(2:end-1,3:end-1);
    localAdrW = reshape(localIndexW,size(localIndexW,1)*size(localIndexW,2),1);

    %%
for kz = 1:1:size(nc,3)-1  % 一层一层迭代

    % Extract known field values for current propagation step

    pb = phi(:,:,kz);

    %%%% Call functions for calculation of diagonals :only for TE,Ex, "BC" = ABC
%%   [ diagC,diagN,diagS,diagE,diagW ]      = diagonalsPade(beta_0,neff,nc(:,:,kz+1),xg(:,:,kz+1),yg(:,:,kz+1),dim_y,dim_xl,dim_yl,grid,gammaBoundaryCondition,POLARIZATION,FX,BC);

%%%%diagC = padeX(n,xg,dim_xl,dim_yl,dim_y,localAdrSlgs,globalAdrSlgs,POLARIZATION,FIELDCOMPONENT) + padeY(n,yg,dim_xl,dim_yl,localAdrSlgs,globalAdrSlgs,POLARIZATION,FIELDCOMPONENT) + beta_0^2 .* (epsilon(globalAdrSlgs) - neff^2);
    a = xg(:,:,kz+1);
    de = a(globalAdrSlgs + dim_y) - a(globalAdrSlgs);% delta e in east direction 东
    dw = a(globalAdrSlgs) - a(globalAdrSlgs - dim_y);% delta w in west direction 西
    % Round off error
    de = 1e-12*round(de*1e12);
    dw = 1e-12*round(dw*1e12);
    epsilon = nc(:,:,kz+1).*nc(:,:,kz+1);% epsilon := n
    aWest = 2./(dw.*(dw+de)) .* (2*epsilon(globalAdrSlgs - dim_y) ./ (epsilon(globalAdrSlgs) + epsilon(globalAdrSlgs - dim_y))); 
    aWestReshaped = zeros(dim_xl*dim_yl,1);
    aWestReshaped(localAdrSlgs) = aWest;
    
    aEast = (2./(de.*(de+dw))) .* (2*epsilon(globalAdrSlgs + dim_y) ./ (epsilon(globalAdrSlgs) + epsilon(globalAdrSlgs + dim_y)));
    aEastReshaped = zeros(dim_xl*dim_yl,1);
    aEastReshaped(localAdrSlgs) = aEast;
    
    aX = -4./(de.*dw) + aEastReshaped + aWestReshaped;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a = yg(:,:,kz+1);
    dn = a(globalAdrSlgs - 1) - a(globalAdrSlgs);
    ds = a(globalAdrSlgs) - a(globalAdrSlgs + 1);
    
    % Round off error
    dn = 1e-12*round(dn*1e12);
    ds = 1e-12*round(ds*1e12); 
    
    aNorth = (2./(dn.*(dn+ds))) .* ones(length(globalAdrSlgs),1);
    aNorthReshaped = zeros(dim_xl*dim_yl,1);
    aNorthReshaped(localAdrSlgs) = aNorth;
    
    a_s = (2./(ds.*(ds+dn))) .* ones(length(globalAdrSlgs),1);
    aSouthReshaped = zeros(dim_xl*dim_yl,1);
    aSouthReshaped(localAdrSlgs) = a_s;
    
    aY = - aNorthReshaped - aSouthReshaped; 

    diagC = aX + aY + beta_0^2 .* (epsilon(globalAdrSlgs) - k_bar^2);
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        diagN = padeN(beta_0,yg(:,:,kz+1),dim_xl,dim_yl,localAdrN,globalAdrN,POLARIZATION,FIELDCOMPONENT);
    a = yg(:,:,kz+1);
    dn = a(globalAdrN - 1) - a(globalAdrN);
    ds = a(globalAdrN) - a(globalAdrN + 1);
    % Round off error
    dn = 1e-12*round(dn*1e12);
    ds = 1e-12*round(ds*1e12);
    
    aNorth = (2./(dn.*(dn+ds))) .* ones(length(globalAdrN),1);
    aNorthReshaped = zeros(dim_xl*dim_yl,1);
    aNorthReshaped(localAdrN) = aNorth;
    diagN = aNorthReshaped;

%%%%       diagS = padeS(beta_0,yg(:,:,kz+1),dim_xl,dim_yl,localAdrS,globalAdrS,POLARIZATION,FIELDCOMPONENT);% aSorthReshaped
    a = yg(:,:,kz+1);
    dn = a(globalAdrN - 1) - a(globalAdrN);
    ds = a(globalAdrN) - a(globalAdrN + 1);
    a_s = (2./(ds.*(ds+dn))) .* ones(length(globalAdrS),1);
    aSouthReshaped = zeros(dim_xl*dim_yl,1);
    aSouthReshaped(localAdrS) = a_s;
    diagS = aSouthReshaped;
    

%%%%     diagE = padeE(beta_0,xg(:,:,kz+1),dim_y,dim_xl,dim_yl,localAdrE,globalAdrE,POLARIZATION,FIELDCOMPONENT);% aEastReshaped
    a = xg(:,:,kz+1);
    de = a(globalAdrE + dim_y) - a(globalAdrE);
    dw = a(globalAdrE) - a(globalAdrE - dim_y);
    
    % Round off error 
    de = 1e-12*round(de*1e12);
    dw = 1e-12*round(dw*1e12);
  
    aEast = (2./(de.*(de+dw))) .* (2*epsilon(globalAdrE + dim_y) ./ (epsilon(globalAdrE) + epsilon(globalAdrE + dim_y)));
    aEastReshaped = zeros(dim_xl*dim_yl,1);
    aEastReshaped(localAdrE) = aEast;
    diagE = aEastReshaped;
    
%%%%     diagW = padeW(beta_0,xg(:,:,kz+1),dim_y,dim_xl,dim_yl,localAdrW,globalAdrW,POLARIZATION,FIELDCOMPONENT);% aWestReshaped
    a = xg(:,:,kz+1);
    de = a(globalAdrE + dim_y) - a(globalAdrE);
    dw = a(globalAdrE) - a(globalAdrE - dim_y);
     % Round off error 
    de = 1e-12*round(de*1e12);
    dw = 1e-12*round(dw*1e12);  
    
    aWest = 2./(dw.*(dw+de)) .* (2*epsilon(globalAdrW - dim_y) ./ (epsilon(globalAdrW) + epsilon(globalAdrW - dim_y)));
    aWestReshaped = zeros(dim_xl*dim_yl,1);
    aWestReshaped(localAdrW) = aWest;
    diagW = aWestReshaped;
    
    
    
    
%% [ diagBC,diagBN,diagBS,diagBE,diagBW ] = diagonalsPade(beta_0,neff,nc(:,:,kz),xg(:,:,kz),yg(:,:,kz),dim_y,dim_xl,dim_yl,grid,gammaBoundaryCondition,POLARIZATION,FX,BC);
    a = xg(:,:,kz);
    de = a(globalAdrSlgs + dim_y) - a(globalAdrSlgs);
    dw = a(globalAdrSlgs) - a(globalAdrSlgs - dim_y);
    % Round off error
    de = 1e-12*round(de*1e12);
    dw = 1e-12*round(dw*1e12);
    epsilon = nc(:,:,kz).*nc(:,:,kz);
    aWest = 2./(dw.*(dw+de)) .* (2*epsilon(globalAdrSlgs - dim_y) ./ (epsilon(globalAdrSlgs) + epsilon(globalAdrSlgs - dim_y))); 
    aWestReshaped = zeros(dim_xl*dim_yl,1);
    aWestReshaped(localAdrSlgs) = aWest;
    
    aEast = (2./(de.*(de+dw))) .* (2*epsilon(globalAdrSlgs + dim_y) ./ (epsilon(globalAdrSlgs) + epsilon(globalAdrSlgs + dim_y)));
    aEastReshaped = zeros(dim_xl*dim_yl,1);
    aEastReshaped(localAdrSlgs) = aEast;
    
    aX = -4./(de.*dw) + aEastReshaped + aWestReshaped;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a = yg(:,:,kz);
    dn = a(globalAdrSlgs - 1) - a(globalAdrSlgs);
    ds = a(globalAdrSlgs) - a(globalAdrSlgs + 1);
    
    % Round off error
    dn = 1e-12*round(dn*1e12);
    ds = 1e-12*round(ds*1e12); 
    
    aNorth = (2./(dn.*(dn+ds))) .* ones(length(globalAdrSlgs),1);
    aNorthReshaped = zeros(dim_xl*dim_yl,1);
    aNorthReshaped(localAdrSlgs) = aNorth;
    
    a_s = (2./(ds.*(ds+dn))) .* ones(length(globalAdrSlgs),1);
    aSouthReshaped = zeros(dim_xl*dim_yl,1);
    aSouthReshaped(localAdrSlgs) = a_s;
    
    aY = - aNorthReshaped - aSouthReshaped; 

    diagBC = aX + aY + beta_0^2 .* (epsilon(globalAdrSlgs) - k_bar^2);
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        diagN = padeN(beta_0,yg(:,:,kz+1),dim_xl,dim_yl,localAdrN,globalAdrN,POLARIZATION,FIELDCOMPONENT);
    a = yg(:,:,kz);
    dn = a(globalAdrN - 1) - a(globalAdrN);
    ds = a(globalAdrN) - a(globalAdrN + 1);
    % Round off error
    dn = 1e-12*round(dn*1e12);
    ds = 1e-12*round(ds*1e12);
    
    aNorth = (2./(dn.*(dn+ds))) .* ones(length(globalAdrN),1);
    aNorthReshaped = zeros(dim_xl*dim_yl,1);
    aNorthReshaped(localAdrN) = aNorth;
    diagBN = aNorthReshaped;

%%%%       diagS = padeS(beta_0,yg(:,:,kz+1),dim_xl,dim_yl,localAdrS,globalAdrS,POLARIZATION,FIELDCOMPONENT);% aSorthReshaped
    a = yg(:,:,kz);
    dn = a(globalAdrN - 1) - a(globalAdrN);
    ds = a(globalAdrN) - a(globalAdrN + 1);
    a_s = (2./(ds.*(ds+dn))) .* ones(length(globalAdrS),1);
    aSouthReshaped = zeros(dim_xl*dim_yl,1);
    aSouthReshaped(localAdrS) = a_s;
    diagBS = aSouthReshaped;
    

%%%%     diagE = padeE(beta_0,xg(:,:,kz+1),dim_y,dim_xl,dim_yl,localAdrE,globalAdrE,POLARIZATION,FIELDCOMPONENT);% aEastReshaped
    a = xg(:,:,kz);
    de = a(globalAdrE + dim_y) - a(globalAdrE);
    dw = a(globalAdrE) - a(globalAdrE - dim_y);
    
    % Round off error 
    de = 1e-12*round(de*1e12);
    dw = 1e-12*round(dw*1e12);
  
    aEast = (2./(de.*(de+dw))) .* (2*epsilon(globalAdrE + dim_y) ./ (epsilon(globalAdrE) + epsilon(globalAdrE + dim_y)));
    aEastReshaped = zeros(dim_xl*dim_yl,1);
    aEastReshaped(localAdrE) = aEast;
    diagBE = aEastReshaped;
    
%%%%     diagW = padeW(beta_0,xg(:,:,kz+1),dim_y,dim_xl,dim_yl,localAdrW,globalAdrW,POLARIZATION,FIELDCOMPONENT);% aWestReshaped
    a = xg(:,:,kz);
    de = a(globalAdrE + dim_y) - a(globalAdrE);
    dw = a(globalAdrE) - a(globalAdrE - dim_y);
     % Round off error 
    de = 1e-12*round(de*1e12);
    dw = 1e-12*round(dw*1e12);  
    
    aWest = 2./(dw.*(dw+de)) .* (2*epsilon(globalAdrW - dim_y) ./ (epsilon(globalAdrW) + epsilon(globalAdrW - dim_y)));
    aWestReshaped = zeros(dim_xl*dim_yl,1);
    aWestReshaped(localAdrW) = aWest;
    diagBW = aWestReshaped;
    
    

    %% Apply multistep method

    for ii = 1:1:1  % Padé(1,1) only requires one multistep

        % Merge diagonals in system matrix

        A = sparse((size(nc,1)-2)*(size(nc,2)-2),(size(nc,1)-2)*(size(nc,2)-2)); %sparse(m,n) 生成 m×n 全零稀疏矩阵。
        A = spdiags(1 + VX(ii)*diagC,0,A); % S = spdiags(Bin,d,A) 将 d 指定的 A 中的对角线替换为 Bin 的列。 spy(A)
        A = spdiags([VX(ii)*diagN(2:end); 0],-1,A);
        A = spdiags([0; VX(ii)*diagS(1:end-1)],1,A);
        A = spdiags([zeros(dim_yl,1); VX(ii)*diagE(1:end-dim_yl)],dim_yl,A);
        A = spdiags([VX(ii)*diagW(dim_yl+1:end); zeros(dim_yl,1)],-dim_yl,A);
%         spy(A)

        % Compute right side of equation

        diagBC = 1 + UX(ii)*diagBC;
        diagBN = [UX(ii)*diagBN(2:end); 0];
        diagBS = [0; UX(ii)*diagBS(1:end-1)];
        diagBE = [zeros(dim_yl,1); UX(ii)*diagBE(1:end-dim_yl)];
        diagBW = [UX(ii)*diagBW(dim_yl+1:end); zeros(dim_yl,1)];

        % Multiply right side with known field values at 'kz'

        diagBC = pb(globalAdrSlgs)          .* diagBC;
        diagBN = pb(globalAdrSlgs - 1)      .* diagBN;
        diagBS = pb(globalAdrSlgs + 1)      .* diagBS;
        diagBE = pb(globalAdrSlgs + dim_y)  .* diagBE;
        diagBW = pb(globalAdrSlgs - dim_y)  .* diagBW;

        % Form vector for right side

        b = sparse(diagBC + diagBN + diagBS + diagBE + diagBW);

        %% Solve SoLE

        [phiSolution,~] = bicgstab(A,b,solver_tolerance);

        % Reshape solution to match field matrix

        pb = zeros(dim_y,dim_x);
        pb(2:end-1,2:end-1) = reshape(phiSolution,dim_yl,dim_xl);

        %% Apply absorber if specified

        if strcmp(BC,'ABC')

            if isnumeric(ABSORBER) && (ABSORBER > 0) && (ABSORBER < 1) && (ABSORBER ~= 0)

                nThreshold             = n_min(kz) + delta_n(kz) * ABSORBER;
                adr_n_threshold         = squeeze(nc(:,:,kz)) <= nThreshold;
                pb(adr_n_threshold)    = 0;

            elseif ~isnumeric(ABSORBER) || (ABSORBER < 0) || (ABSORBER > 1)

                out = 'Invalid specification of absorber: 0 <= absorber < 1.';
                disp(out)
                return

            end

        end

    end

    %% Merge solution in global field matrix

    phi(:,:,kz+1) = pb;

    % Update progress bar

    if strcmp(PROGRESS,'bar') == 1

        stepPercent = 1;  % Defines step size for progress bar update

        if floor(100*kz/(size(nc,3)-1)) > c*stepPercent

          if ishandle(h)    % Check if waitbar has not been closed to prevent error

            minutesPassed      = floor(toc/60);
            secondsPassed      = toc - minutesPassed*60;
            minutesRemaining   = floor((toc/c)*((100/stepPercent)-c)/60);
            secondsRemaining   = (toc/c)*((100/stepPercent)-c) - minutesRemaining*60;

            waitbar(kz/(size(nc,3)-1),h,['Progress: ' num2str(c*stepPercent) '%. Time remaining: ' num2str(floor(minutesRemaining)) 'm ' num2str(ceil(secondsRemaining)) 's.'])

            c = c + 1;

          else

             disp('Process interrupted')
             break

          end
        end

    elseif strcmp(PROGRESS,'cl') == 1

        stepPercent = 10;  % Defines step size for progress bar update

        if floor(100*kz/(size(nc,3)-1)) > c*stepPercent

            minutesPassed      = floor(toc/60);
            secondsPassed      = toc - minutesPassed*60;
            minutesRemaining   = floor((toc/c)*((100/stepPercent)-c)/60);
            secondsRemaining   = (toc/c)*((100/stepPercent)-c) - minutesRemaining*60;

            out = ['   Progress: ' num2str(c*stepPercent) '%. Time remaining: ' num2str(floor(minutesRemaining)) 'm ' num2str(ceil(secondsRemaining)) 's.'];
            disp(out)

            c = c + 1;

        end

    end
end




%% Visualization

Ex_int = squeeze(sum(abs(phi) .* abs(yg)));
figure
surf(z(2:end),x,squeeze(Ex_int(:,2:end))/max(max(Ex_int(:,2:end))))
shading flat
xlabel('z [\mum]')
ylabel('x [\mum]')
title('Normalized absolute E-field (after y-integration) [a.u.]')
axis([z(2) z(end) x(1) x(end) 0 1])









