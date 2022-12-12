%% Create the concentration profile of an open MMI-Element
% Exemplary BPM calculation for a rectangular step-index MMI-element

% Parameter definition
clear
close all
clc

PROGRESS = 'cl';

lambda              = 1330e-9;       % Wavelength
alpha               = 0.5;           % Crank-Nicolson-Constant (statbility Condition) alpha=0.5: oscillation.
delta               = 1e-6;
alpha               = alpha + delta;
solver_tolerance    = 1e-12;         % Tolerance of BiCGM solver

POLARIZATION = 'TE';
FIELDCOMPONENT = 'Ex';
BC = 'TBC';
ABSORBER = 0;                     % Small value needed to prevent diverging field, 0 to deactivate.
    
dx_fine   = .5e-6;                   % Step size for fine step
dy_fine   = .5e-6;                   % Step size for fine step
dx_coarse = 1e-6;                    % Step coarse for fine step               
dy_coarse = 1e-6;                    % Step coarse for fine step
dz        = 1e-6;                    % Propagationstep µm(1e-6)

Excitation.type = 'gauss';           % Type of Excitation
Excitation.visualize = 1;            % Plot Excitation field
Excitation.sigma_x = 5e-6;           % Only needed for fieldtype 'gauss'
Excitation.sigma_y = 3e-6;           % Only needed for fieldtype 'gauss'
Excitation.threshold = 1-1/exp(1);   % Only needed for fieldtype 'full'
Excitation.numberInterpolation = 0;  % Only needed for fieldtype 'modal'
Excitation.mode = 'k_bar';           % Only needed for fieldtype 'modal'


out = 'Loading index profile...';
disp(out)
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading index profile and base vectors
load('Index_GI_MM_Diffused_Waveguide.mat')
x = [0: dx_coarse: x(end)];
y = [0: dy_coarse: y(end)];
% change index profile to Luft profile: 
n = ones(length(y),length(x)); %  !!!!!!!!!!!!!!!!!!!!!!!!!!



% Grid
z = [0:dz:100e-6];                                                                % Base vector
[xg,yg,zg] = meshgrid(x,y,z);                                                       % Base grid

xgb = squeeze(xg(:,:,1));  % Transversal grid of first slice
ygb = squeeze(yg(:,:,1));  % Transversal grid of first slice




n_max = max(max(n));
n_min = min(min(n));
nc = repmat(n,[1 1 length(z)]);
beta_0 = 2*pi/lambda;       % Wave number
k_bar  = (n_max - n_min)/2; % Mean wave number
k_bar  = 1; % Mean wave number only for Luft !!!!!!!!!!!!!!!!!!!!!!!



out = ['  Index profile loaded: ' num2str(toc) 's.'];
disp(out)
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Some possible manipulations of index profile
append = 0;
append_dist = 0e-6;     % um
prepend = 0;
prepend_dist = 0e-6;    % um
upper_cladding = 1;

% Append exit
% Can be used to replicate the last slice

out = 'Prepending and/or appending first and last slice...';
disp(out)
tic

if append == 1

    append_steps = floor(append_dist/dz);

    nc_end  = repmat(squeeze(nc(:,:,end)),[1 1 append_steps]);
    nc      = cat(3,nc,nc_end); 

    z_end  = z(end)+dz:dz:z(end)+append_dist;
    z      = [z z_end]; 

    [xg,yg,zg] = meshgrid(x,y,z);

end

% Prepend input
% Can be used to replicate the first slice

if prepend == 1

    prepend_steps = floor(prepend_dist/dz);

    nc_input = repmat(squeeze(nc(:,:,1)),[1 1 prepend_steps]);
    nc       = cat(3,nc_input,nc);

    z_input = 0:dz:prepend_dist-dz;
    z = [z_input z+prepend_dist];

    [xg,yg,zg] = meshgrid(x,y,z);

end

out = ['  Prepending/appending done: ' num2str(toc) 's.'];
disp(out)
tic

% Attaching upper cladding 连接上覆层
    
if upper_cladding == 1
    
    out = 'Attaching upper cladding...';
    disp(out)
    tic

    % Generating new vectors and grids
    yCladd = y(end)+dy_fine:dy_fine:y(end)+10e-6;
    y      = [y yCladd];
    [xg,yg,zg] = meshgrid(x,y,z);

    % Removing round-off errors
    xg = 1e-12*round(xg*1e12);
    yg = 1e-12*round(yg*1e12);
    zg = 1e-12*round(zg*1e12);
     
    % Attaching upper cladding to existing index profile
    nCladd = 1.522*ones(length(yCladd),length(x),length(z));

    nc = cat(1,nc,nCladd); % 垂直串联这两个矩阵

    out = ['  Upper Cladding attached: ' num2str(toc) 's.'];
    disp(out)

end
 
%% BPM flags and parameters

out = 'Defining variables and BPM flags...';
disp(out)
tic

dim_y   = size(nc,1); % Global dimension
dim_x   = size(nc,2); 

dim_yl   = dim_y - 2; % Local dimensions (without boundary values)
dim_xl   = dim_x - 2;

grid.localIndex = zeros(size(nc,1),size(nc,2));   % Grid variables
grid.globalIndex = zeros(size(nc,1),size(nc,2));

grid.boundaryFlag = zeros(size(nc,1),size(nc,2)); % Boundary flag
grid.boundaryFlag([1 end],[1:end 1:end]) = 1;
grid.boundaryFlag([1:end 1:end],[1 end]) = 1;

grid.localIndex(2:end-1,2:end-1) = reshape(1:1:dim_xl*dim_yl',dim_yl,dim_xl);% 内部添加标签（边界除外）
grid.globalIndex(1:end) = 1:1:length(grid.globalIndex(1:end));% 所有添加标签

out = ['  Defined necessary parameters: ' num2str(toc) 's.'];
disp(out)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[Ex,~,~,~] = FDBPMPade11Semivec(nc,lambda,k_bar,alpha,solver_tolerance,xg,yg,dz,Excitation,'TE',{'Ex';'Hy'},BC,ABSORBER,PROGRESS);
%% Definition of variables and constants

format long
beta_0 = 2*pi/lambda;           % wave number
beta_z = beta_0 * k_bar;         % propagation constant as defined by neff
n_max_b = max(max(nc(:,:,1)));     % maximum value of n
n_min_b = min(min(nc(:,:,1)));     % minimum value of n
delta_n_b = n_max_b - n_min_b;        % maximum refractive index contrast

n_max = squeeze(max(max(nc(:,:,:))));
n_min = squeeze(min(min(nc(:,:,:))));
delta_n = n_max - n_min;

% Check field components input
[FX,~] = checkFieldComponents({'Ex';'Hy'}); % Function checks which field component is to be evaluated

%% BPM specific parameters

dim_y   = size(nc,1); % Number of global elements in y direction
dim_x   = size(nc,2); % Number of global elements in x direction

dim_xl  = size(nc,2) - 2; % Number of local elements in y direction (without the boundary of the structure)
dim_yl  = size(nc,1) - 2; % Number of local elements in x direction (without the boundary of the structure)

phi = zeros(size(nc,1),size(nc,2),size(nc,3)); % field matrix

grid.localIndex = zeros(size(nc,1),size(nc,2)); % Grid that speficies lokal element numbers
grid.localIndex(2:end-1,2:end-1) = reshape(1:1:dim_xl*dim_yl',dim_yl,dim_xl);

grid.globalIndex = zeros(size(nc,1),size(nc,2)); % Grid that speficies global element numbers. Note that in the global case (looking at the outer dimensions of the grid that include the boundary) there is no difference between 'Index' and 'Adress'. For this reason 'Index' is used in matrix form and 'Adr' is used as vector.
grid.globalIndex(1:end) = 1:1:length(grid.globalIndex(1:end));

grid.boundaryFlag = zeros(size(nc,1),size(nc,2));  % Grid that flags boundary elements of the structure
grid.boundaryFlag([1 end],[1:end 1:end]) = 1;
grid.boundaryFlag([1:end 1:end],[1 end]) = 1;

globalIndexSlgs   = grid.globalIndex(2:end-1,2:end-1); % Vector for global adresses of all lokal elements
globalAdrSlgs   = reshape(globalIndexSlgs,size(globalIndexSlgs,1)*size(globalIndexSlgs,2),1);

%% Definition of initial fields depending on chosen excitation
sigma_x = Excitation.sigma_x; % sigma of gaussian distribution as measure for the extention of the gaussian beam in x-direction in [um]
sigma_y = Excitation.sigma_y; % sigma of gaussian distribution as measure for the extention of the gaussian beam in y-direction in [um]


xg1 = squeeze(xg(:,:,1));
yg1 = squeeze(yg(:,:,1));

[r_max,c_max] = find(squeeze(nc(:,:,1)) == max(max(nc(:,:,1)))); % Those parameters assure, that the gaussian beam hits the profile at its maximum preventing mismatch
% xg1 = xg1 - xg1(r_max(1),c_max(1)); % x-grid for the definition of the gaussian beam
% yg1 = yg1 - yg1(r_max(1),c_max(1)); % y-grid for the definition of the gaussian beam

% change parameter to Luft profile: the gaussian beam hits the profile at xg=100.yg=50  
xg1 = xg1 - xg1(50,100); % x-grid for the definition of the gaussian beam
yg1 = yg1 - yg1(50,100); % y-grid for the definition of the gaussian beam

phiInput = 1*exp(-xg1.^2/(2*(sigma_x))^2 -yg1.^2/(2*(sigma_y))^2);  % Definition of gaussian beam
phi(:,:,1) = phiInput;                                              % Merging gaussian beam in global field matrix


figure
surf(xg(:,:,1),yg(:,:,1),abs(phi(:,:,1)))
xlabel('x')
ylabel('y')
title('Exciting field distribution')
shading flat


%% Defining gamma for handling different boundary conditions

if strcmp(BC,'ABC')

    gammaBoundaryCondition = 0;

elseif strcmp(BC,'TBC')

    gammaBoundaryCondition = {};

else

    out = 'Invalid specification of boundary condition. Possible choices are ''ABC'' or ''TBC''.';
    disp(out)
    return

end
%% Generate multistep paramters
[UX,VX] = genMultistepVars11(dz,alpha,beta_z);

tic
c = 1; % Global progress counter (waitbar)

if strcmp(PROGRESS,'bar') == 1

    h = waitbar(0,'','Name',['Computing Padé ' POLARIZATION '-BPM with ' BC ' boundary condition...']);

end

%% Propagation in z direction

for kz = 1:1:size(nc,3)-1
%     kz = 1;

    % Extract known field values for current propagation step

    pb = phi(:,:,kz);

    %% Boundary for TBC

%     if strcmp(BC,'TBC') == 1

        % Define address vector of relevant boundary elements
        kx = 2:1:size(nc,2)-1;
        ky = 2:1:size(nc,1)-1;

        % Initialize field coefficient according to correct direction (NSWE)
        eta_N=zeros(1,size(nc,2)); % eta_N=zeros(1,size(nc,2)-2);
        eta_S=zeros(1,size(nc,2));
        eta_W=zeros(size(nc,1),1);
        eta_E=zeros(size(nc,1),1);

        % Find those values of inner boundary element that is not zero.
        % This avoides the division by zero when calculating eta quotient.
        kx_N=find(pb(2,kx)~=0);
        kx_S=find(pb(end-2,kx)~=0);
        ky_W=find(pb(ky,2)~=0);
        ky_E=find(pb(ky,end-2)~=0);

        % Calculation of the wave numbers 
        eta_N(kx_N)=pb(3,kx_N+1)./pb(2,kx_N+1); % (第三行和第二行的比值 当作 北边透明边界条件传递值)
        eta_S(kx_S)=pb(end-1,kx_S+1)./pb(end-2,kx_S+1);  % (倒数第二行和倒数第三行的比值 当作 南边透明边界条件传递值)
        eta_W(ky_W)=pb(ky_W+1,3)./pb(ky_W+1,2); % (第三列和第二列的比值 当作 西边透明边界条件传递值)
        eta_E(ky_E)=pb(ky_E+1,end-1)./pb(ky_E+1,end-2);% (倒数第二列和倒数第三列的比值 当作 东边透明边界条件传递值)
        
%%%%%%%%  Wave number in x-direction k_x
        kW=zeros(dim_yl,1);
        kW(ky_W)=1/(dx_fine.*1i)*log(eta_W(ky_W));  %%% 为何此处取值dx_fine???
        kE=zeros(dim_yl,1);
        kE(ky_E)=-1/(dx_fine.*1i)*log(eta_E(ky_E));

%%%%%%%%  Wave number in y direction k_y
        kN=zeros(1,dim_xl);
        kN(kx_N)=1/(dy_fine.*1i)*log(eta_N(kx_N));
        kS=zeros(1,dim_xl);
        kS(kx_S)=-1/(dy_fine.*1i)*log(eta_S(kx_S));

        % Change sign of the real part of the wave number if necessary 
        % wave number 的实部必须是负数才能保证无边界反射
        
        W_korrekt=find(real(kW)<0);
        kW_korrekt=kW(W_korrekt);
        Real_kW_korrekt=real(kW_korrekt).*(-1);
        Imag_kW=imag(kW_korrekt).*(i);
        kW(W_korrekt)=Real_kW_korrekt+Imag_kW ;

        E_korrekt=find(real(kE)<0);
        kE_korrekt=kE(E_korrekt);
        Real_kE_korrekt=real(kE_korrekt).*(-1);
        Imag_kE=imag(kE_korrekt).*(i);
        kE(E_korrekt)=Real_kE_korrekt+Imag_kE;

        N_korrekt=find(real(kN)<0);
        kN_korrekt=kN(N_korrekt);
        Real_kN_korrekt=real(kN_korrekt).*(-1);
        Imag_kN=imag(kN_korrekt).*(i);
        kN(N_korrekt)=Real_kN_korrekt+Imag_kN;

        S_korrekt=find(real(kS)<0);
        kS_korrekt=kS(S_korrekt);
        Real_kS_korrekt=real(kS_korrekt).*(-1);
        Imag_kS=imag(kS_korrekt).*(i);
        kS(S_korrekt)=Real_kS_korrekt+Imag_kS;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % phi0=phi1exp(-j*k*delta_x)
        pb(ky,1)  =pb(ky,2).*exp(-1i*dx_fine.*kW);%West       
        pb(ky,end)=pb(ky,end-1).*exp(-1i*dx_fine.*kE);%East
        pb(end,kx)=pb(end-1,kx).*exp(-1i*dy_fine.*kS);%South
        pb(1,kx)  =pb(2,kx).*exp(-1i*dy_fine.*kN);%North

%     end

    %% Call functions for calculation of diagonals

%[ diagBC,diagBN,diagBS,diagBE,diagBW ] = diagonalsPade(beta_0,k_bar,nc(:,:,kz),xg(:,:,kz),yg(:,:,kz),dim_y,dim_xl,dim_yl,grid,gammaBoundaryCondition,POLARIZATION,FX,BC);
%[ diagC,diagN,diagS,diagE,diagW ]      = diagonalsPade(beta_0,k_bar,nc(:,:,kz+1),xg(:,:,kz+1),yg(:,:,kz+1),dim_y,dim_xl,dim_yl,grid,gammaBoundaryCondition,POLARIZATION,FX,BC);
    epsilon = nc(:,:,kz+1).*nc(:,:,kz+1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    globalIndexSlgs   = grid.globalIndex(2:end-1,2:end-1);
    globalAdrSlgs   = reshape(globalIndexSlgs,size(globalIndexSlgs,1)*size(globalIndexSlgs,2),1);

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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        kx = 2:1:size(nc,2)-1;
        ky = 2:1:size(nc,1)-1;
        
        localAdrNTBC = grid.localIndex(2,kx)';
        localAdrETBC = grid.localIndex(ky,end-1)';
        localAdrSTBC = grid.localIndex(end-1,kx)';
        localAdrWTBC = grid.localIndex(ky,2)';

        globalAdrNTBC = grid.globalIndex(2,kx)';
        globalAdrETBC = grid.globalIndex(ky,end-1)';
        globalAdrSTBC = grid.globalIndex(end-1,kx)';
        globalAdrWTBC = grid.globalIndex(ky,2)';

        gammaN = ones(dim_xl*dim_yl,1);
        gammaE = ones(dim_xl*dim_yl,1);
        gammaS = ones(dim_xl*dim_yl,1);
        gammaW = ones(dim_xl*dim_yl,1);
        
        phi_N = find(pb(2,kx)~=0);
        phi_S = find(pb(end-1,kx)~=0);
        phi_W = find(pb(ky,2)~=0);
        phi_E = find(pb(ky,end-1)~=0);
%         gammaN(localAdrNTBC) = gammaBoundaryCondition{1}; 
%         gammaE(localAdrETBC) = gammaBoundaryCondition{2};
%         gammaS(localAdrSTBC) = gammaBoundaryCondition{3};
%         gammaW(localAdrWTBC) = gammaBoundaryCondition{4};

        gamma_N=zeros(dim_xl,1);
        gamma_N(phi_N) =pb(1,phi_N+1)./pb(2,phi_N+1);
        gammaN(localAdrNTBC) = gamma_N;
        gamma_S=zeros(dim_xl,1);
        gamma_S(phi_S) =(pb(end,phi_S+1)./pb(end-1,phi_S+1));
        gammaS(localAdrSTBC) = gamma_S;
        gamma_E=zeros(dim_yl,1);
        gamma_E(phi_E) =(pb(phi_E+1,end)./pb(phi_E+1,end-1));
        gammaE(localAdrETBC) =gamma_E;
        gamma_W=zeros(dim_yl,1);
        gamma_W(phi_W)=(pb(phi_W+1,1)./pb(phi_W+1,2));
        gammaW(localAdrWTBC) =gamma_W;
        
        
        %% Generation of diagonals
% [ diagC,diagN,diagS,diagE,diagW ]      = diagonalsPade(beta_0,neff,nc(:,:,kz+1),xgc(:,:,kz+1),ygc(:,:,kz+1),dim_y,dim_xl,dim_yl,grid,gammaBoundaryCondition,POLARIZATION,FX,BC,pb);
%%       diagNTBC = gammaN .* padeN(n,yg,dim_xl,dim_yl,localAdrNTBC,globalAdrNTBC,POLARIZATION,FIELDCOMPONENT); 
a = yg(:,:,kz+1);
dn = a(globalAdrN - 1) - a(globalAdrN);
ds = a(globalAdrN) - a(globalAdrN + 1);
% Round off error
dn = 1e-12*round(dn*1e12);
ds = 1e-12*round(ds*1e12);
aNorth = (2./(dn.*(dn+ds))) .* ones(length(globalAdrN),1);
aNorthReshaped = zeros(dim_xl*dim_yl,1);
aNorthReshaped(localAdrN) = aNorth;
diagNTBC = gammaN .* aNorthReshaped;



%%       diagSTBC = gammaS .* padeS(nc(:,:,kz+1),ygc(:,:,kz+1),dim_xl,dim_yl,localAdrSTBC,globalAdrSTBC,POLARIZATION,FIELDCOMPONENT);
a = yg(:,:,kz+1);
dn = a(globalAdrS - 1) - a(globalAdrS);
ds = a(globalAdrS) - a(globalAdrS + 1);
% Round off error
dn = 1e-12*round(dn*1e12);
ds = 1e-12*round(ds*1e12);
a_s = (2./(ds.*(ds+dn))) .* ones(length(globalAdrS),1);
aSouthReshaped = zeros(dim_xl*dim_yl,1);
aSouthReshaped(localAdrS) = a_s;
diagSTBC = gammaS .* aSouthReshaped;
        

%%         diagETBC = gammaE .* padeE(nc(:,:,kz+1),ygc(:,:,kz+1),dim_y,dim_xl,dim_yl,localAdrETBC,globalAdrETBC,POLARIZATION,FIELDCOMPONENT);
a = xg(:,:,kz+1);
de = a(globalAdrE + dim_y) - a(globalAdrE);
dw = a(globalAdrE) - a(globalAdrE - dim_y);
% Round off error
de = 1e-12*round(de*1e12);
dw = 1e-12*round(dw*1e12);
epsilon = nc(:,:,kz+1).*nc(:,:,kz+1);

aEast = (2./(de.*(de+dw))) .* (2*epsilon(globalAdrE + dim_y) ./ (epsilon(globalAdrE) + epsilon(globalAdrE + dim_y)));
aEastReshaped = zeros(dim_xl*dim_yl,1);
aEastReshaped(localAdrE) = aEast;
diagETBC = gammaE .* aEastReshaped;
        
        
%%         diagWTBC = gammaW .* padeW(nc(:,:,kz+1),ygc(:,:,kz+1),dim_y,dim_xl,dim_yl,localAdrWTBC,globalAdrWTBC,POLARIZATION,FIELDCOMPONENT); 
a = xg(:,:,kz+1);
de = a(globalAdrW + dim_y) - a(globalAdrW);
dw = a(globalAdrW) - a(globalAdrW - dim_y);
% Round off error
de = 1e-12*round(de*1e12);
dw = 1e-12*round(dw*1e12);
epsilon = nc(:,:,kz+1).*nc(:,:,kz+1);

aWest = 2./(dw.*(dw+de)) .* (2*epsilon(globalAdrW - dim_y) ./ (epsilon(globalAdrW) + epsilon(globalAdrW - dim_y)));
aWestReshaped = zeros(dim_xl*dim_yl,1);
aWestReshaped(localAdrW) = aWest;
diagWTBC = gammaW .* aWestReshaped;
        
%%%%%%%%%%%%
diagTBC = (diagNTBC + diagSTBC + diagETBC + diagWTBC);


%% diagC = padeX(nc(:,:,kz+1),xg(:,:,kz+1),dim_xl,dim_yl,dim_y,localAdrSlgs,globalAdrSlgs,POLARIZATION,{'Ex';'Hy'}) + padeY(nc(:,:,kz+1),yg(:,:,kz+1),dim_xl,dim_yl,localAdrSlgs,globalAdrSlgs,POLARIZATION,{'Ex';'Hy'}) + beta_0^2 .* (epsilon(globalAdrSlgs) - k_bar^2);
a = xg(:,:,kz+1);
de = a(globalAdrSlgs + dim_y) - a(globalAdrSlgs);
dw = a(globalAdrSlgs) - a(globalAdrSlgs - dim_y);
% Round off error
de = 1e-12*round(de*1e12);
dw = 1e-12*round(dw*1e12);

aX = -4./(de.*dw) + aEastReshaped + aWestReshaped;
aY = - aNorthReshaped - aSouthReshaped; 
diagC = aX + aY + beta_0^2 .* (epsilon(globalAdrSlgs) - k_bar^2);

% diagN = padeN(nc(:,:,kz+1),yg(:,:,kz+1),dim_xl,dim_yl,localAdrN,globalAdrN,POLARIZATION,{'Ex';'Hy'});
diagN = aNorthReshaped;

% diagS = padeS(nc(:,:,kz+1),yg(:,:,kz+1),dim_xl,dim_yl,localAdrS,globalAdrS,POLARIZATION,{'Ex';'Hy'});
diagS = aSouthReshaped;
% diagE = padeE(nc(:,:,kz+1),xg(:,:,kz+1),dim_y,dim_xl,dim_yl,localAdrE,globalAdrE,POLARIZATION,{'Ex';'Hy'});
diagE = aEastReshaped;
% diagW = padeW(nc(:,:,kz+1),xg(:,:,kz+1),dim_y,dim_xl,dim_yl,localAdrW,globalAdrW,POLARIZATION,{'Ex';'Hy'});
diagW = aWestReshaped;


diagC = (diagC + diagTBC);
        




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ diagBC,diagBN,diagBS,diagBE,diagBW ] = diagonalsPade(beta_0,neff,nc(:,:,kz),xgc(:,:,kz),ygc(:,:,kz),dim_y,dim_xl,dim_yl,grid,gammaBoundaryCondition,POLARIZATION,FX,BC,pb);
%%         diagNTBC = gammaN .* padeN(nc(:,:,kz),ygc(:,:,kz),dim_xl,dim_yl,localAdrNTBC,globalAdrNTBC,POLARIZATION,FIELDCOMPONENT); 
a = yg(:,:,kz);
dn = a(globalAdrN - 1) - a(globalAdrN);
ds = a(globalAdrN) - a(globalAdrN + 1);
% Round off error
dn = 1e-12*round(dn*1e12);
ds = 1e-12*round(ds*1e12);
aNorth = (2./(dn.*(dn+ds))) .* ones(length(globalAdrN),1);
aNorthReshaped = zeros(dim_xl*dim_yl,1);
aNorthReshaped(localAdrN) = aNorth;
diagNTBC = gammaN .* aNorthReshaped;



%%       diagSTBC = gammaS .* padeS(nc(:,:,kz),ygc(:,:,kz),dim_xl,dim_yl,localAdrSTBC,globalAdrSTBC,POLARIZATION,FIELDCOMPONENT);
a = yg(:,:,kz);
dn = a(globalAdrS - 1) - a(globalAdrS);
ds = a(globalAdrS) - a(globalAdrS + 1);
% Round off error
dn = 1e-12*round(dn*1e12);
ds = 1e-12*round(ds*1e12);
a_s = (2./(ds.*(ds+dn))) .* ones(length(globalAdrS),1);
aSouthReshaped = zeros(dim_xl*dim_yl,1);
aSouthReshaped(localAdrS) = a_s;
diagSTBC = gammaS .* aSouthReshaped;
        

%%         diagETBC = gammaE .* padeE(nc(:,:,kz),ygc(:,:,kz),dim_y,dim_xl,dim_yl,localAdrETBC,globalAdrETBC,POLARIZATION,FIELDCOMPONENT);
a = xg(:,:,kz);
de = a(globalAdrE + dim_y) - a(globalAdrE);
dw = a(globalAdrE) - a(globalAdrE - dim_y);
% Round off error
de = 1e-12*round(de*1e12);
dw = 1e-12*round(dw*1e12);
epsilon = nc(:,:,kz+1).*nc(:,:,kz+1);

aEast = (2./(de.*(de+dw))) .* (2*epsilon(globalAdrE + dim_y) ./ (epsilon(globalAdrE) + epsilon(globalAdrE + dim_y)));
aEastReshaped = zeros(dim_xl*dim_yl,1);
aEastReshaped(localAdrE) = aEast;
diagETBC = gammaE .* aEastReshaped;
        
        
%%         diagWTBC = gammaW .* padeW(nc(:,:,kz),ygc(:,:,kz),dim_y,dim_xl,dim_yl,localAdrWTBC,globalAdrWTBC,POLARIZATION,FIELDCOMPONENT); 
a = xg(:,:,kz);
de = a(globalAdrW + dim_y) - a(globalAdrW);
dw = a(globalAdrW) - a(globalAdrW - dim_y);
% Round off error
de = 1e-12*round(de*1e12);
dw = 1e-12*round(dw*1e12);
epsilon = nc(:,:,kz+1).*nc(:,:,kz+1);

aWest = 2./(dw.*(dw+de)) .* (2*epsilon(globalAdrW - dim_y) ./ (epsilon(globalAdrW) + epsilon(globalAdrW - dim_y)));
aWestReshaped = zeros(dim_xl*dim_yl,1);
aWestReshaped(localAdrW) = aWest;
diagWTBC = gammaW .* aWestReshaped;
        
%%%%%%%%%%%%
diagTBC = (diagNTBC + diagSTBC + diagETBC + diagWTBC);


%% diagC = padeX(nc(:,:,kz),xg(:,:,kz),dim_xl,dim_yl,dim_y,localAdrSlgs,globalAdrSlgs,POLARIZATION,{'Ex';'Hy'}) + padeY(nc(:,:,kz+1),yg(:,:,kz+1),dim_xl,dim_yl,localAdrSlgs,globalAdrSlgs,POLARIZATION,{'Ex';'Hy'}) + beta_0^2 .* (epsilon(globalAdrSlgs) - k_bar^2);
a = xg(:,:,kz);
de = a(globalAdrSlgs + dim_y) - a(globalAdrSlgs);
dw = a(globalAdrSlgs) - a(globalAdrSlgs - dim_y);
% Round off error
de = 1e-12*round(de*1e12);
dw = 1e-12*round(dw*1e12);
aX = -4./(de.*dw) + aEastReshaped + aWestReshaped;
aY = - aNorthReshaped - aSouthReshaped; 
diagBC = aX + aY + beta_0^2 .* (epsilon(globalAdrSlgs) - k_bar^2);

% diagN = padeN(nc(:,:,kz),yg(:,:,kz),dim_xl,dim_yl,localAdrN,globalAdrN,POLARIZATION,{'Ex';'Hy'});
diagBN = aNorthReshaped;

% diagS = padeS(nc(:,:,kz),yg(:,:,kz),dim_xl,dim_yl,localAdrS,globalAdrS,POLARIZATION,{'Ex';'Hy'});
diagBS = aSouthReshaped;
% diagE = padeE(nc(:,:,kz),xg(:,:,kz),dim_y,dim_xl,dim_yl,localAdrE,globalAdrE,POLARIZATION,{'Ex';'Hy'});
diagBE = aEastReshaped;
% diagW = padeW(nc(:,:,kz),xg(:,:,kz),dim_y,dim_xl,dim_yl,localAdrW,globalAdrW,POLARIZATION,{'Ex';'Hy'});
diagBW = aWestReshaped;


diagBC = (diagC + diagTBC);  
    




    %% Apply multistep method

    for ii = 1:1:1  % Padé(1,1) only requires one multistep

        % Merge diagonals in system matrix

        A = sparse((size(nc,1)-2)*(size(nc,2)-2),(size(nc,1)-2)*(size(nc,2)-2));
        A = spdiags(1 + VX(ii)*diagC,0,A);
%         A = spdiags(1 + diagC,0,A);
        A = spdiags([VX(ii)*diagN(2:end); 0],-1,A);
        A = spdiags([0; VX(ii)*diagS(1:end-1)],1,A);
        A = spdiags([zeros(dim_yl,1); VX(ii)*diagE(1:end-dim_yl)],dim_yl,A);
        A = spdiags([VX(ii)*diagW(dim_yl+1:end); zeros(dim_yl,1)],-dim_yl,A);
        
%         spy(A);
        
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
end
spy(A);

%% Visualization
Ex_int = squeeze(sum(abs(phi) .* abs(yg)));
figure
surf(z(2:end),x,squeeze(Ex_int(:,2:end))/max(max(Ex_int(:,2:end))))
shading flat
xlabel('z [\mum]')
ylabel('x [\mum]')
title('Normalized absolute E-field (after y-integration) [a.u.]')
axis([z(2) z(end) x(1) x(end) 0 1])

%% Visualization - slice
slice = 101;
test = phi(:,:,slice);

figure
surf(xg(:,:,1),yg(:,:,1),abs(test))
xlabel('x')
ylabel('y')
title(['field distribution at ', num2str(slice) ,'. slice'])
shading flat

figure
surf(xg(:,:,1),yg(:,:,1),abs(phi(:,:,1)))
xlabel('x')
ylabel('y')
title('Exciting field distribution')
shading flat
%% test

% test1 =  A(1,:);
% test2 =  A(2,:);
% test3 =  A(3,:);
% test4 =  A(4,:);
% test154 =  A(154,:);
% test10000 =  A(10000,:);
% test20000 =  A(20000,:);
% test30000 =  A(30000,:);
% test40000 =  A(40000,:);
% test54001 =  A(54001,:);