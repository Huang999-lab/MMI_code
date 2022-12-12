%% Exemplary BPM calculation for a rectangular step-index MMI-element
%% The external material at the beginning of the input and output channels is glass
%% Parameter: L_x L_y L_coup L_MMI
clear
close all
clc

% Create the concentration profile of an open MMI-Element
% Parameter definition
PROGRESS = 'cl';

lambda              = 1330e-9;       % Wavelength
alpha               = 0.5;           % Crank-Nicolson-Constant
delta               = 0; % or some small value 1e-6
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
dz        = .5e-6;                    % Propagationstep µm(1e-6)

Excitation.type = 'gauss';           % Type of Excitation
Excitation.visualize = 0;            % Plot Excitation field
Excitation.sigma_x = 3e-6;           % Only needed for fieldtype 'gauss'
Excitation.sigma_y = 3e-6;           % Only needed for fieldtype 'gauss'
Excitation.threshold = 1-1/exp(1);   % Only needed for fieldtype 'full'
Excitation.numberInterpolation = 0;  % Only needed for fieldtype 'modal'
Excitation.mode = 'k_bar';           % Only needed for fieldtype 'modal'

% Some possible manipulations of index profile
append = 0;
append_dist = 0e-6;     % um
prepend = 0;
prepend_dist = 0e-6;    % um
upper_cladding = 0;

%% Generation of the refractive index profile/ symmetrical coupling area/Symmetrischer Koppelbereich
% coupling area 1000µm
% querschnitt: y*x = 50µm*100µm
out = 'Loading index profile...';
disp(out)
tic

beta_0 = 2*pi/lambda;       % Wave number
% Loading index profile and base vectors
L_x = 150e-6; % Länge der querschnitt in x
L_y = 100e-6;
n_glass = 1.522;
n_1 = 1.530;
neff = (n_glass+n_1)/2;

% Grid
x = linspace(-L_x/2,L_x/2,151);    
y = linspace(-L_y/2,L_y/2,101);
% Removing round-off errors
x = 1e-12*round(x*1e12);
y = 1e-12*round(y*1e12);
    
z = [0:dz:2000e-6];   %  Länge der gesamte z                                                              % Base vector
[xg,yg,zg] = meshgrid(x,y,z);                                                       % Base grid

xgb = squeeze(xg(:,:,1));  % Transversal grid of first slice
ygb = squeeze(yg(:,:,1));  % Transversal grid of first slice

%% coupling Input area
distance_y_center = 0;
distance_x_center = -50;
Excitation.coup_center = [y((length(y)+1)/2 + distance_y_center),x((length(x)+1)/2 + distance_x_center)]; % (y,x)

L_coup = 200e-6;
z_coup = [0:dz:L_coup];
rad_coup = 25e-6;
Excitation.gausscoup = sqrt((xg(:,:,1)-Excitation.coup_center(2)).^2+(yg(:,:,1)-Excitation.coup_center(1)).^2)<=rad_coup;

n = ones(length(y),length(x))*n_glass;
n(Excitation.gausscoup) = n_1;
nc = repmat(n,[1 1 length(z_coup)]);

n_dim = size(n);
%% MMI area
% rectangle: length=width
L_MMI = 1800e-6;
z_MMI = [L_coup+dz:dz:L_coup+L_MMI];

x_width_MMI = 140;
y_width_MMI = 90;

n = ones(length(y),length(x))*n_glass; 
n(n_dim(1)/2 - y_width_MMI/2:n_dim(1)/2 + y_width_MMI/2, n_dim(2)/2 - x_width_MMI/2:n_dim(2)/2 + x_width_MMI/2) = n_1;
a = repmat(n,[1 1 length(z_MMI)]);
nc = cat(3,nc,a);

%% decoupling Output area
% L_decoup = 200e-6;
% z_decoup = [z_MMI(end)+dz:dz:z(end)];
% rad_decoup = 5e-6;
% 
% distance_x_center = 8;
% distance_y_center = 8;
% top_center = [y((length(y)+1)/2 + distance_y_center),x((length(x)+1)/2)]; % (y,x)
% butten_center = [y((length(y)+1)/2 - distance_y_center),x((length(x)+1)/2)]; % (y,x)
% 
% n = ones(length(y),length(x))*n_luft; 
% B1=sqrt((xg(:,:,1)-top_center(2)).^2+(yg(:,:,1)-top_center(1)).^2)<=rad_decoup;
% B2=sqrt((xg(:,:,1)-butten_center(2)).^2+(yg(:,:,1)-butten_center(1)).^2)<=rad_decoup;
% n(B1) = n_1;
% n(B2) = n_1;
% 
% a = repmat(n,[1 1 length(z_decoup)]);
% nc = cat(3,nc,a);


out = ['  Index profile loaded: ' num2str(toc) 's.'];
disp(out)
tic
% figure
% surf(xg(:,:,1),yg(:,:,1),nc(:,:,1))
% xlabel('x')
% ylabel('y')
% zlabel('n')
% title('MMI section slice.1')
% % x-z Visualization
% figure
% n_int = squeeze(sum(abs(nc) .* abs(yg)));
% surf(z(2:end),x,squeeze(n_int(:,2:end))/max(max(n_int(:,2:end))))
% 
% xlabel('z [\mum]')
% ylabel('x [\mum]')
% zlabel('n')
% title('symmetrical: x-z cross section after y-integration')
%% Append exit
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

%% Prepend input
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

%% Attaching upper cladding
    
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

    nc = cat(1,nc,nCladd);

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

grid.localIndex(2:end-1,2:end-1) = reshape(1:1:dim_xl*dim_yl',dim_yl,dim_xl);
grid.globalIndex(1:end) = 1:1:length(grid.globalIndex(1:end));

out = ['  Defined necessary parameters: ' num2str(toc) 's.'];
disp(out)

%% Execute BPM

out = [POLARIZATION '-BPM with ' BC ' boundary condition started for graded-index multimode diffused waveguide at 1330nm wavelength.'];
disp(out);

[Ex,~,~,~] = FDBPMPade11Semivec(nc,lambda,neff,alpha,solver_tolerance,xg,yg,dz,Excitation,'TE',{'Ex';'Hy'},BC,ABSORBER,PROGRESS);



%% Visualization after y Intergration
Ex_int = squeeze(sum(abs(Ex) .* abs(yg)));
figure
surf(z(2:end)*1e6,x*1e6,squeeze(Ex_int(:,2:end))/max(max(Ex_int(:,2:end))))
colorbar;
shading flat
xlabel('z [\mum]')
ylabel('x [\mum]')
title('Normalized absolute E-field (after y-integration) [a.u.]')
% axis([z(2) z(end) x(1) x(end) 0 1])

%% Visualization at maximums value
Ex_max = squeeze(abs(max(max(Ex))));

figure
plot(z*1e6,Ex_max,'o')
xlabel('z [\mum]')
ylabel('Ex [a.u.]')
title('Maximal E-field at each slice')
%% Visualization at slice
slice = 1;
figure
surf(xg(:,:,1)*1e6,yg(:,:,1)*1e6,abs(Ex(:,:,slice)))
xlabel('x [\mum]')
ylabel('y [\mum]')
title(['field distribution at ',num2str(slice),' slice'])
shading flat %每个网格线段和面具有恒定颜色，该颜色由该线段的端点或该面的角边处具有最小索引的颜色值确定


%% animation
% for i = 1:20:length(z)
%     slicestep=int2str(i);
%     imagesc(abs(Ex(:,:,i)));
% %     caxis([0 0.6]); % 设置颜色图范围
%     colorbar;
%     axis xy;
%     title(['E_x, slice = ',slicestep]);
%     xlabel('x coordinate');
%     ylabel('y coordinate');
%     
%     pause(0.1)
% end