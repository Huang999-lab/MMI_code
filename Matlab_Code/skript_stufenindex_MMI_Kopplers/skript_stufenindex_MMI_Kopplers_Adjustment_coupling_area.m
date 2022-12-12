%% Adjustment of the coupling area/Anpassung des Koppelbereichs
dx_coarse = 1e-6;                    % Step coarse for fine step               
dy_coarse = 1e-6;                    % Step coarse for fine step
dz        = 1e-6;                    % Propagationstep µm(1e-6)
L_x = 99e-6; % Länge der querschnitt in x
L_y = 99e-6;
n_glass = 1.522;
n_1 = 1.530;

%% Grid
x = linspace(-L_x/2,L_x/2,100);    
y = linspace(-L_y/2,L_y/2,100);
% x = [0: dx_coarse: 99*dx_coarse];
% y = [0: dy_coarse: 99*dy_coarse];

z = [0:dz:199*dz];                                                                % Base vector
[xg,yg,zg] = meshgrid(x,y,z);                                                       % Base grid

n = ones(length(y),length(x))*n_glass; 
n_dim = size(n);

%% coupling area
% rectangle: length=width=10*dx_coarse

% len_Ein = 10;
% n(n_dim(1)/2 - len_Ein:n_dim(1)/2 + len_Ein, n_dim(2)/2 - len_Ein:n_dim(2)/2 + len_Ein) = 1.530;

% Circle: radius=10*dx_coarse
z_coup = [0:dz:49*dz];

rad_coup = 10*dx_coarse;
a = xg(:,:,1);
center_button_x = a(50,50);
a = yg(:,:,1);
center_button_y = a(20,20);

B=sqrt((xg(:,:,1)-center_button_x).^2+(yg(:,:,1)-center_button_y).^2)<=rad_coup;
n(B)=n_1;
nc = repmat(n,[1 1 length(z_coup)]);

figure
surf(xg(:,:,1),yg(:,:,1),nc(:,:,1))
xlabel('x')
ylabel('y')
zlabel('n')
title('Coupling cross section (slice 1-50)')
%% MMI area
% rectangle: length=width=10*dx_coarse
z_MMI = [50*dz:dz:149*dz];

n = ones(length(y),length(x))*n_glass; 
len_MMI = 80;
n(n_dim(1)/2 - len_MMI/2:n_dim(1)/2 + len_MMI/2, n_dim(2)/2 - len_MMI/2:n_dim(2)/2 + len_MMI/2) = n_1;
n = repmat(n,[1 1 length(z_MMI)]);
nc = cat(3,nc,n);

figure
surf(xg(:,:,1),yg(:,:,1),nc(:,:,11))
xlabel('x')
ylabel('y')
zlabel('n')
title('MMI section (slice 51-150)')

%% decoupling area
% rectangle: length=width=10*dx_coarse

% len_Ein = 10;
% n(n_dim(1)/2 - len_Ein:n_dim(1)/2 + len_Ein, n_dim(2)/2 - len_Ein:n_dim(2)/2 + len_Ein) = 1.530;


% Circle: radius=10*dx_coarse
n = ones(length(y),length(x))*n_glass; 
z_decoup = [150*dz:dz:199*dz];

rad_decoup_top = 10*dx_coarse;
a = xg(:,:,1);
center_top_x = a(50,50);
a = yg(:,:,1);
center_top_y = a(80,80);

rad_decoup_button = 10*dx_coarse;
a = xg(:,:,1);
center_button_x = a(50,50);
a = yg(:,:,1);
center_button_y = a(20,20);

B=sqrt((xg(:,:,1)-center_top_x).^2+(yg(:,:,1)-center_top_y).^2)<=rad_decoup_top;
n(B)=n_1;
B=sqrt((xg(:,:,1)-center_button_x).^2+(yg(:,:,1)-center_button_y).^2)<=rad_decoup_button;
n(B)=n_1;

n = repmat(n,[1 1 length(z_decoup)]);
nc = cat(3,nc,n);

figure
surf(xg(:,:,1),yg(:,:,1),nc(:,:,191))
xlabel('x')
ylabel('y')
zlabel('n')
title('Decoupling cross section (slice 151-199)')

%% x-z Visualization
figure
n_int = squeeze(sum(abs(nc) .* abs(yg)));
surf(z(2:end),x,squeeze(n_int(:,2:end))/max(max(n_int(:,2:end))))

xlabel('z [\mum]')
ylabel('x [\mum]')
zlabel('n')
title('Adjustment coupling: x-z cross section after y-integration')
