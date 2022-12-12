%% Mesh Definition, odd grid
N = 49;%x
M = 49;%y
delta_x = 1;
delta_y = 1;
[x,y] = ndgrid(1:delta_x:N,1:delta_y:M);

%%set boundary condition
u = zeros(N/delta_x,M/delta_y);
% set sigma/eps = 1
sigma = 1;
eps = 1;

% boundary condition
V_0 = 5;
u(1,:) = V_0;
u(:,1) = V_0;
u(end,:) = V_0;
u(:,end) = V_0;
B = 1; % Band
% u((M+1)/2-B,(M+1)/2-B:(M+1)/2+B) = sigma/eps;
% u((M+1)/2+B,(M+1)/2-B:(M+1)/2+B) = sigma/eps;
% u((M+1)/2-B:(M+1)/2+B,(M+1)/2-B) = sigma/eps;
% u((M+1)/2-B:(M+1)/2+B,(M+1)/2+B) = sigma/eps;

u((M+1)/2-B:(M+1)/2+B,(M+1)/2-B:(M+1)/2+B) = sigma/eps;


%% Gauss-Seidel iteration
% delta_u(i,j) = (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - 4*u(i,j))/h^2 = 0;
for iter = 1:100
    % u_1:u_40 forwards G-S
    for j = 2:(M+1)/2-B-1
        for i = 2:N-1
            u(i,j) = (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1))/4;
        end
    end
    
    for j = (M+1)/2-B:(M+1)/2+B
        for i = 2:(M+1)/2-B-1
            u(i,j) = (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1))/4;
        end
        
        for i = (M+1)/2+B+1:N-1
            u(i,j) = (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1))/4;
        end
    end
    
    for j = (M+1)/2+B+1:N-1
        for i = 2:N-1
            u(i,j) = (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1))/4;
        end
    end
    
    
    % u_1:u_40 backwards G-S
    for j = (M+1)/2-B-1:-1:2
        for i = N-1:-1:2
            u(i,j) = (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1))/4;
        end
    end
    
    for j = (M+1)/2+B:-1:(M+1)/2-B
        for i = (M+1)/2-B-1:-1:2
            u(i,j) = (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1))/4;
        end
        
        for i = N-1:-1:(M+1)/2+B+1
            u(i,j) = (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1))/4;
        end
    end
    
    for j = N-1:-1:(M+1)/2+B+1
        for i = N-1:-1:2
            u(i,j) = (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1))/4;
        end
    end
end

figure
mesh(x,y,u)
xlabel('x')
ylabel('y')
zlabel('u')
title('Gauss-Seidel iteration solotion')

%% Matrix form (M-2)*(N-2) x (M-2)*(N-2) Matrix
diag_block = eye(M-2)*(-2*(1/delta_x^2 + 1/delta_y^2));
diag_block = diag_block + diag(ones(M-3,1)/delta_y^2,1);
diag_block = diag_block + diag(ones(M-3,1)/delta_y^2,-1);
Matrix_1 = kron(eye(N-2), diag_block); % assum no Bound inside --> (M-2)*(N-2) x (M-2)*(N-2) Matrix

Matrix_1 = Matrix_1 + diag(ones((M-3)*(N-2),1) / delta_x^2, N-2);
Matrix_1 = Matrix_1 + diag(ones((M-3)*(N-2),1) / delta_x^2, -(N-2));
spy(Matrix_1)


f = zeros(size(Matrix_1,1),1);% (M-2)*(N-2) vector
% boundary condition

% -1/h^2*10
fa_index = [1;N-2;(N-2)*(M-3)+1;size(Matrix_1,1)];%four corners
f(fa_index) = -1/delta_x^2*10;


% -1/h^2*sigma/eps
innen_index = find(u~=5);
box_index = find(u(innen_index) == sigma/eps);% innen_index_index
f(box_index) = 2/delta_x^2*sigma/eps; %% not sure !!!!!!!!!!!!!!!!     -2/delta_x^2*sigma/eps

f1_index = box_index+1;
f1_index = setdiff(f1_index,box_index);%box bottom

f2_index = box_index-1;
f2_index = setdiff(f2_index,box_index);%box top

f3_index = box_index-(N-2);
f3_index = setdiff(f3_index,box_index);%box left

f4_index = box_index+N-2;
f4_index = setdiff(f4_index,box_index);%box right
fb_index = [f1_index;f2_index; f3_index; f4_index];
f(fb_index) = -1/delta_x^2*sigma/eps;




% -1/h^2*5
f6_index = [2:N-3]';%bound left
f7_index = [length(innen_index)-(N-4):length(innen_index)-1]';%bound right
f8_index = [N-1:N-2:length(innen_index)-2*(N-2)+1]';%bound top
f9_index = [2*(N-2):N-2:length(innen_index)-(N-2)]';%bound bottom
fc_index = [f6_index;f7_index; f8_index; f9_index];
f(fc_index) = -1/delta_x^2*5;

% 0
test_index = [fa_index; box_index ;fb_index;fc_index];
fd_index = setdiff(1:length(innen_index),test_index);%0 condition
f(fd_index) = 0;

% solution
sol_2 = Matrix_1\f;
sol_2 = reshape(sol_2,(N-2),(M-2));
sol_2 = [sol_2,ones((N-2),1)*V_0];
sol_2 = [ones((N-2),1)*V_0,sol_2];
sol_2 = [sol_2;ones(N,1)'*V_0];
sol_2 = [ones(N,1)'*V_0;sol_2];

figure
mesh(x,y,sol_2)
xlabel('x')
ylabel('y')
zlabel('u')
title('Matrix solotion')
%% Matrix form --- 40 x 40 Matrix

% v_index = find(~u); % points index Linear System
% box_index = find(u==sigma/eps);
% 
% % -1/h^2*sigma/eps
% f1_index = box_index+1;
% f1_index = setdiff(f1_index,box_index);%box bottom
% 
% f2_index = box_index-1;
% f2_index = setdiff(f2_index,box_index);%box top
% 
% f3_index = box_index-N;
% f3_index = setdiff(f3_index,box_index);%box left
% 
% f4_index = box_index+N;
% f4_index = setdiff(f4_index,box_index);%box right
% fa_index = [f1_index;f2_index; f3_index; f4_index];
% 
% % -1/h^2*10
% fb_index = [v_index(1);v_index(N-2);v_index(end-(N-2)+1); v_index(end)];%four corners
% 
% 
% % -1/h^2*5
% f6_index = [v_index(2):v_index(N-3)]';%bound left
% f7_index = [v_index(end-(N-4)):v_index(end-1)]';%bound right
% f8_index = [v_index(N-1):N:v_index(end-2*(N-2)+1)]';%bound top
% f9_index = [v_index(2*(N-2)):N:v_index(end-N+3)]';%bound bottom
% fc_index = [f6_index;f7_index; f8_index; f9_index];
% 
% % 0
% test_index = [fa_index;fb_index; fc_index];
% fd_index = setdiff(v_index,test_index);%0 condition


%% Matrix form
% point = M*N -(2*M+2*N-4) - (2*B+1)^2;
% A = zeros(point);
% 
% for i = 1:point
%     A(i,i) = -2*(1/h^2 + 1/h^2); % A(i,i) = -4
% end
% 
% 
% 
% 
% % first 14 points
% for i = 1:N-3 % 1:6 --> point i at j Column vector
%     for j = 1:(M+1)/2-B-2  % 1:(M+1)/2-B-2 --> 1:2
%         A(i+(j-1)*(N-2), i+(j-1)*(N-2)+1) = 1/h^2 ; % the bottom coordinate of the selected coordinate (point to multiply by u)
%         A(i+(j-1)*(N-2)+1, i+(j-1)*(N-2)) = 1/h^2 ; % the top coordinate of the selected coordinate (point to multiply by u)
%     end
% end
% 
% for i = 1:N-3 % 1:6 --> point i at j Column vector
%     for j = 1:(M+1)/2-B-2  % 1:(M+1)/2-B-2 --> 1:2
%         A(i+(j-1)*(N-2), i+j*(N-2)) = 1/h^2 ; % the right coordinate of the selected coordinate (point to multiply by u)
%         A(i+j*(N-2), i+(j-1)*(N-2)) = 1/h^2 ; % the left coordinate of the selected coordinate (point to multiply by u)
%     end
% end

% points 15-18: