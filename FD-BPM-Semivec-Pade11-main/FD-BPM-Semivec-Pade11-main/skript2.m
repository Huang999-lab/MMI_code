% for kz = 1:1:size(nc,3)-1  % 一层一层迭代
kz = 1;
    % Extract known field values for current propagation step

    pb = phi(:,:,kz);

    %%%% Call functions for calculation of diagonals :only for TE,Ex, "BC" = ABC
%%   [ diagC,diagN,diagS,diagE,diagW ]      = diagonalsPade(beta_0,neff,nc(:,:,kz+1),xg(:,:,kz+1),yg(:,:,kz+1),dim_y,dim_xl,dim_yl,grid,gammaBoundaryCondition,POLARIZATION,FX,BC);

%%%%diagC = padeX(n,xg,dim_xl,dim_yl,dim_y,localAdrSlgs,globalAdrSlgs,POLARIZATION,FIELDCOMPONENT) + padeY(n,yg,dim_xl,dim_yl,localAdrSlgs,globalAdrSlgs,POLARIZATION,FIELDCOMPONENT) + beta_0^2 .* (epsilon(globalAdrSlgs) - neff^2);
    a = xg(:,:,kz+1);
    de = a(globalAdrSlgs + dim_y) - a(globalAdrSlgs);% delta e in east direction 右东
    dw = a(globalAdrSlgs) - a(globalAdrSlgs - dim_y);% delta w in west direction 左西
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
    dn = a(globalAdrSlgs - 1) - a(globalAdrSlgs);% delta e in north direction 上北
    ds = a(globalAdrSlgs) - a(globalAdrSlgs + 1);% delta e in south direction 下南
    
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
    
    
    
%%        diagN = padeN(beta_0,yg(:,:,kz+1),dim_xl,dim_yl,localAdrN,globalAdrN,POLARIZATION,FIELDCOMPONENT);
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

%%      diagS = padeS(beta_0,yg(:,:,kz+1),dim_xl,dim_yl,localAdrS,globalAdrS,POLARIZATION,FIELDCOMPONENT);% aSorthReshaped
    a = yg(:,:,kz+1);
    dn = a(globalAdrN - 1) - a(globalAdrN);
    ds = a(globalAdrN) - a(globalAdrN + 1);
    a_s = (2./(ds.*(ds+dn))) .* ones(length(globalAdrS),1);
    aSouthReshaped = zeros(dim_xl*dim_yl,1);
    aSouthReshaped(localAdrS) = a_s;
    diagS = aSouthReshaped;
    

%%    diagE = padeE(beta_0,xg(:,:,kz+1),dim_y,dim_xl,dim_yl,localAdrE,globalAdrE,POLARIZATION,FIELDCOMPONENT);% aEastReshaped
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
    
%%     diagW = padeW(beta_0,xg(:,:,kz+1),dim_y,dim_xl,dim_yl,localAdrW,globalAdrW,POLARIZATION,FIELDCOMPONENT);% aWestReshaped
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

%     % Update progress bar
% 
%     if strcmp(PROGRESS,'bar') == 1
% 
%         stepPercent = 1;  % Defines step size for progress bar update
% 
%         if floor(100*kz/(size(nc,3)-1)) > c*stepPercent
% 
%           if ishandle(h)    % Check if waitbar has not been closed to prevent error
% 
%             minutesPassed      = floor(toc/60);
%             secondsPassed      = toc - minutesPassed*60;
%             minutesRemaining   = floor((toc/c)*((100/stepPercent)-c)/60);
%             secondsRemaining   = (toc/c)*((100/stepPercent)-c) - minutesRemaining*60;
% 
%             waitbar(kz/(size(nc,3)-1),h,['Progress: ' num2str(c*stepPercent) '%. Time remaining: ' num2str(floor(minutesRemaining)) 'm ' num2str(ceil(secondsRemaining)) 's.'])
% 
%             c = c + 1;
% 
%           else
% 
%              disp('Process interrupted')
%              
% 
%           end
%         end
% 
%     elseif strcmp(PROGRESS,'cl') == 1
% 
%         stepPercent = 10;  % Defines step size for progress bar update
% 
%         if floor(100*kz/(size(nc,3)-1)) > c*stepPercent
% 
%             minutesPassed      = floor(toc/60);
%             secondsPassed      = toc - minutesPassed*60;
%             minutesRemaining   = floor((toc/c)*((100/stepPercent)-c)/60);
%             secondsRemaining   = (toc/c)*((100/stepPercent)-c) - minutesRemaining*60;
% 
%             out = ['   Progress: ' num2str(c*stepPercent) '%. Time remaining: ' num2str(floor(minutesRemaining)) 'm ' num2str(ceil(secondsRemaining)) 's.'];
%             disp(out)
% 
%             c = c + 1;
% 
%         end

%     end
% end
