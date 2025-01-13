addpath(genpath('./'));                     % add the whole directory to path, if not already done

%% SET PARAMETERS
c0 = 1;                                     % speed of light m/s (normalized to 1)
lambda0 = 2;                                % central wavelength (um)

skip = 4;                                   % number of iteration frames between plots (higher->faster, lower->more plots)
display_plots = false;                       % plotting during the run?


alpha = 5e2;                                % step size in permittivity (~1e2-1e4 works well)
a = 3;                                     % smooth-max weight factor (see paper)
beta = 0.5;                                 % ratio of electron speed to speed of light
N = 700;                                   % number of iterations

in_material = false;                        % evaluate E_max in material? or in surrounding regions. (NOTE: it doesn't work well, I would suggest just evaluating in optimization region)
starting = 0;                               % 0 -> vacuum, 1 -> random, 2 -> midway epsilon

grids_in_lam = 50;                         % number of grid points in a free space wavelength
% gap_nm_values = 200:10:1300;                % gap size in nm variations with step of 10
% gap_nm_values = 200:100:1300;
gap_nm_values = [800];
L = 1.0;                                   % size of optimization region (um)
% NOTE: if this ^ is too big and the epsilon is too large, the simulations
% can diverge.  This is because there are many degrees of freedom and
% resonance can occur very strongly. Need to try different values and see
% what works.
npml = 10;                                  % number of PML (absorbing region) points (need > 10 at least)

% relative permittivity of material region.  uncomment to select
eps = 3.4363^2;     % Si 2um
%eps = 1.4381^2;      % fused silica 2um
% eps = 1.9834^2;     % Si3N4
%eps = 1.9^2;        % GaOx

nmax = sqrt(eps);    % refractive index of material region

gamma = 0.9;                             % 'momentum term', see paper.  Set between 0-1, can speed up simulation in some cases

%% SET OTHER CONSTANTS (DON'T CHANGE)
dlx = lambda0/grids_in_lam;                 % grid size along electron trajectory axis
dly  = dlx;                                 % spacing in the perpendicular direction

G_best_values = [];                         % Array to store G_best for each gap size
G_best_times_gap_values = [];               % Array to store G_best * gap for each gap size

for gap_nm = gap_nm_values
    pos_src = floor(npml+grids_in_lam/4);       % number of grid points between left edge and source
    spc_pts = floor(grids_in_lam/4);            % number of grid points between source and structure
    gap_pts = floor(gap_nm/1000/dlx);           % number of grid points in the gap
    Lpts = round(L/dlx);                        % number of points in the optimization region
    
    Nx = ceil(lambda0*beta/dlx);                % number of grid points in x
    Ny = gap_pts+2*(pos_src + Lpts + spc_pts);  % number of grid points perpendicular to trajectory
    
    nx = floor(Nx/2);
    ny = floor(Ny/2);
    
    % First compute G maximization, then do G/E_max maximization (for comparison)
    for min_G_Emax = 0
        
        %% This section defines the input parameters that my FDFD code needs to run.
        %  see the FDFD.m code or FDFD_TFSF.m for a more detailed explanation.
        
        ER  = ones(Nx,Ny);                      % relative permittivity grid map
        MuR = ones(Nx,Ny);                      % relative permeability grid map
        ER_best = ones(Nx,Ny);                  % storing the best permittivity map
        A_best = 0;
        
        b = zeros(Nx,Ny);                       % TFSF map.  read up on total-field scattered-field if you are interested.
        b(:, pos_src:pos_src + spc_pts + Lpts + gap_pts + Lpts + spc_pts) = 1;  % define the total field region on the grid
        kinc = [0,1];                           % plane wave incident direction (perp. to electron)
        
        RES = [dlx,dly];                        % grid resolution vector
        BC = [-1,-1];                           % boundary condition vector (periodic if -1)
        NPML = [0,0,npml,npml];                 % PML cells on the boundaries (x-,x+,y-,y+)
        Pol= 'Hz';                              % field polarization
        spc = spc_pts*dly;                      % space between source and objects in um
        gap = gap_pts*dly;                      % gap size in um
        
        xs = dlx*(1:Nx);                        % constant to compute the eta object.  x-pos along gap.
        
        delta_device = zeros(Nx,Ny);            % delta_device is 0 where the permittivity doesn't change.  otherwise it is 1 in the optimization region.
        delta_device(1:Nx, pos_src + spc_pts : pos_src + spc_pts + Lpts) = 1;
        delta_device(1:Nx, pos_src + spc_pts + Lpts + gap_pts : pos_src + spc_pts + Lpts + gap_pts + Lpts) = 1;
        delta_device_vec = delta_device(:);     % vector version of delta_device (matlab likes this better)
        
        % define the eta vector field.  see the paper for more details.
        eta = zeros(Nx,Ny);
        eta(:,ny) = 1/Nx*exp(2*pi*1i*dlx*(0:Nx-1)/lambda0/beta);
        eta_vec = eta(:);
        
        % define stating permittivity based on what value the 'starting variable'
        % holds
        for i = (1:Nx)
            for j = (1:Ny)
                if (delta_device(i,j) == 1)
                    if (starting == 1)
                        ER(i,j) = rand*(eps-1)+1;
                    elseif (starting == 2)
                        ER(i,j) = eps/2+0.5;
                    else
                    end
                end
            end
        end
        
        % run the simulation with accelerator input (plane wave) but all empty space
        [fields, ~] = FDFD_TFSF(ones(Nx,Ny),MuR,RES,NPML,BC,lambda0,Pol,b,kinc);
        
        % get the fields and the E0 (normalization)
        Ex = fields.Ex;
        Ey = fields.Ey;
        E0 = sqrt(abs(Ex(nx, ny))^2 + abs(Ey(nx, ny))^2);
        
        % define variables to store the iteration progress
        G_best = 0;                 % best gradient
        Gs = zeros(N,1);            % gradients over iteration
        E_maxs = zeros(N,1);        % max E-fields over iteration
        G_by_Es = zeros(N,1);       % G/E over iteration, computed directly
        G_by_Sa = zeros(N,1);       % G/E over iteration, computed with smooth-max
        
        phis = zeros(N,1);          % phase of the maximum accelerating input plane wave
        phi = 0;                    % assume input light phase of 0 to start
        AVM_prev = zeros(Nx,Ny);    % store previous sensitivity information for momentum update
        
        if display_plots
            figure(1);                  % open a figure to plot
        end
        
        if ~min_G_Emax
            display('working on gradient maximized structure');
        else
            display('working on acceleration factor maximized structure');
        end
        upd = textprogressbar(N);
        
        for j = (1:N)
            
            upd(j);
            % original simulation (structure in accelerator mode)
            [fields, extra] = FDFD_TFSF(ER,MuR,RES,NPML,BC,lambda0,Pol,b,kinc);
            % get fields
            Ex = fields.Ex/E0;
            Ey = fields.Ey/E0;
            
            % compute gradient (see paper)
            g = sum(sum(eta.*Ex));
            G = real(g);
            
            % get phase
            phis(j) = angle(g);
            
            % get numerical spatial derivative operators
            DEY = extra.derivatives.DEY;
            DEX = extra.derivatives.DEX;
            
            % turn the permittivity map into a vecor.  Then compute some
            % quantities for later.
            ER_vec = ER(:);
            delta_ER = ER > eps/2;
            delta_ER_vec = delta_ER(:);
            chi = delta_device.*(ER - ones(Nx,Ny));
            
            % compute operators from maxwell's eqs. turning Mz into Jx and Jy (Hz into Ex and Ey)
            Ox = -1i*lambda0/2/pi/c0*spdiags(1./ER_vec,0,Nx*Ny,Nx*Ny)*DEY;
            Oy =  1i*lambda0/2/pi/c0*spdiags(1./ER_vec,0,Nx*Ny,Nx*Ny)*DEX;
            
            % create the adjoint vector corresponding to eta (again, see paper)
            eta_aj = [eta_vec; zeros(Nx*Ny,1)];
            
            % if you are evaluating E_max in the material, compute |E| there,
            % otherwise, compute |E| in the full optimization region.
            if (in_material)
                E_abs = (chi/(eps-1)).*sqrt(abs(Ex).^2 + abs(Ey).^2);
            else
                E_abs = delta_device.*sqrt(abs(Ex).^2 + abs(Ey).^2);
            end
            
            % compute auxiliary vectors for later.  (too complicated to explain
            % here.  ask me in person if you're interested).
            x_abs = E_abs(:);
            alpha_vec = exp(x_abs*a);
            alpha_T_1 = sum(alpha_vec);
            Sa = sum(alpha_vec.*x_abs)/alpha_T_1;
            
            X_vec = conj(Ex(:))./x_abs;
            Y_vec = conj(Ey(:))./x_abs;
            
            
            x = [Ex(:); Ey(:)];
            P = [speye(Nx*Ny) speye(Nx*Ny)];
            z = conj(x./[x_abs;x_abs]);
            z(isnan(z)) = 0;
            z(isinf(z)) = 0;
            spdiagz = spdiags(z,0,Nx*Ny*2,Nx*Ny*2);
            
            R = (P*spdiagz);
            %R(isnan(R)) =  0;
            %R = [diag(z(1:Nx*Ny)) diag(z(Nx*Ny+1:end))];
            S = real(1/alpha_T_1*(speye(Nx*Ny) + a*spdiags(x_abs,0,Nx*Ny,Nx*Ny) - a*sum(alpha_vec.*x_abs)/alpha_T_1*speye(Nx*Ny)));
            sigma = transpose(alpha_vec)*S*R;
            sigma(isnan(sigma)) = 0;
            b_aj1 = transpose(G/Sa^2*sigma);
            b_aj2 = -eta_aj/Sa;
            
            % construct final adjoint source
            if (min_G_Emax)
                b_aj = b_aj1 + b_aj2;
            else
                b_aj = -eta_aj;
            end
            b_aj = reshape(Ox*b_aj(1:Nx*Ny) + Oy*b_aj(Nx*Ny+1:end),[Nx,Ny]);
            b_aj(isnan(b_aj)) = 0 ;
            
            % get the factored form of the system operator.
            AF = extra.AF;
            
            % run simulation with adjoint source.
            [fields_aj, ~] = FDFD_fast(ER,MuR,RES,NPML,BC,lambda0,Pol,b_aj,AF);
            
            % get fields
            x_aj = fields_aj.x/E0;
            Ex_aj = reshape(x_aj(1:Nx*Ny),[Nx,Ny]);
            Ey_aj = reshape(x_aj(Nx*Ny+1:end),[Nx,Ny]);
            
            % compute sensitivity information
            AVM = -real((Ex.*Ex_aj.*delta_device + Ey.*Ey_aj.*delta_device));
            
            % record relevant variables in the arrays
            E_max = max(max((E_abs)));
            E_maxs(j) = E_max;
            Gs(j) = G;
            G_by_Es(j) = G/E_max;
            G_by_Sa(j) = G/Sa;
            
            % update permittivity
            ER = ER + alpha*AVM + alpha*gamma*AVM_prev;
            
            % update the previous sensitivity map
            AVM_prev = AVM;
            
            % if permittivity out of bounds, reset inside the correct bounds.
            ER(ER < 1) = 1;
            ER(ER > eps) = eps;
            
            % record best permittivity if applicable
            if (G > G_best)
                G_best = G;
                ER_best = ER;
            end
            
            % plot stuff without too much hastle, display % done
            if (display_plots && mod(j,skip) == 0)
                %    perc_done = j/N*100
                clf;
                subplot(2,2,1);
                disp = [];
                for k = (1:5)
                    disp = [disp; real(ER)];
                end
                
                imagesc(disp,[1,eps])
                colormap(flipud(gray))
                title('relative permittivity')
                set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal')
                set(gca,'FontSize',22,'fontWeight','normal')
                colorbar()
                
                subplot(2,2,2);
                colorbar()
                plot(Gs(1:j),'k');
                xlabel('iteration number')
                ylabel('gradient (E_0)')
                title('acceleration gradient at \phi = 0')
                set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal')
                set(gca,'FontSize',22,'fontWeight','normal')
                
                subplot(2,2,3);
                plot((1:j),G_by_Es(1:j));
                hold all;
                plot((1:j),G_by_Sa(1:j));
                xlabel('iteration number')
                ylabel('G/|E|max')
                title('acceleration factor')
                legend({'actual','smooth-max'})
                set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal')
                set(gca,'FontSize',22,'fontWeight','normal')
                
                
                subplot(2,2,4); hold all;
                plot((1:j),phis(1:j));
                plot((1:j),zeros(j,1));
                xlabel('iteration number');
                ylabel('\phi');
                legend({'computed','\phi=0 (target)'})
                title('acceleration phase (\phi)')
                set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal')
                set(gca,'FontSize',22,'fontWeight','normal')
                
                pause(0.001);
            end
        end
        
        
        %% POST PROCESSING STUFF
        
        % create final field display
        field_disp = [];
        ER_disp = [];
        for i = (1:5)
            field_disp = [field_disp;Ex];
            ER_disp = [ER_disp;  (ER-ones(Nx,Ny))*10000];
        end
        
        % plot movie
        NT = 0;         % number of time steps
        %figure(2);
        for t = (1:NT)
            clf;
            colormap(redbluecmap)
            imagesc(transpose(ER_disp + real(field_disp*exp(-1i*t/40))),[-5,5]); pause(0.0001);
        end
        
        % force the permittivity distribution binary for ER and ER_best
        eps_avg = (eps+1)/2;
        ER(ER<eps_avg) = 1;
        ER(ER>=eps_avg) = eps;
        ER_best(ER_best<eps_avg) = 1;
        ER_best(ER_best>=eps_avg) = eps;
        
        % do another simulation of the binary distribution for ER
        % [fields, extra] = FDFD_TFSF(ER,MuR,RES,NPML,BC,lambda0,Pol,b,kinc);
        % Ex = fields.Ex/E0;
        % Ey = fields.Ey/E0;
        
        % % compute the gradient
        % g_ER = sum(sum(eta.*Ex));
        % G_ER = abs(g_ER); % why abs?
        
        % do another simulation of the binary distribution for ER_best
        [fields_best, extra_best] = FDFD_TFSF(ER_best,MuR,RES,NPML,BC,lambda0,Pol,b,kinc);
        Ex_best = fields_best.Ex/E0;
        Ey_best = fields_best.Ey/E0;
        
        % compute the gradient for ER_best
        g_best = sum(sum(eta.*Ex_best));
        G_best = abs(g_best); % why abs?
        
        % save both gradients into the same file
        timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
        fname = sprintf('result/final_acceleration_gradients_gap_%d_%s_gap_%d.txt', gap_nm, timestamp, gap_nm);
        fileID = fopen(fname, 'w');
        % fprintf(fileID, 'Gradient for ER: %f\n', G_ER);
        fprintf(fileID, 'G_best: %f\n', G_best);
        fprintf(fileID, 'g_best: %f + %fi\n', real(g_best), imag(g_best));
        fclose(fileID);
        fprintf('File saved as: %s\n', fname);
        
        % store G_best value for this gap size
        G_best_values = [G_best_values; G_best];
        G_best_times_gap_values = [G_best_times_gap_values; G_best * gap_nm];
        
        % display and save final structure
        % finalFig = figure('Name','Final Structure','Visible','on');
        
        % binary distributionにした後のERで，繰り返し連結用の変数を初期化
        % disp_final = [];
        % for k = 1:5
        %     % ER を縦方向に 5 回連結
        %     disp_final = [disp_final; real(ER)];
        % end
        
        % % 繰り返した配列を可視化
        % imagesc(disp_final, [1, eps]);  % 2値化後なので [1, eps] の範囲
        % colormap(flipud(gray));
        % axis equal tight;
        % title('Final Binarized Structure');
        % colorbar();
        
        % % タイムスタンプ入りの画像ファイル名 (PNG 等)
        % figName = sprintf('result/final_structure_%s_gap_%d.png', timestamp, gap_nm);
        % saveas(finalFig, figName);  % 画像保存
        
        % display and save best structure
        if display_plots
            bestFig = figure('Name','Best Structure','Visible','on');
        else
            bestFig = figure('Name','Best Structure','Visible','off');
        end
        
        % 繰り返し連結用の変数を初期化
        disp_best = [];
        for k = 1:5
            % ER_best を縦方向に 5 回連結
            disp_best = [disp_best; real(ER_best)];
        end
        
        % 繰り返した配列を可視化
        imagesc(disp_best, [1, eps]);
        colormap(flipud(gray));
        axis equal tight;
        title(sprintf('Best Structure (gap = %d nm)', gap_nm));
        colorbar();
        
        % タイムスタンプ入りの画像ファイル名 (PNG 等)
        figNameBest = sprintf('result/best_structure_gap_%d_%s_gap_%d.png', gap_nm, timestamp, gap_nm);
        saveas(bestFig, figNameBest);  % 画像保存
    end
end

% Plot G_best for each gap size and save the figure
figure;
plot(gap_nm_values, G_best_values, '-o');
xlabel('Gap size (nm)');
ylabel('G\_best');
title('G\_best for each gap size');
grid on;
saveas(gcf, sprintf('result/G_best_vs_gap_size_%s.png', timestamp));

% Plot G_best * gap for each gap size and save the figure
figure;
plot(gap_nm_values, G_best_times_gap_values, '-o');
xlabel('Gap size (nm)');
ylabel('Power = G\_best * Gap');
title('Power = G\_best * Gap for each gap size');
grid on;
saveas(gcf, sprintf('result/G_best_times_gap_vs_gap_size_%s.png', timestamp));

