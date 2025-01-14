addpath(genpath('./'));                     % add the whole directory to path, if not already done

%% SET PARAMETERS
c0 = 1;                                     % speed of light m/s (normalized to 1)
lambda0 = 2;                                % central wavelength (um)

skip = 4;                                   % number of iteration frames between plots (higher->faster, lower->more plots)
display_plots = false;                      % plotting during the run? (false にするとiteration中の表示を行わない)

alpha = 5e2;                                % step size in permittivity (~1e2-1e4 works well)
a = 3;                                      % smooth-max weight factor (see paper)
beta = 0.5;                                 % ratio of electron speed to speed of light
N = 800;                                    % number of iterations

in_material = false;                        % evaluate E_max in material? or in surrounding regions.
starting = 0;                               % 0 -> vacuum, 1 -> random, 2 -> midway epsilon

grids_in_lam = 50;                          % number of grid points in a free space wavelength
npml = 10;                                  % number of PML (absorbing region) points (need > 10 at least)

% relative permittivity of material region.  uncomment to select
eps = 3.4363^2;     % Si 2um
%eps = 1.4381^2;    % fused silica 2um
%eps = 1.9834^2;    % Si3N4
%eps = 1.9^2;       % GaOx

nmax = sqrt(eps);                           % refractive index of material region

gamma = 0.9;                                % 'momentum term', see paper. 0-1

%% 新たに追加: gap を変化させるための配列
gap_nm_values = 100:10:1300;

%% 各 gap に対する最終的な G_best を格納する配列
G_best_values            = [];
G_best_times_gap_times2  = [];  % G_best * gap * 2
Gsum_times_gap_values    = [];  % (G1_best + G2_best)*gap

%% 出力フォルダ名を設定
output_folder_name = 'result/2_channel_step_10_jan13';

%% ループ開始
for gap_nm = gap_nm_values
    %% SET OTHER CONSTANTS (DON'T CHANGE)
    dlx = lambda0/grids_in_lam;                 % grid size along electron trajectory axis
    dly = dlx;                                  % spacing in the perpendicular direction
    
    % gap_nm から grid point に換算
    gap_pts = floor(gap_nm/1000/dlx);           % number of grid points in the gap
    
    % ここでは「2つのギャップ + 中央ギャップ (gap_gap_nm)」のようにしていたコードを
    % そのまま残していますが，適宜変更してください．
    gap_gap_nm = 200;                           % 例として固定 (2つのギャップの間のギャップ)
    gap_gap_pts = floor(gap_gap_nm/1000/dlx);   % number of grid points in the gap between the two gaps
    
    L = 1.0;                                    % size of optimization region (um)
    Lpts = round(L/dlx);                        % number of points in the optimization region
    
    pos_src = floor(npml+grids_in_lam/4);       % number of grid points between left edge and source
    spc_pts = floor(grids_in_lam/4);            % number of grid points between source and structure
    
    Nx = ceil(lambda0*beta/dlx);
    % 2つのギャップ + 中央 gap_gap_pts + 上下2つの最適化領域 + PML の外の領域 など
    Ny = 2*gap_pts + 2*(pos_src + Lpts + spc_pts) + gap_gap_pts;
    
    nx = floor(Nx/2);
    ny1 = floor(gap_pts/2 + pos_src + Lpts + spc_pts);
    ny2 = floor(gap_pts + gap_pts/2 + gap_gap_pts + pos_src + Lpts + spc_pts);
    
    % First compute G maximization, then do G/E_max maximization (for comparison)
    % for min_G_Emax = (0:1)
    for min_G_Emax = 0
        
        %% This section defines the input parameters that my FDFD code needs to run.
        ER  = ones(Nx,Ny);
        MuR = ones(Nx,Ny);
        ER_best = ones(Nx,Ny);
        A_best = 0;
        
        b = zeros(Nx,Ny);
        % define the total field region on the grid
        b(:, pos_src:pos_src + spc_pts + Lpts + gap_pts + gap_gap_pts + gap_pts + Lpts + spc_pts) = 1;
        kinc = [0,1];
        RES = [dlx,dly];
        BC = [-1,-1];
        NPML = [0,0,npml,npml];
        Pol= 'Hz';
        spc = spc_pts*dly;
        gap = gap_pts*dly;
        
        xs = dlx*(1:Nx);
        
        delta_device = zeros(Nx,Ny);
        delta_device(1:Nx, pos_src + spc_pts : pos_src + spc_pts + Lpts) = 1;
        delta_device(1:Nx, pos_src + spc_pts + Lpts + gap_pts : pos_src + spc_pts + Lpts + gap_pts + gap_gap_pts) = 1;
        delta_device(1:Nx, pos_src + spc_pts + Lpts + gap_pts + gap_gap_pts + gap_pts : pos_src + spc_pts + Lpts + gap_pts + gap_gap_pts + gap_pts + Lpts) = 1;
        delta_device_vec = delta_device(:);
        
        % define the eta vector fields for the two channels
        eta1 = zeros(Nx,Ny);
        eta1(:,ny1) = 1/Nx*exp(2*pi*1i*dlx*(0:Nx-1)/lambda0/beta);
        eta1_vec = eta1(:);
        
        eta2 = zeros(Nx,Ny);
        eta2(:,ny2) = 1/Nx*exp(2*pi*1i*dlx*(0:Nx-1)/lambda0/beta);
        eta2_vec = eta2(:);
        
        % define stating permittivity
        for i = (1:Nx)
            for j = (1:Ny)
                if (delta_device(i,j) == 1)
                    if (starting == 1)
                        ER(i,j) = rand*(eps-1)+1;
                    elseif (starting == 2)
                        ER(i,j) = eps/2+0.5;
                    else
                        % starting=0 -> vacuum
                    end
                end
            end
        end
        
        % run the simulation with accelerator input (plane wave) but all empty space
        [fields, ~] = FDFD_TFSF(ones(Nx,Ny),MuR,RES,NPML,BC,lambda0,Pol,b,kinc);
        
        % get the fields and the E0 (normalization)
        Ex = fields.Ex;
        Ey = fields.Ey;
        E0 = sqrt(abs(Ex(nx, ny1))^2 + abs(Ey(nx, ny1))^2);
        
        % define variables to store the iteration progress
        G_best_local = 0;          % best gradient in this run
        Gs = zeros(N,1);
        E_maxs = zeros(N,1);
        G_by_Es = zeros(N,1);
        G_by_Sa = zeros(N,1);
        
        phis = zeros(N,1);
        phi = 0;
        AVM_prev = zeros(Nx,Ny);
        
        %---- 変更点3: display_plots が false なら iteration中のウィンドウ表示は行わない
        if display_plots
            figure(1);  % open a figure to plot
        end
        
        if ~min_G_Emax
            display('working on gradient maximized structure');
        else
            display('working on acceleration factor maximized structure');
        end
        upd = textprogressbar(N);
        
        for j = (1:N)
            
            upd(j);
            % original simulation
            [fields, extra] = FDFD_TFSF(ER,MuR,RES,NPML,BC,lambda0,Pol,b,kinc);
            Ex = fields.Ex/E0;
            Ey = fields.Ey/E0;
            
            % compute gradients
            g1 = sum(sum(eta1.*Ex));
            G1 = real(g1);
            g2 = sum(sum(eta2.*Ex));
            G2 = real(g2);
            g = g1 + g2;
            G = real(g);
            
            % get phase
            phis(j) = angle(g);
            
            % get numerical spatial derivative operators
            DEY = extra.derivatives.DEY;
            DEX = extra.derivatives.DEX;
            
            ER_vec = ER(:);
            chi = delta_device.*(ER - ones(Nx,Ny));
            
            Ox = -1i*lambda0/2/pi/c0*spdiags(1./ER_vec,0,Nx*Ny,Nx*Ny)*DEY;
            Oy =  1i*lambda0/2/pi/c0*spdiags(1./ER_vec,0,Nx*Ny,Nx*Ny)*DEX;
            
            eta1_aj = [eta1_vec; zeros(Nx*Ny,1)];
            eta2_aj = [eta2_vec; zeros(Nx*Ny,1)];
            
            % E_max は in_material を考慮
            if (in_material)
                E_abs = (chi/(eps-1)).*sqrt(abs(Ex).^2 + abs(Ey).^2);
            else
                E_abs = delta_device.*sqrt(abs(Ex).^2 + abs(Ey).^2);
            end
            
            x_abs = E_abs(:);
            alpha_vec = exp(x_abs*a);
            alpha_T_1 = sum(alpha_vec);
            Sa = sum(alpha_vec.*x_abs)/alpha_T_1;
            
            x = [Ex(:); Ey(:)];
            z = conj(x./[x_abs;x_abs]);
            z(isnan(z)) = 0;
            z(isinf(z)) = 0;
            spdiagz = spdiags(z,0,Nx*Ny*2,Nx*Ny*2);
            P = [speye(Nx*Ny) speye(Nx*Ny)];
            
            S = real(1/alpha_T_1*(speye(Nx*Ny) + a*spdiags(x_abs,0,Nx*Ny,Nx*Ny) ...
                - a*sum(alpha_vec.*x_abs)/alpha_T_1*speye(Nx*Ny)));
            sigma = transpose(alpha_vec)*S*(P*spdiagz);
            sigma(isnan(sigma)) = 0;
            
            b_aj1 = transpose(G/Sa^2 * sigma);
            b_aj2 = -eta1_aj/Sa - eta2_aj/Sa;
            
            if (min_G_Emax)
                b_aj = b_aj1 + b_aj2;
            else
                b_aj = b_aj2;
            end
            b_aj = reshape(Ox*b_aj(1:Nx*Ny) + Oy*b_aj(Nx*Ny+1:end),[Nx,Ny]);
            b_aj(isnan(b_aj)) = 0 ;
            
            AF = extra.AF;
            [fields_aj, ~] = FDFD_fast(ER,MuR,RES,NPML,BC,lambda0,Pol,b_aj,AF);
            
            x_aj = fields_aj.x/E0;
            Ex_aj = reshape(x_aj(1:Nx*Ny),[Nx,Ny]);
            Ey_aj = reshape(x_aj(Nx*Ny+1:end),[Nx,Ny]);
            
            AVM = -real((Ex.*Ex_aj.*delta_device + Ey.*Ey_aj.*delta_device));
            
            % record relevant variables
            E_max = max(max(E_abs));
            E_maxs(j) = E_max;
            Gs(j) = G;
            G_by_Es(j) = G/E_max;
            G_by_Sa(j) = G/Sa;
            
            % update permittivity
            ER = ER + alpha*AVM + alpha*gamma*AVM_prev;
            AVM_prev = AVM;
            
            ER(ER < 1) = 1;
            ER(ER > eps) = eps;
            
            if (G > G_best_local)
                G_best_local = G;
                ER_best = ER;
            end
            
            %---- 変更点3: plotting during iteration は display_plots が true の時だけ
            if display_plots && mod(j,skip) == 0
                clf;
                subplot(2,2,1);
                disp_map = [];
                for k_ = (1:5)
                    disp_map = [disp_map; real(ER)];
                end
                imagesc(disp_map,[1,eps])
                colormap(flipud(gray))
                title('relative permittivity')
                set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal')
                set(gca,'FontSize',22,'fontWeight','normal')
                colorbar()
                
                subplot(2,2,2);
                plot(Gs(1:j),'k');
                xlabel('iteration number')
                ylabel('power (G)')
                title('acceleration gradient at \phi = 0')
                set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal')
                set(gca,'FontSize',22,'fontWeight','normal')
                colorbar()
                
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
        
        % force binary
        eps_avg = (eps+1)/2;
        ER(ER<eps_avg) = 1;
        ER(ER>=eps_avg) = eps;
        ER_best(ER_best<eps_avg) = 1;
        ER_best(ER_best>=eps_avg) = eps;
        
        % do another simulation of the binary distribution for ER_best
        [fields_best, extra_best] = FDFD_TFSF(ER_best,MuR,RES,NPML,BC,lambda0,Pol,b,kinc);
        Ex_best = fields_best.Ex/E0;
        Ey_best = fields_best.Ey/E0;
        
        % calculate g1_best and g2_best after the loop
        g1_best = sum(sum(eta1.*Ex_best));
        g2_best = sum(sum(eta2.*Ex_best));
        
        G1_best = real(g1_best);
        G2_best = real(g2_best);
        
        g_best = g1_best + g2_best;
        % 最終的に abs で取るかは元のコードの通り
        G_best_local = abs(g_best);
        
        % 結果を保存
        timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
        fname = sprintf('%s/final_acceleration_gradients_gap_%d_%s.txt', output_folder_name, gap_nm, timestamp);
        fileID = fopen(fname, 'w');
        fprintf(fileID, 'G_best (abs): %f\n', G_best_local);
        fprintf(fileID, 'G1_best: %f\n', G1_best);
        fprintf(fileID, 'G2_best: %f\n', G2_best);
        fprintf(fileID, 'g_best: %f + %fi\n', real(g_best), imag(g_best));
        fprintf(fileID, 'g1_best: %f + %fi\n', real(g1_best), imag(g1_best));
        fprintf(fileID, 'g2_best: %f + %fi\n', real(g2_best), imag(g2_best));
        fclose(fileID);
        fprintf('File saved as: %s\n', fname);
        
        % best structure の可視化・保存
        if display_plots
            bestFig = figure('Name','Best Structure','Visible','on');
        else
            bestFig = figure('Name','Best Structure','Visible','off');
        end
        
        disp_best = [];
        for k_ = 1:5
            disp_best = [disp_best; real(ER_best)];
        end
        imagesc(disp_best, [1, eps]);
        colormap(flipud(gray));
        axis equal tight;
        title(sprintf('Best Structure (gap = %d nm)', gap_nm));
        colorbar();
        
        figNameBest = sprintf('%s/best_structure_gap_%d_%s.png', output_folder_name, gap_nm, timestamp);
        saveas(bestFig, figNameBest);
        
        %---- ここで今回の gap に対する G_best_local, G1_best, G2_best を記録
        G_best_values           = [G_best_values; G_best_local];
        % 2チャネル分 => G_best * gap_nm * 2
        G_best_times_gap_times2 = [G_best_times_gap_times2; G_best_local * gap_nm * 2];
        % (G1_best + G2_best)*gap_nm
        Gsum_times_gap_values   = [Gsum_times_gap_values; (G1_best + G2_best)*gap_nm];
        
    end % end of min_G_Emax loop
    
end % end of gap_nm loop


%% gap を x軸として，以下の3種をプロット
% (a) G_best vs gap
figure;
plot(gap_nm_values, G_best_values, '-o');
xlabel('gap (nm)');
ylabel('G\_best');
title('G\_best vs. gap (2-channel)');
grid on;
saveas(gcf, sprintf('%s/G_best_vs_gap_multi_channel_%s.png', output_folder_name, datestr(now,'yyyy-mm-dd_HHMMSS')));

% (b) G_best * gap * 2 vs gap
figure;
plot(gap_nm_values, G_best_times_gap_times2, '-o');
xlabel('gap (nm)');
ylabel('G\_best * gap * 2');
title('G\_best * gap * 2 vs. gap (2-channel)');
grid on;
saveas(gcf, sprintf('%s/G_best_times_gap_times2_vs_gap_multi_channel_%s.png', output_folder_name, datestr(now,'yyyy-mm-dd_HHMMSS')));

% (c) (G1_best + G2_best) * gap vs gap
figure;
plot(gap_nm_values, Gsum_times_gap_values, '-o');
xlabel('gap (nm)');
ylabel('(G1\_best + G2\_best) * gap');
title('(G1\_best + G2\_best) * gap vs. gap (2-channel)');
grid on;
saveas(gcf, sprintf('%s/Gsum_times_gap_vs_gap_multi_channel_%s.png', output_folder_name, datestr(now,'yyyy-mm-dd_HHMMSS')));

