addpath(genpath('./'));                     % add the whole directory to path, if not already done

%% SET PARAMETERS
c0 = 1;                                     % speed of light m/s (normalized to 1)
lambda0 = 2;                                % central wavelength (um)

skip = 4;                                   % number of iteration frames between plots (higher->faster, lower->more plots)
display_plots = false;                      % plotting during the run?

alpha = 5e2;                                % step size in permittivity (~1e2-1e4 works well)
a = 3;                                      % smooth-max weight factor (see paper)
beta = 0.5;                                 % ratio of electron speed to speed of light

in_material = false;                        % evaluate E_max in material? or in surrounding regions.
starting = 0;                               % 0 -> vacuum, 1 -> random, 2 -> midway epsilon

grids_in_lam = 100;                         % number of grid points in a free space wavelength
gap_nm_values = 40:40:1000;                 % gap size in nm variations
% gap_nm_values = 400;                 % gap size in nm variations

N = 4000;                                   % number of iterations
% N = 100;                                   % number of iterations
parpool('local', 10);
timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
L = 0.4;                                    % size of optimization region (um)
npml = 10;                                  % number of PML points (need > 10)

% relative permittivity of material region. Uncomment to select
eps = 3.4363^2;     % Si 2um
%nmax = sqrt(eps);   % refractive index (not used explicitly here)
gamma = 0.9;        % momentum term

output_folder_name = 'result/single_channel_feb16';

%% SET OTHER CONSTANTS (DON'T CHANGE)
dlx = lambda0/grids_in_lam;    % grid size along electron trajectory axis
dly = dlx;                     % spacing in the perpendicular direction

% 結果保存用配列
num_gap = length(gap_nm_values);
G_best_values = zeros(num_gap, 1);
G_best_times_gap_values = zeros(num_gap, 1);
E_max_values = zeros(num_gap, 1);   % 各gapにおける最適化後のE_max

% parforループで各gapについて最適化を実施
parfor idx = 1:num_gap
    gap_nm = gap_nm_values(idx);
    
    pos_src = floor(npml+grids_in_lam/4);       % grid points from left edge to source
    spc_pts = floor(grids_in_lam/4);            % grid points between source and structure
    gap_pts = floor(gap_nm/1000/dlx);           % grid points corresponding to gap
    Lpts = round(L/dlx);                        % grid points in optimization region
    
    Nx = ceil(lambda0*beta/dlx);                % grid points in x-direction
    Ny = gap_pts + 2*(pos_src + Lpts + spc_pts);  % grid points in y-direction
    nx = floor(Nx/2);
    ny = floor(Ny/2);
    
    % 1チャネルの場合、min_G_Emaxのループは1回のみ（[0]として扱う）
    for min_G_Emax = 0  %%% ループ本体（1回のみ）
        
        %% 初期設定：FDFD用パラメータ
        ER = ones(Nx,Ny);            % relative permittivity grid
        MuR = ones(Nx,Ny);           % relative permeability grid
        ER_best = ones(Nx,Ny);       % best found structure
        A_best = 0;
        
        % TFSF領域の設定
        b = zeros(Nx,Ny);
        b(:, pos_src:pos_src + spc_pts + Lpts + gap_pts + Lpts + spc_pts) = 1;
        kinc = [0,1];              % 入射平面波方向
        
        RES = [dlx,dly];           % grid resolution
        BC = [-1,-1];              % boundary conditions
        NPML = [0,0,npml,npml];     % PML設定
        Pol = 'Hz';                % polarization
        spc = spc_pts*dly;         % sourceと構造間隔（um）
        gap = gap_pts*dly;         % gapサイズ（um）
        
        xs = dlx*(1:Nx);
        
        % 最適化領域のマスク
        delta_device = zeros(Nx,Ny);
        delta_device(1:Nx, pos_src + spc_pts : pos_src + spc_pts + Lpts) = 1;
        delta_device(1:Nx, pos_src + spc_pts + Lpts + gap_pts : pos_src + spc_pts + Lpts + gap_pts + Lpts) = 1;
        delta_device_vec = delta_device(:);
        
        % etaベクトル場の定義（加速子の入力としての重み付け）
        eta = zeros(Nx,Ny);
        eta(:,ny) = 1/Nx * exp(2*pi*1i*dlx*(0:Nx-1)/lambda0/beta);
        eta_vec = eta(:);
        
        % 初期ERの設定
        for i = 1:Nx
            for j = 1:Ny
                if delta_device(i,j) == 1
                    if starting == 1
                        ER(i,j) = rand*(eps-1)+1;
                    elseif starting == 2
                        ER(i,j) = eps/2+0.5;
                    end
                end
            end
        end
        
        % 基準電場 E0 の計算（全空間でのシミュレーション）
        [fields, ~] = FDFD_TFSF(ones(Nx,Ny), MuR, RES, NPML, BC, lambda0, Pol, b, kinc);
        Ex = fields.Ex;
        Ey = fields.Ey;
        E0 = sqrt(abs(Ex(nx, ny))^2 + abs(Ey(nx, ny))^2);
        
        % 反復最適化の初期化
        G_best = 0;    % 最良勾配
        Gs = zeros(N,1);
        E_maxs = zeros(N,1);
        G_by_Es = zeros(N,1);
        G_by_Sa = zeros(N,1);
        phis = zeros(N,1);
        AVM_prev = zeros(Nx,Ny);
        
        if display_plots
            figure(1);
        end
        
        if ~min_G_Emax
            display('working on gradient maximized structure');
        else
            display('working on acceleration factor maximized structure');
        end
        upd = textprogressbar(N);
        
        for j = 1:N
            upd(j);
            [fields, extra] = FDFD_TFSF(ER, MuR, RES, NPML, BC, lambda0, Pol, b, kinc);
            Ex = fields.Ex / E0;
            Ey = fields.Ey / E0;
            
            % 勾配計算
            g = sum(sum(eta .* Ex));
            G = real(g);
            phis(j) = angle(g);
            
            DEY = extra.derivatives.DEY;
            DEX = extra.derivatives.DEX;
            ER_vec = ER(:);
            chi = delta_device .* (ER - ones(Nx,Ny));
            if in_material
                E_abs = (chi/(eps-1)) .* sqrt(abs(Ex).^2 + abs(Ey).^2);
            else
                E_abs = delta_device .* sqrt(abs(Ex).^2 + abs(Ey).^2);
            end
            E_abs_vec = E_abs(:);
            E_maxs(j) = max(E_abs_vec);
            
            Ox = -1i*lambda0/(2*pi*c0) * spdiags(1./ER_vec,0,Nx*Ny,Nx*Ny) * DEY;
            Oy =  1i*lambda0/(2*pi*c0) * spdiags(1./ER_vec,0,Nx*Ny,Nx*Ny) * DEX;
            eta_aj = [eta_vec; zeros(Nx*Ny,1)];
            % ここでは1チャネルなので adjoint source は -eta_aj
            b_aj = -eta_aj;
            b_aj = reshape(Ox * b_aj(1:Nx*Ny) + Oy * b_aj(Nx*Ny+1:end), [Nx,Ny]);
            b_aj(isnan(b_aj)) = 0;
            
            AF = extra.AF;
            [fields_aj, ~] = FDFD_fast(ER, MuR, RES, NPML, BC, lambda0, Pol, b_aj, AF);
            x_aj = fields_aj.x / E0;
            Ex_aj = reshape(x_aj(1:Nx*Ny), [Nx,Ny]);
            Ey_aj = reshape(x_aj(Nx*Ny+1:end), [Nx,Ny]);
            
            % 感度情報の計算
            AVM = -real((Ex .* Ex_aj .* delta_device + Ey .* Ey_aj .* delta_device));
            ER = ER + alpha*AVM + alpha*gamma*AVM_prev;
            AVM_prev = AVM;
            
            ER(ER < 1) = 1;
            ER(ER > eps) = eps;
            
            if G > G_best
                G_best = G;
                ER_best = ER;
            end
            Gs(j) = G;
            
            if display_plots && mod(j, skip) == 0
                clf;
                subplot(2,2,1);
                disp = [];
                for kk = 1:5
                    disp = [disp; real(ER)];
                end
                imagesc(disp, [1, eps]);
                colormap(flipud(gray));
                title('relative permittivity');
                xlabel('pixel'); ylabel('pixel');  %%% NEW: 軸ラベルを追加
                set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal');
                set(gca,'FontSize',22,'fontWeight','normal');
                colorbar();
                
                subplot(2,2,2);
                plot(Gs(1:j),'k');
                xlabel('iteration number');
                ylabel('gradient (E_0)');
                title('acceleration gradient at \phi = 0');
                set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal');
                set(gca,'FontSize',22,'fontWeight','normal');
                colorbar();
                
                subplot(2,2,3);
                plot(1:j, G_by_Es(1:j));
                hold on;
                plot(1:j, G_by_Sa(1:j));
                xlabel('iteration number');
                ylabel('G/|E|max');
                title('acceleration factor');
                legend({'actual','smooth-max'});
                set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal');
                set(gca,'FontSize',22,'fontWeight','normal');
                
                subplot(2,2,4); hold on;
                plot(1:j, phis(1:j));
                plot(1:j, zeros(j,1));
                xlabel('iteration number');
                ylabel('\phi');
                legend({'computed','\phi=0 (target)'});
                title('acceleration phase (\phi)');
                set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal');
                set(gca,'FontSize',22,'fontWeight','normal');
                
                pause(0.001);
            end
        end
        
        %% POST PROCESSING
        % 二値化処理
        eps_avg = (eps+1)/2;
        ER(ER < eps_avg) = 1;
        ER(ER >= eps_avg) = eps;
        ER_best(ER_best < eps_avg) = 1;
        ER_best(ER_best >= eps_avg) = eps;
        
        % 最適化後構造による再シミュレーション
        [fields_best, extra_best] = FDFD_TFSF(ER_best, MuR, RES, NPML, BC, lambda0, Pol, b, kinc);
        Ex_best = fields_best.Ex / E0;
        Ey_best = fields_best.Ey / E0;
        g_best = sum(sum(eta .* Ex_best));
        G_best = abs(g_best);
        E_abs_best = delta_device .* sqrt(abs(Ex_best).^2 + abs(Ey_best).^2);
        E_max_best = max(E_abs_best(:));
        
        % 個別ファイル出力（txt形式での個別出力は残す）
        fname = sprintf('%s/final_acceleration_gradients_gap_%d_%s_gap_%d.txt', output_folder_name, gap_nm, timestamp, gap_nm);
        fileID = fopen(fname, 'w');
        fprintf(fileID, 'G_best: %f\n', G_best);
        fprintf(fileID, 'g_best: %f + %fi\n', real(g_best), imag(g_best));
        fprintf(fileID, 'E_max: %f\n', E_max_best);
        fprintf(fileID, 'L: %f\n', L);
        fclose(fileID);
        fprintf('File saved as: %s\n', fname);
        
        % 結果保存用配列に格納
        G_best_values(idx) = G_best;
        G_best_times_gap_values(idx) = G_best * gap_nm;
        E_max_values(idx) = E_max_best;
        
        % Best Structure の画像出力（5回縦連結して表示）
        if display_plots
            bestFig = figure('Name','Best Structure','Visible','on');
        else
            bestFig = figure('Name','Best Structure','Visible','off');
        end
        disp_best = [];
        for kk = 1:5
            disp_best = [disp_best; real(ER_best)];
        end
        imagesc(disp_best, [1, eps]);
        colormap(flipud(gray));
        axis equal tight;
        colorbar();
        xlabel('pixel'); ylabel('pixel');   %%% NEW: 軸ラベル追加
        title(sprintf('Best Structure (gap = %d nm)', gap_nm));
        figNameBest = sprintf('%s/best_structure_gap_%d_%s_gap_%d.png', output_folder_name, gap_nm, timestamp, gap_nm);
        saveas(bestFig, figNameBest);
        close(bestFig);
        
        %%% NEW: 電場分布の画像出力（最適化後構造での電場分布を5回連結して出力）
        [fields_best_for_plot, ~] = FDFD_TFSF(ER_best, MuR, RES, NPML, BC, lambda0, Pol, b, kinc);
        Ex_best_plot = fields_best_for_plot.Ex;
        Ey_best_plot = fields_best_for_plot.Ey;
        E_best_magnitude = sqrt(abs(Ex_best_plot).^2 + abs(Ey_best_plot).^2);
        efieldFig = figure('Name','Electric Field Magnitude','Visible','off');
        disp_efield = [];
        for kk = 1:5
            disp_efield = [disp_efield; E_best_magnitude];
        end
        imagesc(disp_efield);
        axis equal tight;
        colormap jet;
        colorbar();
        xlabel('pixel'); ylabel('pixel');   %%% NEW: 軸ラベル追加
        title(sprintf('Electric Field Magnitude (gap = %d nm)', gap_nm));
        efieldFigName = sprintf('%s/electric_field_gap_%d_%s_gap_%d.png', output_folder_name, gap_nm, timestamp, gap_nm);
        saveas(efieldFig, efieldFigName);
        close(efieldFig);
    end
end

% 各gapにおけるG_best, G_best*gap, E_maxをプロットして保存
figure;
plot(gap_nm_values, G_best_values, '-o');
xlabel('gap (nm)');
ylabel('abs(g)');
title('abs(g) vs gap');
grid on;
saveas(gcf, sprintf('%s/G_best_vs_gap_size_%s.png', output_folder_name, timestamp));

figure;
plot(gap_nm_values, G_best_times_gap_values, '-o');
xlabel('gap (nm)');
ylabel('abs(g) * gap');
title('abs(g) * gap vs gap');
grid on;
saveas(gcf, sprintf('%s/G_best_times_gap_vs_gap_size_%s.png', output_folder_name, timestamp));

%%% NEW: 全gapのサマリーをCSV形式で出力
summary_filename = sprintf('%s/summary_all_%s.csv', output_folder_name, timestamp);
fid_summary = fopen(summary_filename, 'w');
fprintf(fid_summary, 'gap_nm,G_best,abs(g)*gap,E_max\n');
for idx = 1:num_gap
    fprintf(fid_summary, '%d,%f,%f,%f\n', gap_nm_values(idx), G_best_values(idx), G_best_times_gap_values(idx), E_max_values(idx));
end
fclose(fid_summary);
fprintf('Summary CSV saved as: %s\n', summary_filename);

delete(gcp('nocreate'));
