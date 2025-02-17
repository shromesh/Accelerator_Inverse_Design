addpath(genpath('./'));  % add the whole directory to path, if not already done

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
% gap_nm_values = 40:40:1000;                 % gap size in nm variations
% gap_nm_values = 360;                 % gap size in nm variations
gap_nm_values = 520;                 % gap size in nm variations

N = 4000;                                   % number of iterations
parpool('local', 10);
timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
L = 0.4;                                    % size of optimization region (um)
npml = 10;                                  % number of PML points (need > 10)

% relative permittivity of material region. Uncomment one.
eps = 3.4363^2;     % Si 2um
% nmax = sqrt(eps);
gamma = 0.9;        % momentum term

output_folder_name = 'result/single_channel_feb17_x_field';

%% SET OTHER CONSTANTS (DON'T CHANGE)
dlx = lambda0/grids_in_lam;    % grid size along electron trajectory axis
dly = dlx;

num_gap = length(gap_nm_values);

%%% NEW %%% 事前計算：設計領域（src, spc 領域を除く）の左右端インデックスを各 gap 用に計算
% pos_src, spc_pts, Lpts は gap_nm に依存せず共通
pos_src = floor(npml+grids_in_lam/4);
spc_pts = floor(grids_in_lam/4);
Lpts = round(L/dlx);
% 設計領域は、ソース領域(spc_pts)を除いた部分なので
y_design_start_all = pos_src + spc_pts;  % これは全 gap で同じ
y_design_end_all = zeros(num_gap,1); % いらんかったかも．GPTの出力参照
for i = 1:num_gap
    gap_pts = floor(gap_nm_values(i)/1000/dlx);
    % 設計領域の右端は、src+spcから Lpts, gap_pts, Lpts 分だけ延長した位置
    y_design_end_all(i) = pos_src + spc_pts + Lpts + gap_pts + Lpts;
end
%%% NEW %%%

%% SETUP FOR OPTIMIZATION (保存用配列)
G_best_values = zeros(num_gap, 1);
G_best_times_gap_values = zeros(num_gap, 1);
E_max_values = zeros(num_gap, 1);

% parfor ループで各 gap について最適化
parfor idx = 1:num_gap
    gap_nm = gap_nm_values(idx);
    
    pos_src = floor(npml+grids_in_lam/4);
    spc_pts = floor(grids_in_lam/4);
    gap_pts = floor(gap_nm/1000/dlx);
    Lpts = round(L/dlx);
    
    Nx = ceil(lambda0*beta/dlx);
    Ny = gap_pts + 2*(pos_src + Lpts + spc_pts);
    nx = floor(Nx/2);
    ny = floor(Ny/2);
    
    % 初期設定
    ER = ones(Nx,Ny);
    MuR = ones(Nx,Ny);
    ER_best = ones(Nx,Ny);
    b = zeros(Nx,Ny);
    b(:, pos_src:pos_src+spc_pts+Lpts+gap_pts+Lpts+spc_pts) = 1;
    kinc = [0,1];
    
    RES = [dlx,dly];
    BC = [-1,-1];
    NPML = [0,0,npml,npml];
    Pol = 'Hz';
    spc = spc_pts*dly;
    gap = gap_pts*dly;
    
    xs = dlx*(1:Nx);
    
    % 最適化領域のマスク
    delta_device = zeros(Nx,Ny);
    delta_device(1:Nx, pos_src+spc_pts : pos_src+spc_pts+Lpts) = 1;
    delta_device(1:Nx, pos_src+spc_pts+Lpts+gap_pts : pos_src+spc_pts+Lpts+gap_pts+Lpts) = 1;
    delta_device_vec = delta_device(:);
    
    % eta ベクトル場の定義（加速子入力としての重み付け）
    eta = zeros(Nx,Ny);
    eta(:,ny) = 1/Nx * exp(2*pi*1i*dlx*(0:Nx-1)/lambda0/beta);
    eta_vec = eta(:);
    
    % 初期ERの設定
    for i = 1:Nx
        for j = 1:Ny
            if delta_device(i,j)==1
                if starting==1
                    ER(i,j) = rand*(eps-1)+1;
                elseif starting==2
                    ER(i,j) = eps/2+0.5;
                end
            end
        end
    end
    
    % 基準電場 E0 の計算（全空間シミュレーション）
    [fields, ~] = FDFD_TFSF(ones(Nx,Ny), MuR, RES, NPML, BC, lambda0, Pol, b, kinc);
    Ex = fields.Ex;
    Ey = fields.Ey;
    E0 = sqrt(abs(Ex(nx,ny))^2 + abs(Ey(nx,ny))^2);
    
    % 反復最適化初期化
    G_best = 0;
    Gs = zeros(N,1);
    E_maxs = zeros(N,1);
    phis = zeros(N,1);
    AVM_prev = zeros(Nx,Ny);
    
    for j = 1:N
        [fields, extra] = FDFD_TFSF(ER, MuR, RES, NPML, BC, lambda0, Pol, b, kinc);
        Ex = fields.Ex/E0;
        Ey = fields.Ey/E0;
        
        % 勾配計算（x方向成分 Ex を重視）
        g = sum(sum(eta.*Ex));
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
        
        Ox = -1i*lambda0/(2*pi*c0)*spdiags(1./ER_vec,0,Nx*Ny,Nx*Ny)*DEY;
        Oy =  1i*lambda0/(2*pi*c0)*spdiags(1./ER_vec,0,Nx*Ny,Nx*Ny)*DEX;
        eta_aj = [eta_vec; zeros(Nx*Ny,1)];
        % 1チャネルなので、adjoint source は -eta_aj
        b_aj = -eta_aj;
        b_aj = reshape(Ox*b_aj(1:Nx*Ny)+Oy*b_aj(Nx*Ny+1:end), [Nx,Ny]);
        b_aj(isnan(b_aj)) = 0;
        
        AF = extra.AF;
        [fields_aj, ~] = FDFD_fast(ER, MuR, RES, NPML, BC, lambda0, Pol, b_aj, AF);
        x_aj = fields_aj.x/E0;
        Ex_aj = reshape(x_aj(1:Nx*Ny), [Nx,Ny]);
        Ey_aj = reshape(x_aj(Nx*Ny+1:end), [Nx,Ny]);
        
        AVM = -real((Ex.*Ex_aj.*delta_device + Ey.*Ey_aj.*delta_device));
        ER = ER + alpha*AVM + alpha*gamma*AVM_prev;
        AVM_prev = AVM;
        
        ER(ER<1) = 1;
        ER(ER>eps) = eps;
        
        if G > G_best
            G_best = G;
            ER_best = ER;
        end
        Gs(j) = G;
    end
    
    %% POST PROCESSING
    eps_avg = (eps+1)/2;
    ER_best(ER_best<eps_avg)=1;
    ER_best(ER_best>=eps_avg)=eps;
    % 保存用構造
    ER_best_1D{idx} = ER_best;
    
    [fields_best, extra_best] = FDFD_TFSF(ER_best, MuR, RES, NPML, BC, lambda0, Pol, b, kinc);
    Ex_best = fields_best.Ex/E0;
    Ey_best = fields_best.Ey/E0;
    g1_best = sum(sum(eta.*Ex_best));  % ※ single channel では g1 としてまとめている
    g_best  = g1_best;
    G_best_final = abs(g_best);
    
    E_abs = delta_device .* sqrt(abs(Ex_best).^2 + abs(Ey_best).^2);
    E_max = max(E_abs(:));
    
    G_best_values(idx) = G_best_final;
    G_best_times_gap_values(idx) = G_best_final * gap_nm;
    E_max_values(idx) = E_max;
    
    %%% NEW: 電場分布再計算（最適化後構造）
    [fields_best_for_plot, ~] = FDFD_TFSF(ER_best, MuR, RES, NPML, BC, lambda0, Pol, b, kinc);
    Ex_best_plot = fields_best_for_plot.Ex;
    Ey_best_plot = fields_best_for_plot.Ey;
    E_best_magnitude = sqrt(abs(Ex_best_plot).^2+abs(Ey_best_plot).^2);
    E_field_magnitude_1D{idx} = E_best_magnitude;
    %%% NEW: x方向電場 (Ex) の結果も保存
    Ex_field_1D{idx} = Ex_best_plot;
    
end  % end of parfor

%% ファイル出力＆画像保存部
% 再定義
dlx = lambda0/grids_in_lam;
for idx = 1:num_gap
    gap_nm = gap_nm_values(idx);
    
    %%% 個別テキスト出力
    fname = sprintf('%s/final_acceleration_gradients_gap_%d_%s.txt', output_folder_name, gap_nm, timestamp);
    fileID = fopen(fname, 'w');
    fprintf(fileID, 'G_best: %f\n', G_best_values(idx));
    fprintf(fileID, 'E_max: %f\n', E_max_values(idx));
    fclose(fileID);
    fprintf('File saved as: %s\n', fname);
    
    % 構造画像出力（5回連結）
    ER_best_k = ER_best_1D{idx};
    bestFig = figure('Name','Best Structure','Visible','off');
    disp_best = [];
    for kk = 1:5
        disp_best = [disp_best; real(ER_best_k)];
    end
    imagesc(disp_best, [1, eps]);
    colormap(flipud(gray));
    axis equal tight;
    colorbar();
    xlabel('pixel'); ylabel('pixel');
    title(sprintf('Best Structure (gap = %d nm)', gap_nm));
    figNameBest = sprintf('%s/best_structure_gap_%d_%s.png', output_folder_name, gap_nm, timestamp);
    saveas(bestFig, figNameBest);
    close(bestFig);
    
    %%% NEW: 電場大きさ画像出力（5回連結）【設計領域のみ表示】
    %%%% CROPPED: 設計領域は pos_src+spc_pts から (pos_src+spc_pts+Lpts+gap_pts+Lpts) まで
    y_design_start = pos_src + spc_pts;
    gap_pts = floor(gap_nm/1000/dlx);  % 再計算
    y_design_end = pos_src + spc_pts + Lpts + gap_pts + Lpts;
    
    E_best_magnitude = E_field_magnitude_1D{idx};
    E_best_magnitude_design = E_best_magnitude(:, y_design_start:y_design_end);
    efieldFig = figure('Name','Electric Field Magnitude','Visible','off');
    disp_efield = [];
    for kk = 1:5
        disp_efield = [disp_efield; E_best_magnitude_design];
    end
    imagesc(disp_efield);
    axis equal tight;
    colormap jet;
    colorbar();
    xlabel('pixel'); ylabel('pixel');
    title(sprintf('Electric Field Magnitude (gap = %d nm)', gap_nm));
    efieldFigName = sprintf('%s/electric_field_gap_%d_%s.png', output_folder_name, gap_nm, timestamp);
    saveas(efieldFig, efieldFigName);
    close(efieldFig);
    
    %%% NEW: Ex (x方向電場) の出力（5回連結）【設計領域のみ表示】
    Ex_best_plot = Ex_field_1D{idx};
    Ex_best_plot_design = Ex_best_plot(:, y_design_start:y_design_end);
    exFig = figure('Name','Electric Field Ex','Visible','off');
    disp_ex = [];
    for kk = 1:5
        disp_ex = [disp_ex; real(Ex_best_plot_design)];
    end
    imagesc(disp_ex);
    axis equal tight;
    colormap jet;
    colorbar();
    xlabel('pixel'); ylabel('pixel');
    title(sprintf('Electric Field Ex (gap = %d nm)', gap_nm));
    exFigName = sprintf('%s/electric_field_Ex_gap_%d_%s.png', output_folder_name, gap_nm, timestamp);
    saveas(exFig, exFigName);
    close(exFig);
end

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
