addpath(genpath('./'));                     % add the whole directory to path, if not already done

%% SET PARAMETERS
c0 = 1;                                     % speed of light m/s (normalized to 1)
lambda0 = 2;                                % central wavelength (um)

skip = 4;                                   % number of iteration frames between plots (higher->faster, lower->more plots)
display_plots = false;                      % plotting during the run? (false にするとiteration中の表示を行わない)

alpha = 5e2;                                % step size in permittivity (~1e2-e4 works well)
a = 3;                                      % smooth-max weight factor (see paper)
beta = 0.5;                                 % ratio of electron speed to speed of light

in_material = false;                        % evaluate E_max in material? or in surrounding regions.
starting = 0;                               % 0 -> vacuum, 1 -> random, 2 -> midway epsilon

grids_in_lam = 100;                         % number of grid points in a free space wavelength

%% 新たに追加: gap を変化させるための配列
gap_nm_values = 40:40:1000;
% gap_nm_values = [300, 400];

%% gap_gap を変化させるための配列
gap_gap_nm_values = 100:200:1000;
% gap_gap_nm_values = [300, 400];

N = 4000;                                   % number of iterations
% N = 100;                                   % number of iterations

parpool('local', 10);
timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
npml = 10;                                  % number of PML (absorbing region) points (need > 10 at least)

% relative permittivity of material region.  uncomment to select
eps = 3.4363^2;     % Si 2um
% eps = 1.4381^2;    % fused silica 2um
% eps = 1.9834^2;    % Si3N4
% eps = 1.9^2;       % GaOx

gamma = 0.9;                                % 'momentum term', see paper. 0-1

%% 出力フォルダ名を設定
output_folder_name = 'result/double_channel_feb16';

% -------------------------------------------------------------
% gap_nm_values, gap_gap_nm_values の長さ
ngap    = length(gap_nm_values);
ngapgap = length(gap_gap_nm_values);
nComb   = ngap * ngapgap; % 全組み合わせ数

% --- 事前に Nx, Ny を計算 (最初のパラメータセットで代表させる) ---
dlx_pre = lambda0/grids_in_lam;
gap_pts_pre = floor(gap_nm_values(1)/1000/dlx_pre);
gap_gap_pts_pre = floor(gap_gap_nm_values(1)/1000/dlx_pre);
Lpts_pre = round(0.4/dlx_pre);
pos_src_pre = floor(npml+grids_in_lam/4);
spc_pts_pre = floor(grids_in_lam/4);
Nx_pre = ceil(lambda0*beta/dlx_pre);
Ny_pre = 2*gap_pts_pre + 2*(pos_src_pre + Lpts_pre + spc_pts_pre) + gap_gap_pts_pre;
% -------------------------------------------------------------

% -------------------------------------------------------------
%  1D の配列として用意（結果をすべて格納）
G_best_values_1D             = zeros(nComb, 1);
G_best_abs_sums_1D           = zeros(nComb, 1);
G_best_values_times_gap_1D   = zeros(nComb, 1);
G_best_abs_sums_times_gap_1D = zeros(nComb, 1);
G_best_values_final_1D       = zeros(nComb, 1);
G1_best_1D                   = zeros(nComb, 1);
G2_best_1D                   = zeros(nComb, 1);
g_best_complex_1D            = complex(zeros(nComb, 1), zeros(nComb, 1));
g1_best_complex_1D           = complex(zeros(nComb, 1), zeros(nComb, 1));
g2_best_complex_1D           = complex(zeros(nComb, 1), zeros(nComb, 1));
E_max_1D                     = zeros(nComb, 1);
abs_g1_plus_g2_times_gap_1D  = zeros(nComb, 1);
abs_g_best_times_gap_1D      = zeros(nComb, 1);
ER_best_1D                   = cell(nComb, 1); % 構造格納用
abs_g1_best_1D               = zeros(nComb, 1);
abs_g2_best_1D               = zeros(nComb, 1);
abs_sum_g_times_gap_1D       = zeros(nComb, 1);

%%% NEW %%% % 電場分布（最適化後構造で計算）の結果を格納する cell 配列
E_field_magnitude_1D = cell(nComb, 1);
% -------------------------------------------------------------

% -------------------------------------------------------------
% ループ開始 (parfor)
parfor k = 1:nComb % 1次元の parfor ループに変更
    % --- 1次元インデックス k から (idx, jGapGap) を復元 ---
    idx     = floor((k-1)/ngapgap) + 1;
    jGapGap = mod(k-1, ngapgap) + 1;
    
    gap_nm     = gap_nm_values(idx);
    gap_gap_nm = gap_gap_nm_values(jGapGap);
    
    %% SET OTHER CONSTANTS
    dlx = lambda0/grids_in_lam;
    dly = dlx;
    gap_pts = floor(gap_nm/1000/dlx);
    gap_gap_pts = floor(gap_gap_nm/1000/dlx);
    L = 0.4;
    Lpts = round(L/dlx);
    pos_src = floor(npml+grids_in_lam/4);
    spc_pts = floor(grids_in_lam/4);
    Nx = ceil(lambda0*beta/dlx);
    Ny = 2*gap_pts + 2*(pos_src + Lpts + spc_pts) + gap_gap_pts;
    nx = floor(Nx/2);
    ny1 = floor(gap_pts/2 + pos_src + Lpts + spc_pts);
    ny2 = floor(gap_pts + gap_pts/2 + gap_gap_pts + pos_src + Lpts + spc_pts);
    
    ER  = ones(Nx,Ny);
    MuR = ones(Nx,Ny);
    ER_best = ones(Nx,Ny);
    b = zeros(Nx,Ny);
    b(:, pos_src:pos_src + spc_pts + Lpts + gap_pts + gap_gap_pts + gap_pts + Lpts + spc_pts) = 1;
    kinc = [0,1];
    RES = [dlx,dly];
    BC = [-1,-1];
    NPML = [0,0,npml,npml];
    Pol= 'Hz';
    spc = spc_pts*dly;
    gap = gap_pts*dly;
    xs = dlx*(1:Nx);
    
    % デバイス配置用のマスク
    delta_device = zeros(Nx,Ny);
    delta_device(1:Nx, pos_src + spc_pts : pos_src + spc_pts + Lpts) = 1;
    delta_device(1:Nx, pos_src + spc_pts + Lpts + gap_pts : pos_src + spc_pts + Lpts + gap_pts + gap_gap_pts) = 1;
    delta_device(1:Nx, pos_src + spc_pts + Lpts + gap_pts + gap_gap_pts + gap_pts : pos_src + spc_pts + Lpts + gap_pts + gap_gap_pts + gap_pts + Lpts) = 1;
    delta_device_vec = delta_device(:);
    
    % eta1, eta2 の設定
    eta1 = zeros(Nx,Ny);
    eta1(:,ny1) = 1/Nx*exp(2*pi*1i*dlx*(0:Nx-1)/lambda0/beta);
    eta1_vec = eta1(:);
    
    eta2 = zeros(Nx,Ny);
    eta2(:,ny2) = 1/Nx*exp(2*pi*1i*dlx*(0:Nx-1)/lambda0/beta);
    eta2_vec = eta2(:);
    
    % ER の初期化
    for i = 1:Nx
        for j = 1:Ny
            if delta_device(i,j) == 1
                if starting == 1
                    ER(i,j) = rand*(eps-1)+1;
                elseif starting == 2
                    ER(i,j) = eps/2+0.5;
                else
                    % starting = 0 -> vacuum
                end
            end
        end
    end
    
    % 基準電場 E0 の計算
    [fields, ~] = FDFD_TFSF(ones(Nx,Ny),MuR,RES,NPML,BC,lambda0,Pol,b,kinc);
    Ex = fields.Ex;
    Ey = fields.Ey;
    E0 = sqrt(abs(Ex(nx, ny1))^2 + abs(Ey(nx, ny1))^2);
    
    G_best_local = 0;
    AVM_prev = zeros(Nx,Ny);
    
    if display_plots
        figure(1);
    end
    
    % 最適化反復ループ
    display('working on gradient maximized structure');
    upd = textprogressbar(N);
    Gs     = zeros(N,1);
    E_maxs = zeros(N,1);
    phis   = zeros(N,1);
    
    for jj = 1:N
        upd(jj);
        [fields, extra] = FDFD_TFSF(ER,MuR,RES,NPML,BC,lambda0,Pol,b,kinc);
        Ex = fields.Ex/E0;
        Ey = fields.Ey/E0;
        g1 = sum(sum(eta1.*Ex));
        g2 = sum(sum(eta2.*Ex));
        g  = g1 + g2;
        phis(jj) = angle(g);
        
        DEY = extra.derivatives.DEY;
        DEX = extra.derivatives.DEX;
        ER_vec = ER(:);
        
        % 電界振幅 (E_max の算出用)
        chi = delta_device.*(ER - ones(Nx,Ny));
        if in_material
            E_abs = (chi/(eps-1)).*sqrt(abs(Ex).^2 + abs(Ey).^2);
        else
            E_abs = delta_device.*sqrt(abs(Ex).^2 + abs(Ey).^2);
        end
        E_abs_vec = E_abs(:);
        E_maxs(jj) = max(E_abs_vec);
        
        Ox = -1i*lambda0/(2*pi*c0)*spdiags(1./ER_vec,0,Nx*Ny,Nx*Ny)*DEY;
        Oy =  1i*lambda0/(2*pi*c0)*spdiags(1./ER_vec,0,Nx*Ny,Nx*Ny)*DEX;
        eta1_aj = [eta1_vec; zeros(Nx*Ny,1)];
        eta2_aj = [eta2_vec; zeros(Nx*Ny,1)];
        b_aj = - (eta1_aj + eta2_aj);
        b_aj = reshape(Ox*b_aj(1:Nx*Ny) + Oy*b_aj(Nx*Ny+1:end),[Nx,Ny]);
        b_aj(isnan(b_aj)) = 0;
        
        AF = extra.AF;
        [fields_aj, ~] = FDFD_fast(ER,MuR,RES,NPML,BC,lambda0,Pol,b_aj,AF);
        x_aj = fields_aj.x/E0;
        Ex_aj = reshape(x_aj(1:Nx*Ny),[Nx,Ny]);
        Ey_aj = reshape(x_aj(Nx*Ny+1:end),[Nx,Ny]);
        
        AVM = -real((Ex.*Ex_aj.*delta_device + Ey.*Ey_aj.*delta_device));
        ER = ER + alpha*AVM + alpha*gamma*AVM_prev;
        AVM_prev = AVM;
        
        ER(ER < 1) = 1;
        ER(ER > eps) = eps;
        
        if abs(g) > G_best_local
            G_best_local = abs(g);
            ER_best = ER;
        end
        Gs(jj) = real(g);
        
        if display_plots && mod(jj,skip)==0
            clf;
            subplot(2,2,1);
            disp_map = [];
            for k_ = 1:5
                disp_map = [disp_map; real(ER)];
            end
            imagesc(disp_map,[1,eps])
            colormap(flipud(gray))
            title('relative permittivity')
            set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal')
            set(gca,'FontSize',22,'fontWeight','normal')
            colorbar()
            
            subplot(2,2,2);
            plot(Gs(1:jj),'k');
            xlabel('iteration number')
            ylabel('gradient (E_0)')
            title('acceleration gradient at \phi = 0')
            set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal')
            set(gca,'FontSize',22,'fontWeight','normal')
            colorbar()
            
            subplot(2,2,3); hold on;
            plot(1:jj, phis(1:jj));
            plot(1:jj, zeros(jj,1));
            xlabel('iteration number');
            ylabel('\phi');
            legend({'computed','\phi=0 (target)'});
            title('acceleration phase (\phi)');
            set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal')
            set(gca,'FontSize',22,'fontWeight','normal')
            pause(0.001);
        end
    end
    
    %% POST PROCESSING STUFF
    eps_avg = (eps+1)/2;
    ER_best(ER_best < eps_avg) = 1;
    ER_best(ER_best >= eps_avg) = eps;
    ER_best_1D{k} = ER_best;
    
    [fields_best, extra_best] = FDFD_TFSF(ER_best,MuR,RES,NPML,BC,lambda0,Pol,b,kinc);
    Ex_best = fields_best.Ex/E0;
    Ey_best = fields_best.Ey/E0;
    g1_best = sum(sum(eta1.*Ex_best));
    g2_best = sum(sum(eta2.*Ex_best));
    g_best  = g1_best + g2_best;
    G_best_local_final = abs(g_best);
    G1_best = abs(g1_best);
    G2_best = abs(g2_best);
    
    E_abs = delta_device .* sqrt(abs(Ex_best).^2 + abs(Ey_best).^2);
    E_max = max(E_abs(:));
    
    %%%% 結果を 1D 配列へ格納
    G_best_values_1D(k)             = G_best_local_final;
    G_best_abs_sums_1D(k)           = (abs(g1_best) + abs(g2_best));
    G_best_values_times_gap_1D(k)   = G_best_local_final * gap_nm;
    G_best_abs_sums_times_gap_1D(k) = (abs(g1_best) + abs(g2_best)) * gap_nm;
    G_best_values_final_1D(k)       = G_best_local_final;
    G1_best_1D(k)                   = G1_best;
    G2_best_1D(k)                   = G2_best;
    g_best_complex_1D(k)            = g_best;
    g1_best_complex_1D(k)           = g1_best;
    g2_best_complex_1D(k)           = g2_best;
    E_max_1D(k)                     = E_max;
    abs_g1_plus_g2_times_gap_1D(k)  = (abs(g1_best) + abs(g2_best)) * gap_nm;
    abs_g_best_times_gap_1D(k)      = abs(g_best) * gap_nm;
    abs_g1_best_1D(k)               = abs(g1_best);
    abs_g2_best_1D(k)               = abs(g2_best);
    abs_sum_g_times_gap_1D(k)       = (abs(g1_best) + abs(g2_best)) * gap_nm;
    
    %%% NEW %%% 電場分布の計算を parfor 内で実施し，結果を保存（最適化後構造での電場分布）
    [fields_best_for_plot, ~] = FDFD_TFSF(ER_best, MuR, [dlx, dlx], [0,0,npml,npml], BC, lambda0, Pol, b, kinc);
    Ex_best_plot = fields_best_for_plot.Ex;
    Ey_best_plot = fields_best_for_plot.Ey;
    E_best_magnitude = sqrt(abs(Ex_best_plot).^2 + abs(Ey_best_plot).^2);
    E_field_magnitude_1D{k} = E_best_magnitude;
    
end % end of parfor

% -------------------------------------------------------------
%  (gap, gap_gap) ごとの結果をファイル書き出し & 画像保存
% -------------------------------------------------------------
% dlx はすべての組み合わせで同じなのでここで再定義
dlx = lambda0/grids_in_lam;

for k = 1:nComb
    idx     = floor((k-1)/ngapgap) + 1;
    jGapGap = mod(k-1, ngapgap) + 1;
    gap_nm     = gap_nm_values(idx);
    gap_gap_nm = gap_gap_nm_values(jGapGap);
    
    %%% 個別テキスト出力
    fname = sprintf('%s/final_acceleration_gradients_gap_%d_gapgap_%d_%s.txt', ...
        output_folder_name, gap_nm, gap_gap_nm, timestamp);
    fileID = fopen(fname, 'w');
    fprintf(fileID, 'G_best (abs): %f\n',   G_best_values_final_1D(k));
    fprintf(fileID, 'G1_best (abs): %f\n',  G1_best_1D(k));
    fprintf(fileID, 'G2_best (abs): %f\n',  G2_best_1D(k));
    fprintf(fileID, 'g_best (complex) = %.4f + %.4fi\n', ...
        real(g_best_complex_1D(k)), imag(g_best_complex_1D(k)));
    fprintf(fileID, 'g1_best (complex) = %.4f + %.4fi\n', ...
        real(g1_best_complex_1D(k)), imag(g1_best_complex_1D(k)));
    fprintf(fileID, 'g2_best (complex) = %.4f + %.4fi\n', ...
        real(g2_best_complex_1D(k)), imag(g2_best_complex_1D(k)));
    fprintf(fileID, 'E_max: %f\n',          E_max_1D(k));
    fprintf(fileID, '(abs(g1)+abs(g2))*gap: %f\n', abs_g1_plus_g2_times_gap_1D(k));
    fprintf(fileID, 'abs(g1+g2)*gap: %f\n', abs_g_best_times_gap_1D(k));
    fclose(fileID);
    fprintf('File saved as: %s\n', fname);
    
    %%% 構造 (ER_best) の 2D イラスト出力（軸ラベルを "pixel" に設定） - 5回繰り返し
    ER_best_k = ER_best_1D{k};
    bestFig = figure('Name','Best Structure','Visible','off');
    disp_best = [];
    for kk_ = 1:5
        disp_best = [disp_best; real(ER_best_k)]; %#ok<AGROW>
    end
    imagesc(disp_best, [1, eps]);
    colormap(flipud(gray));
    axis equal tight;
    colorbar();
    xlabel('pixel');
    ylabel('pixel');
    title(sprintf('Best Structure (gap = %d nm, gap\\_gap = %d nm)', gap_nm, gap_gap_nm));
    
    figNameBest = sprintf('%s/best_structure_gap_%d_gapgap_%d_%s.png', ...
        output_folder_name, gap_nm, gap_gap_nm, timestamp);
    saveas(bestFig, figNameBest);
    close(bestFig);
    
    %%% NEW %%% 電場分布の出力（5回繰り返し）
    E_best_magnitude = E_field_magnitude_1D{k};
    efieldFig = figure('Name','Electric Field Magnitude','Visible','off');
    disp_efield = [];
    for kk_ = 1:5
        disp_efield = [disp_efield; E_best_magnitude];  %#ok<AGROW>
    end
    imagesc(disp_efield);
    axis equal tight;
    colormap jet;
    colorbar();
    xlabel('pixel');
    ylabel('pixel');
    title(sprintf('Electric Field Magnitude (gap = %d nm, gap\\_gap = %d nm)', gap_nm, gap_gap_nm));
    
    efieldFigName = sprintf('%s/electric_field_gap_%d_gapgap_%d_%s.png', ...
        output_folder_name, gap_nm, gap_gap_nm, timestamp);
    saveas(efieldFig, efieldFigName);
    close(efieldFig);
end

%%% NEW %%%
% -------------------------------------------------------------
%  全組み合わせをまとめた CSV ファイルを出力
% -------------------------------------------------------------
summary_filename = sprintf('%s/summary_all_%s.csv', output_folder_name, timestamp);
fid_summary = fopen(summary_filename, 'w');
% ヘッダ行 (カンマ区切り)
fprintf(fid_summary, 'gap_nm,gap_gap_nm,G_best_final,G1_best,G2_best,real_g_best,imag_g_best,E_max,(abs(g1)+abs(g2))*gap,abs(g1+g2)*gap\n');

for k = 1:nComb
    idx     = floor((k-1)/ngapgap) + 1;
    jGapGap = mod(k-1, ngapgap) + 1;
    gap_nm     = gap_nm_values(idx);
    gap_gap_nm = gap_gap_nm_values(jGapGap);
    
    fprintf(fid_summary, '%d,%d,%f,%f,%f,%f,%f,%f,%f,%f\n', ...
        gap_nm, gap_gap_nm, ...
        G_best_values_final_1D(k), ...
        G1_best_1D(k), ...
        G2_best_1D(k), ...
        real(g_best_complex_1D(k)), ...
        imag(g_best_complex_1D(k)), ...
        E_max_1D(k), ...
        abs_g1_plus_g2_times_gap_1D(k), ...
        abs_g_best_times_gap_1D(k) );
end
fclose(fid_summary);
fprintf('Summary file saved as: %s\n', summary_filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (新規) 各種プロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. abs(g1) vs gap for each gap_gap
figure('Name','abs(g1) vs gap for each gap_gap');
hold on;
for jGapGap = 1:ngapgap
    k_vec = (0:ngap-1)*ngapgap + jGapGap;
    plot(gap_nm_values, abs_g1_best_1D(k_vec), '-o', ...
        'DisplayName', sprintf('gap\\_gap = %d nm', gap_gap_nm_values(jGapGap)));
end
legend('show');
xlabel('gap (nm)');
ylabel('abs(g1)');
title('abs(g1) vs gap for each gap\_gap');
grid on;
saveas(gcf, sprintf('%s/abs_g1_vs_gap_for_each_gapgap_%s.png', output_folder_name, timestamp));

% 2. abs(g2) vs gap for each gap_gap
figure('Name','abs(g2) vs gap for each gap_gap');
hold on;
for jGapGap = 1:ngapgap
    k_vec = (0:ngap-1)*ngapgap + jGapGap;
    plot(gap_nm_values, abs_g2_best_1D(k_vec), '-o', ...
        'DisplayName', sprintf('gap\\_gap = %d nm', gap_gap_nm_values(jGapGap)));
end
legend('show');
xlabel('gap (nm)');
ylabel('abs(g2)');
title('abs(g2) vs gap for each gap\_gap');
grid on;
saveas(gcf, sprintf('%s/abs_g2_vs_gap_for_each_gapgap_%s.png', output_folder_name, timestamp));

% 3. (abs(g1)+abs(g2))*gap vs gap for each gap_gap
figure('Name','(abs(g1)+abs(g2))*gap vs gap for each gap_gap');
hold on;
for jGapGap = 1:ngapgap
    k_vec = (0:ngap-1)*ngapgap + jGapGap;
    plot(gap_nm_values, abs_sum_g_times_gap_1D(k_vec), '-o', ...
        'DisplayName', sprintf('gap\\_gap = %d nm', gap_gap_nm_values(jGapGap)));
end
legend('show');
xlabel('gap (nm)');
ylabel('(abs(g1)+abs(g2))*gap');
title('(abs(g1)+abs(g2))*gap vs gap for each gap\_gap');
grid on;
saveas(gcf, sprintf('%s/abs_sum_g_times_gap_vs_gap_for_each_gapgap_%s.png', output_folder_name, timestamp));

% 4. (abs(g1)+abs(g2) with best gap_gap) * gap vs gap
G_abs_sums_best_for_each_gap = zeros(ngap,1);
idx_best_for_each_gap = zeros(ngap,1);
abs_sum_g_best_times_gap_for_each_gap = zeros(ngap,1);

for iGap = 1:ngap
    k_vec = (iGap-1)*ngapgap + (1:ngapgap);
    [G_abs_sums_best_for_each_gap(iGap), localBestIdx] = max(G_best_abs_sums_1D(k_vec));
    idx_best_for_each_gap(iGap) = localBestIdx;
    abs_sum_g_best_times_gap_for_each_gap(iGap) = max(abs_sum_g_times_gap_1D(k_vec));
end
best_gapgap_for_each_gap = gap_gap_nm_values(idx_best_for_each_gap);

figure('Name','(abs(g1)+abs(g2) with best gap_gap)*gap vs gap');
plot(gap_nm_values, abs_sum_g_best_times_gap_for_each_gap, '-o');
xlabel('gap (nm)');
ylabel('(abs(g1)+abs(g2))_{best gap\_gap} * gap');
title('(abs(g1)+abs(g2)) with best gap\_gap * gap vs gap');
grid on;
saveas(gcf, sprintf('%s/abs_sum_g_best_times_gap_vs_gap_%s.png', output_folder_name, timestamp));

% 各 gap_gap を凡例として，gap vs (abs(g1)+abs(g2)) を1次元プロット
figure('Name','abs(g1)+abs(g2) vs gap for each gap_gap');
hold on;
for jGapGap = 1:ngapgap
    k_vec = (0:ngap-1)*ngapgap + jGapGap;
    plot(gap_nm_values, G_best_abs_sums_1D(k_vec), '-o', ...
        'DisplayName', sprintf('gap\\_gap = %d nm', gap_gap_nm_values(jGapGap)));
end
legend('show');
xlabel('gap (nm)');
ylabel('abs(g1)+abs(g2)');
title('abs(g1)+abs(g2) vs gap for each gap\_gap');
grid on;
saveas(gcf, sprintf('%s/abs_g1_plus_abs_g2_vs_gap_for_each_gapgap_%s.png', output_folder_name, timestamp));

% 各 gap で最大となる (abs(g1)+abs(g2)) を抽出して1次元プロット
figure('Name','abs(g1)+abs(g2) with best gap_gap vs gap');
plot(gap_nm_values, G_abs_sums_best_for_each_gap, '-o');
xlabel('gap (nm)');
ylabel('abs(g1)+abs(g2) with best gap\_gap');
title('abs(g1)+abs(g2) with best gap\_gap vs gap');
grid on;
saveas(gcf, sprintf('%s/abs_g1_plus_abs_g2_best_vs_gap_%s.png', output_folder_name, timestamp));

delete(gcp('nocreate'));
