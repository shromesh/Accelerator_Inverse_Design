addpath(genpath('./'));                     % add the whole directory to path, if not already done

%% SET PARAMETERS
c0 = 1;                                     % speed of light m/s (normalized to 1)
lambda0 = 2;                                % central wavelength (um)

skip = 4;                                   % number of iteration frames between plots (higher->faster, lower->more plots)
display_plots = false;                      % plotting during the run? (false にするとiteration中の表示を行わない)

alpha = 5e2;                                % step size in permittivity (~1e2-1e4 works well)
a = 3;                                      % smooth-max weight factor (see paper)
beta = 0.5;                                 % ratio of electron speed to speed of light
N = 2000;                                    % number of iterations

in_material = false;                        % evaluate E_max in material? or in surrounding regions.
starting = 0;                               % 0 -> vacuum, 1 -> random, 2 -> midway epsilon

grids_in_lam = 100;                         % number of grid points in a free space wavelength
npml = 10;                                  % number of PML (absorbing region) points (need > 10 at least)

% relative permittivity of material region.  uncomment to select
eps = 3.4363^2;     % Si 2um
% eps = 1.4381^2;    % fused silica 2um
% eps = 1.9834^2;    % Si3N4
% eps = 1.9^2;       % GaOx

gamma = 0.9;                                % 'momentum term', see paper. 0-1

%% 新たに追加: gap を変化させるための配列
gap_nm_values = 200:20:1000;
% gap_nm_values = [200, 300];

%% gap_gap を変化させるための配列
gap_gap_nm_values = 300:200:1000;
% gap_gap_nm_values = [300, 400];

%% 出力フォルダ名を設定
output_folder_name = 'result/double_channel_gap_step_20_gapgap_step_200_Jan19';

% -------------------------------------------------------------
% 2D で結果を保持するために，配列の長さを取得
ngap   = length(gap_nm_values);
ngapgap = length(gap_gap_nm_values);

% G_best_valuesなどを 2次元配列化
G_best_values_2D             = zeros(ngap, ngapgap);
G_best_abs_sums_2D           = zeros(ngap, ngapgap);
G_best_values_times_gap_2D   = zeros(ngap, ngapgap);
G_best_abs_sums_times_gap_2D = zeros(ngap, ngapgap);

% -------------------------------------------------------------
% ループ開始
iGap = 0;
for gap_nm = gap_nm_values
    iGap = iGap + 1;
    jGapGap = 0;
    for gap_gap_nm = gap_gap_nm_values
        jGapGap = jGapGap + 1;
        
        %% SET OTHER CONSTANTS (DON'T CHANGE)
        dlx = lambda0/grids_in_lam;                 % grid size along electron trajectory axis
        dly = dlx;                                  % spacing in the perpendicular direction
        
        % gap_nm から grid point に換算
        gap_pts = floor(gap_nm/1000/dlx);           % number of grid points in the gap
        
        % gap_gap_nm から grid point に換算
        gap_gap_pts = floor(gap_gap_nm/1000/dlx);   % number of grid points in the gap between the two gaps
        
        L = 0.4;                                    % size of optimization region (um)
        Lpts = round(L/dlx);                        % number of points in the optimization region
        
        pos_src = floor(npml+grids_in_lam/4);       % number of grid points between left edge and source
        spc_pts = floor(grids_in_lam/4);            % number of grid points between source and structure
        
        Nx = ceil(lambda0*beta/dlx);
        % 2つのギャップ + 中央 gap_gap_pts + 上下2つの最適化領域 + PML の外の領域 など
        Ny = 2*gap_pts + 2*(pos_src + Lpts + spc_pts) + gap_gap_pts;
        
        nx = floor(Nx/2);
        ny1 = floor(gap_pts/2 + pos_src + Lpts + spc_pts);
        ny2 = floor(gap_pts + gap_pts/2 + gap_gap_pts + pos_src + Lpts + spc_pts);
        
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
        
        % define starting permittivity
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
        AVM_prev = zeros(Nx,Ny);
        
        if display_plots
            figure(1);  % open a figure to plot
        end
        
        display('working on gradient maximized structure');
        upd = textprogressbar(N);
        
        Gs     = zeros(N,1);
        E_maxs = zeros(N,1);
        phis   = zeros(N,1);
        
        for jj = (1:N)
            
            upd(jj);
            % original simulation
            [fields, extra] = FDFD_TFSF(ER,MuR,RES,NPML,BC,lambda0,Pol,b,kinc);
            Ex = fields.Ex/E0;
            Ey = fields.Ey/E0;
            
            % compute gradients
            g1 = sum(sum(eta1.*Ex));
            g2 = sum(sum(eta2.*Ex));
            g  = g1 + g2;
            G  = real(g);
            
            % get phase
            phis(jj) = angle(g);
            
            % get numerical spatial derivative operators
            DEY = extra.derivatives.DEY;
            DEX = extra.derivatives.DEX;
            
            ER_vec = ER(:);
            chi = delta_device.*(ER - ones(Nx,Ny));
            
            % in_material を考慮した E_abs
            if (in_material)
                E_abs = (chi/(eps-1)).*sqrt(abs(Ex).^2 + abs(Ey).^2);
            else
                E_abs = delta_device.*sqrt(abs(Ex).^2 + abs(Ey).^2);
            end
            
            E_abs_vec = E_abs(:);
            E_maxs(jj) = max(E_abs_vec);
            
            % アジュゲート場計算に必要な諸々
            Ox = -1i*lambda0/(2*pi*c0)*spdiags(1./ER_vec,0,Nx*Ny,Nx*Ny)*DEY;
            Oy =  1i*lambda0/(2*pi*c0)*spdiags(1./ER_vec,0,Nx*Ny,Nx*Ny)*DEX;
            
            eta1_aj = [eta1_vec; zeros(Nx*Ny,1)];
            eta2_aj = [eta2_vec; zeros(Nx*Ny,1)];
            
            % ここでは簡略的に、従来のアルゴリズム通りに
            b_aj = - (eta1_aj + eta2_aj);
            b_aj = reshape(Ox*b_aj(1:Nx*Ny) + Oy*b_aj(Nx*Ny+1:end),[Nx,Ny]);
            b_aj(isnan(b_aj)) = 0 ;
            
            AF = extra.AF;
            [fields_aj, ~] = FDFD_fast(ER,MuR,RES,NPML,BC,lambda0,Pol,b_aj,AF);
            x_aj = fields_aj.x/E0;
            Ex_aj = reshape(x_aj(1:Nx*Ny),[Nx,Ny]);
            Ey_aj = reshape(x_aj(Nx*Ny+1:end),[Nx,Ny]);
            
            AVM = -real((Ex.*Ex_aj + Ey.*Ey_aj).*delta_device);
            
            % update permittivity
            ER = ER + alpha*AVM + alpha*gamma*AVM_prev;
            AVM_prev = AVM;
            
            % permittivity の上下限クリップ
            ER(ER < 1) = 1;
            ER(ER > eps) = eps;
            
            % ベスト更新
            if (abs(g) > G_best_local)
                G_best_local = abs(g);
                ER_best = ER;
            end
            
            Gs(jj) = real(g);
            
            % ---- プロット (iteration中)
            if display_plots && mod(jj,skip)==0
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
                plot(Gs(1:jj),'k');
                xlabel('iteration number')
                ylabel('power (G)')
                title('acceleration gradient at \phi = 0')
                set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal')
                set(gca,'FontSize',22,'fontWeight','normal')
                colorbar()
                
                subplot(2,2,3); hold all;
                plot((1:jj), phis(1:jj));
                plot((1:jj), zeros(jj,1));
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
        ER_best(ER_best<eps_avg) = 1;
        ER_best(ER_best>=eps_avg) = eps;
        
        % do another simulation of the binary distribution for ER_best
        [fields_best, extra_best] = FDFD_TFSF(ER_best,MuR,RES,NPML,BC,lambda0,Pol,b,kinc);
        Ex_best = fields_best.Ex/E0;
        Ey_best = fields_best.Ey/E0;
        
        % calculate g1_best and g2_best after the loop
        g1_best = sum(sum(eta1.*Ex_best));
        g2_best = sum(sum(eta2.*Ex_best));
        g_best  = g1_best + g2_best;
        
        G_best_local_final = abs(g_best);
        G1_best = abs(g1_best);
        G2_best = abs(g2_best);
        
        % E_max の計算 (binary最終構造で)
        E_abs = delta_device.*sqrt(abs(Ex_best).^2 + abs(Ey_best).^2);
        E_max = max(E_abs(:));
        
        % テキスト出力
        timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
        fname = sprintf('%s/final_acceleration_gradients_gap_%d_gapgap_%d_%s.txt', ...
            output_folder_name, gap_nm, gap_gap_nm, timestamp);
        fileID = fopen(fname, 'w');
        fprintf(fileID, 'G_best (abs): %f\n', G_best_local_final);
        fprintf(fileID, 'G1_best (abs): %f\n', G1_best);
        fprintf(fileID, 'G2_best (abs): %f\n', G2_best);
        fprintf(fileID, 'g_best (complex) = %.4f + %.4fi\n', real(g_best), imag(g_best));
        fprintf(fileID, 'g1_best (complex) = %.4f + %.4fi\n', real(g1_best), imag(g1_best));
        fprintf(fileID, 'g2_best (complex) = %.4f + %.4fi\n', real(g2_best), imag(g2_best));
        fprintf(fileID, 'E_max: %f\n', E_max);
        fprintf(fileID, '(abs(g1)+abs(g2))*gap: %f\n', (abs(g1_best) + abs(g2_best)) * gap_nm);
        fprintf(fileID, 'abs(g1+g2)*gap: %f\n', abs(g_best) * gap_nm);
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
        title(sprintf('Best Structure (gap = %d nm, gap\\_gap = %d nm)', gap_nm, gap_gap_nm));
        colorbar();
        
        figNameBest = sprintf('%s/best_structure_gap_%d_gapgap_%d_%s.png', ...
            output_folder_name, gap_nm, gap_gap_nm, timestamp);
        saveas(bestFig, figNameBest);
        
        % ----------------------------------
        % 結果を2次元配列に格納
        G_best_values_2D(iGap, jGapGap)             = G_best_local_final;
        G_best_abs_sums_2D(iGap, jGapGap)           = (abs(g1_best) + abs(g2_best));
        G_best_values_times_gap_2D(iGap, jGapGap)   = G_best_local_final * gap_nm;
        G_best_abs_sums_times_gap_2D(iGap, jGapGap) = (abs(g1_best) + abs(g2_best)) * gap_nm;
        
    end % end of gap_gap_nm loop
    
end % end of gap_nm loop

% -------------------------------------------------------------
% (新) ここで「gap_nm vs gap_gap_nm」で色としてG_bestなどを表示 (2次元ヒートマップ)
% imagesc でも scatter でもOK。ここでは imagesc 例を示す。
% gap_nm_values (縦軸) と gap_gap_nm_values (横軸) を使う。
% imagesc は (x,y,Z) の順序に注意し，axis xy で上が大きい方にする。

% 1. abs(g1+g2)
figure('Name','abs(g1+g2)');
imagesc(gap_gap_nm_values, gap_nm_values, G_best_values_2D);
set(gca, 'YDir', 'normal');  % 上下反転を防ぐ
colorbar();
xlabel('gap\_gap (nm)');
ylabel('gap (nm)');
title('G\_best = abs(g1+g2)');
saveas(gcf, sprintf('%s/abs_g1_plus_g2_%s.png', output_folder_name, timestamp));

% 2. abs(g1) + abs(g2)
figure('Name','abs(g1) + abs(g2)');
imagesc(gap_gap_nm_values, gap_nm_values, G_best_abs_sums_2D);
set(gca, 'YDir', 'normal');
colorbar();
xlabel('gap\_gap (nm)');
ylabel('gap (nm)');
title('abs(g1)+abs(g2)');
saveas(gcf, sprintf('%s/abs_g1_plus_abs_g2_%s.png', output_folder_name, timestamp));

% 3. abs(g1+g2)*gap_nm
figure('Name','abs(g1+g2)*gap');
imagesc(gap_gap_nm_values, gap_nm_values, G_best_values_times_gap_2D);
set(gca, 'YDir', 'normal');
colorbar();
xlabel('gap\_gap (nm)');
ylabel('gap (nm)');
title('abs(g1+g2)*gap');
saveas(gcf, sprintf('%s/abs_g1_plus_g2_times_gap_%s.png', output_folder_name, timestamp));

% 4. (abs(g1)+abs(g2))*gap_nm
figure('Name','(abs(g1)+abs(g2))*gap');
imagesc(gap_gap_nm_values, gap_nm_values, G_best_abs_sums_times_gap_2D);
set(gca, 'YDir', 'normal');
colorbar();
xlabel('gap\_gap (nm)');
ylabel('gap (nm)');
title('(abs(g1)+abs(g2))*gap');
saveas(gcf, sprintf('%s/abs_g1_plus_abs_g2_times_gap_%s.png', output_folder_name, timestamp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (新規) gap_gapをlegendとして、gap vs G=abs(g1)+abs(g2)を1次元プロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','gap vs abs(g1)+abs(g2) for each gap_gap');
hold on;
for jGapGap = 1:ngapgap
    plot(gap_nm_values, G_best_abs_sums_2D(:, jGapGap), '-o', ...
        'DisplayName', sprintf('gap\\_gap = %d nm', gap_gap_nm_values(jGapGap)));
end
legend('show');  % 凡例を表示
xlabel('Gap size (nm)');
ylabel('abs(g1)+abs(g2)');
title('abs(g1)+abs(g2) vs gap for each gap\_gap');
grid on;

% 結果を保存（例: PNG 形式）
saveas(gcf, sprintf('%s/abs_g1_plus_abs_g2_vs_gap_for_each_gapgap_%s.png', ...
    output_folder_name, timestamp));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (新規) 各 gap で最大となる (abs(g1)+abs(g2)) を抽出して1次元プロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- 各gapについて、gap_gapを変化させた中で最もG = abs(g1)+abs(g2)が大きい値を取り出す
%    G_best_abs_sums_2D(iGap, jGapGap) = (|g1| + |g2|) の2次元配列
[G_abs_sums_best_for_each_gap, idx_best_for_each_gap] = max(G_best_abs_sums_2D, [], 2);

% 参考: どの gap_gap で最大になったか知りたい場合は
best_gapgap_for_each_gap = gap_gap_nm_values(idx_best_for_each_gap);

% -- gap vs (最も大きい abs(g1)+abs(g2)) を1次元プロット
figure('Name','Max of abs(g1)+abs(g2) vs gap');
plot(gap_nm_values, G_abs_sums_best_for_each_gap, '-o');
xlabel('Gap size (nm)');
ylabel('max_{gap\\_gap}( abs(g1)+abs(g2) )');
title('Max of abs(g1)+abs(g2) vs gap');
grid on;

% 結果を保存（例: PNG 形式）
saveas(gcf, sprintf('%s/abs_g1_plus_abs_g2_best_vs_gap_%s.png', ...
    output_folder_name, timestamp));
