addpath(genpath('./'));                     % add the whole directory to path, if not already done

%% SET PARAMETERS
c0 = 1;                                     % speed of light m/s (normalized to 1)
lambda0 = 2;                                % central wavelength (um)

N = 4000;                                   % number of iterations
display_plots = true;                       % プロット＆動画を作るなら true に

skip = 30;                                  % 30 ループに1回描画・動画に書き込み
alpha = 5e2;
a = 3;
beta = 0.5;
in_material = false;
starting = 0;
grids_in_lam = 100;
gap_nm_values = 360;                        % gap size in nm
timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
L = 0.4;
npml = 10;
eps = 3.4363^2;                             % 例: Si (2um)
nmax = sqrt(eps);
gamma = 0.9;
output_folder_name = 'result/single_channel_step_40_jan23_parallel_grids_100_L04';

dlx = lambda0/grids_in_lam;
dly = dlx;
G_best_values = zeros(length(gap_nm_values), 1);
G_best_times_gap_values = zeros(length(gap_nm_values), 1);

for idx = 1:length(gap_nm_values)
    gap_nm = gap_nm_values(idx);
    
    pos_src = floor(npml + grids_in_lam/4);
    spc_pts = floor(grids_in_lam/4);
    gap_pts = floor(gap_nm/1000/dlx);
    Lpts    = round(L/dlx);
    
    Nx = ceil(lambda0*beta/dlx);
    Ny = gap_pts + 2*(pos_src + Lpts + spc_pts);
    
    nx = floor(Nx/2);
    ny = floor(Ny/2);
    
    for min_G_Emax = 0
        ER  = ones(Nx,Ny);
        MuR = ones(Nx,Ny);
        ER_best = ones(Nx,Ny);
        
        b = zeros(Nx,Ny);
        b(:, pos_src:pos_src + spc_pts + Lpts + gap_pts + Lpts + spc_pts) = 1;
        kinc = [0,1];
        RES = [dlx,dly];
        BC = [-1,-1];
        NPML = [0,0,npml,npml];
        Pol= 'Hz';
        spc = spc_pts*dly;
        gap = gap_pts*dly;
        
        delta_device = zeros(Nx,Ny);
        delta_device(1:Nx, pos_src + spc_pts : pos_src + spc_pts + Lpts) = 1;
        delta_device(1:Nx, ...
            pos_src + spc_pts + Lpts + gap_pts : pos_src + spc_pts + Lpts + gap_pts + Lpts) = 1;
        
        eta = zeros(Nx,Ny);
        eta(:, ny) = 1/Nx * exp(2*pi*1i*dlx*(0:Nx-1)/lambda0/beta);
        eta_vec = eta(:);
        
        % 初期値設定
        for i = 1:Nx
            for j = 1:Ny
                if delta_device(i,j) == 1
                    if starting == 1
                        ER(i,j) = rand * (eps - 1) + 1;
                    elseif starting == 2
                        ER(i,j) = eps/2 + 0.5;
                    else
                        % starting == 0 なら 1 のまま
                    end
                end
            end
        end
        
        % 空間全体真空で FDFD を回して E0 を取得
        [fields, ~] = FDFD_TFSF(ones(Nx,Ny), MuR, RES, NPML, BC, lambda0, Pol, b, kinc);
        Ex = fields.Ex;
        Ey = fields.Ey;
        E0 = sqrt(abs(Ex(nx, ny))^2 + abs(Ey(nx, ny))^2);
        
        G_best = 0;
        Gs     = zeros(N,1);
        E_maxs = zeros(N,1);
        G_by_Es = zeros(N,1);
        G_by_Sa = zeros(N,1);
        phis   = zeros(N,1);
        AVM_prev = zeros(Nx,Ny);
        
        %＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
        % 動画保存 (display_plots == true) → MP4 出力
        if display_plots
            %--- 図1: 構造の 2D プロット用 ---
            figure_handle_1 = figure(1);
            videoFileName_1 = sprintf('%s/video_structure_gap_%d_%s.mp4',...
                output_folder_name, gap_nm, timestamp);
            v1 = VideoWriter(videoFileName_1,'MPEG-4');
            v1.FrameRate = 5;
            v1.Quality   = 95;
            open(v1);
            
            %--- 図2: gradient グラフ用 ---
            figure_handle_2 = figure(2);
            videoFileName_2 = sprintf('%s/video_gradient_gap_%d_%s.mp4',...
                output_folder_name, gap_nm, timestamp);
            v2 = VideoWriter(videoFileName_2,'MPEG-4');
            v2.FrameRate = 5;
            v2.Quality   = 95;
            open(v2);
        end
        %＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
        
        for j = 1:N
            [fields, extra] = FDFD_TFSF(ER, MuR, RES, NPML, BC, lambda0, Pol, b, kinc);
            Ex = fields.Ex / E0;
            Ey = fields.Ey / E0;
            g = sum(sum(eta.*Ex));
            G = real(g);
            phis(j) = angle(g);
            
            DEY = extra.derivatives.DEY;
            DEX = extra.derivatives.DEX;
            ER_vec = ER(:);
            chi = delta_device .* (ER - 1);
            
            Ox = -1i * lambda0/(2*pi*c0) * spdiags(1./ER_vec, 0, Nx*Ny, Nx*Ny) * DEY;
            Oy =  1i * lambda0/(2*pi*c0) * spdiags(1./ER_vec, 0, Nx*Ny, Nx*Ny) * DEX;
            
            eta_aj = [eta_vec; zeros(Nx*Ny,1)];
            
            if in_material
                E_abs = (chi/(eps-1)) .* sqrt(abs(Ex).^2 + abs(Ey).^2);
            else
                E_abs = delta_device .* sqrt(abs(Ex).^2 + abs(Ey).^2);
            end
            
            x_abs = E_abs(:);
            alpha_vec = exp(x_abs*a);
            alpha_T_1 = sum(alpha_vec);
            Sa = sum(alpha_vec.*x_abs) / alpha_T_1;
            
            x = [Ex(:); Ey(:)];
            z = conj(x ./ [x_abs; x_abs]);
            z(isnan(z)) = 0;
            z(isinf(z)) = 0;
            
            spdiagz = spdiags(z, 0, Nx*Ny*2, Nx*Ny*2);
            P = [speye(Nx*Ny), speye(Nx*Ny)];
            R = P * spdiagz;
            
            S = real(1/alpha_T_1*(speye(Nx*Ny) + a*spdiags(x_abs,0,Nx*Ny,Nx*Ny) ...
                - a*sum(alpha_vec.*x_abs)/alpha_T_1*speye(Nx*Ny)));
            sigma = transpose(alpha_vec)*S*R;
            sigma(isnan(sigma)) = 0;
            
            b_aj1 = transpose(G/Sa^2 * sigma);
            b_aj2 = -eta_aj / Sa;
            
            if min_G_Emax
                b_aj = b_aj1 + b_aj2;
            else
                b_aj = -eta_aj;
            end
            
            b_aj = reshape(Ox*b_aj(1:Nx*Ny) + Oy*b_aj(Nx*Ny+1:end), [Nx,Ny]);
            b_aj(isnan(b_aj)) = 0;
            
            AF = extra.AF;
            [fields_aj, ~] = FDFD_fast(ER,MuR,RES,NPML,BC,lambda0,Pol,b_aj,AF);
            x_aj = fields_aj.x / E0;
            Ex_aj = reshape(x_aj(1:Nx*Ny), [Nx,Ny]);
            Ey_aj = reshape(x_aj(Nx*Ny+1:end), [Nx,Ny]);
            
            AVM = -real(Ex.*Ex_aj.*delta_device + Ey.*Ey_aj.*delta_device);
            
            E_abs_val = max(E_abs(:));
            E_maxs(j) = E_abs_val;
            Gs(j) = G;
            G_by_Es(j) = G / E_abs_val;
            G_by_Sa(j) = G / Sa;
            
            ER = ER + alpha*AVM + alpha*gamma*AVM_prev;
            AVM_prev = AVM;
            ER(ER < 1)   = 1;
            ER(ER > eps) = eps;
            
            if G > G_best
                G_best = G;
                ER_best = ER;
            end
            
            %＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            % skip ループごとに 2つのウィンドウをそれぞれ描画 → 動画に書き込み
            if display_plots && mod(j, skip) == 0
                %--- (1) 図1: 構造の2Dプロットを出力 ---
                % 図1: relative permittivity の描画
                figure(1);
                clf;
                disp_mat = [];
                for kFrame = 1:5
                    disp_mat = [disp_mat; real(ER)];
                end
                imagesc(disp_mat, [1, eps]);
                colormap(flipud(gray));
                title('relative permittivity');
                xlabel('pixel');   % ← x 軸に "pixel"
                ylabel('pixel');   % ← y 軸に "pixel"
                colorbar();
                set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal');
                set(gca,'FontSize',22,'fontWeight','normal');
                drawnow;
                
                frame1 = getframe(gcf);
                writeVideo(v1, frame1);
                
                %--- (2) 図2: 勾配 G の履歴を出力 ---
                figure(2);  % figure_handle_2 をアクティブに
                clf;
                plot(Gs(1:j),'k');
                xlabel('iteration number')
                ylabel('G')  % "gradient (E0)" から "G" に変更
                title('acceleration gradient G')
                set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','normal')
                set(gca,'FontSize',22,'fontWeight','normal')
                grid on;
                drawnow;
                
                frame2 = getframe(gcf);
                writeVideo(v2, frame2);
            end
            %＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
        end
        
        %＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
        % 動画ファイルを閉じる (2つ分)
        if display_plots
            close(v1);
            close(v2);
        end
        %＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
        
        %% POST PROCESSING
        eps_avg = (eps + 1)/2;
        ER(ER < eps_avg)    = 1;
        ER(ER >= eps_avg)   = eps;
        ER_best(ER_best < eps_avg)   = 1;
        ER_best(ER_best >= eps_avg)  = eps;
        
        [fields_best, extra_best] = FDFD_TFSF(ER_best, MuR, RES, NPML, BC, lambda0, Pol, b, kinc);
        Ex_best = fields_best.Ex / E0;
        Ey_best = fields_best.Ey / E0;
        
        g_best = sum(sum(eta .* Ex_best));
        G_best = abs(g_best);
        
        E_abs_best = delta_device .* sqrt(abs(Ex_best).^2 + abs(Ey_best).^2);
        E_max_best = max(E_abs_best(:));
        
        fname = sprintf('%s/final_acceleration_gradients_gap_%d_%s_gap_%d.txt', ...
            output_folder_name, gap_nm, timestamp, gap_nm);
        fileID = fopen(fname, 'w');
        fprintf(fileID, 'G_best: %f\n', G_best);
        fprintf(fileID, 'g_best: %f + %fi\n', real(g_best), imag(g_best));
        fprintf(fileID, 'E_max: %f\n', E_max_best);
        fprintf(fileID, 'L: %f\n', L);
        fclose(fileID);
        
        G_best_values(idx) = G_best;
        G_best_times_gap_values(idx) = G_best * gap_nm;
    end
end

% ... (以降はお好みで最終結果の処理/プロットなど)
