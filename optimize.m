addpath(genpath('./'));  % add the whole directory to path

%% (A) パラメータ設定
c0 = 1;                  % (規格化) 光速
lambda0 = 2;            % 中心波長 (um)

skip = 4;               % 何iterationごとにプロットを更新するか
display_plots = true;   % iteration中に可視化を行うか

alpha = 5e2;            % 感度分布に対する誘電率更新のステップサイズ
a = 3;                  % smooth-max の重み
beta = 0.5;             % 電子速度 (v) / 光速度 (c0)
N = 1400;               % 反復回数 (iteration)

in_material = false;    % E_max を材料内のみで評価するかどうか
starting = 0;           % 初期の誘電率分布: 0->真空, 1->乱数, 2->中間値(eps/2+0.5)

grids_in_lam = 100;     % 1波長あたり何グリッド置くか
npml = 10;              % PMLセル数 (端の吸収境界条件)

% 材料の相対誘電率 eps
eps = (3.4363)^2;       % 例: Si(2um帯) のn=3.4363程度
%eps = 1.4381^2;        % fused silica などに切り替え可
%eps = 1.9834^2;        % Si3N4
%eps = 1.9^2;           % GaOx

gamma = 0.9;            % モメンタム項 (0〜1で設定)

%% -------------------------------
% ここから「チャネル周りの幾何設定」関連
% -------------------------------
%
% [設定例]
%   - n_channels 個のチャネル
%   - 各チャネルに対して gap_nm (nm) の幅
%   - チャネル間に gap_gap_nm (nm) の幅
%   - 最適化領域 (L) をチャネル数 + 1 個ブロックとして並べる
%
n_channels = 2;         % 例: チャネル数を自由に変更
gap_nm      = 200;      % 各チャネルの gap 幅 (nm)
gap_gap_nm  = 300;      % チャネル間ギャップ (nm)  (チャネルが2つ以上の場合に使用)
L = 0.4;                % 各最適化領域の高さ (um)

% [注意] gap_nm_values や n_channels_values のような配列を用いて
%  複数パラメータをループする場合は，2チャネルコードでやっていたように
%  for ループを作ってください．ここでは簡単に1つの gap_nm だけを例示します．
gap_nm_values = [gap_nm];  % 複数試すなら [200, 400, 600, ...] など

%% 結果保存用フォルダ名
output_folder_name = 'result/n_channels_example';
if ~exist(output_folder_name, 'dir')
    mkdir(output_folder_name);
end

%% G_best 等を記録する変数
G_best_values = [];

%% ----- gap_nm_values などをループする場合の例 -----
for gap_nm = gap_nm_values
    
    % グリッドサイズ (dx, dy)
    dlx = lambda0 / grids_in_lam;
    dly = dlx;
    
    % gap_nm, gap_gap_nm からグリッド数に変換
    gap_pts     = floor(gap_nm     /1000/dlx);
    gap_gap_pts = floor(gap_gap_nm /1000/dlx);
    
    % 最適化領域1つあたりのグリッド数 (Lpts)
    Lpts = round(L/dlx);
    
    % 2チャネルコードに合わせて設定していたもの
    pos_src = floor(npml + grids_in_lam/4);  % ソースから左端までの距離
    spc_pts = floor(grids_in_lam/4);         % ソースから構造までの空き
    
    % ----- 全体の Nx, Ny を定義 -----
    % x方向は変わらず: Nx = ceil(lambda0*beta/dlx)
    Nx = ceil(lambda0 * beta / dlx);
    
    % y方向 (Ny) は
    %   - (pos_src + Lpts + spc_pts) というブロックが上下にある (2倍)
    %   - n_channels 個の gap_pts
    %   - (n_channels - 1) 個の gap_gap_pts
    % を足し合わせます (2チャネルの時にやっていた計算の一般化)．
    Ny = 2*(pos_src + Lpts + spc_pts) ...
        + n_channels       * gap_pts ...
        + (n_channels - 1) * gap_gap_pts;
    
    % 参考用に中心付近の x座標, y座標を取得
    nx = floor(Nx/2);   % x方向の中央grid
    % 各チャネルの y座標を後で計算して ny_1, ny_2, ... として格納します
    
    %% (B) min_G_Emax = 0 or 1 の2種類を試す (2チャネルコードの名残)
    %  実際には片方だけでもOK
    for min_G_Emax = 0  % ここではデモのために0だけ実行
        
        % 相対誘電率マップを初期化
        ER  = ones(Nx, Ny);
        MuR = ones(Nx, Ny);
        
        ER_best = ER;   % ベストな構造を記憶するために
        G_best_local = 0;
        
        % TFSF領域 b(x,y) の定義
        b = zeros(Nx,Ny);
        b(:, pos_src : pos_src + spc_pts + Lpts + gap_pts*(n_channels) + gap_gap_pts*(n_channels-1) + Lpts + spc_pts) = 1;
        
        kinc = [0, 1];    % y 方向から入射 (2チャネルコードと同じ)
        RES  = [dlx, dly];
        BC   = [-1, -1];
        NPML = [0,0, npml, npml];
        Pol  = 'Hz';
        
        % delta_device の定義 (最適化可能領域)
        %  2チャネル時は「3ブロック」の最適化領域を定義していました．
        %  一般化では (n_channels + 1) ブロックとし，
        %  それぞれのブロックの開始/終了 y を順次足していく形で定義します．
        delta_device = zeros(Nx,Ny);
        
        y_start = pos_src + spc_pts;
        % 各チャネル + 1 の数だけブロックを作成
        for i_block = 1 : (n_channels + 1)
            % ブロック i_block の y 範囲
            y_block_start = y_start;
            y_block_end   = y_start + Lpts - 1;  % -1 は「含む」形にするため
            
            delta_device(:, y_block_start:y_block_end) = 1;
            
            % ブロック後に gap を挟む
            y_start = y_block_end + 1;  % ブロック終わりの次
            if i_block <= n_channels
                % gap_pts を足す
                y_start = y_start + gap_pts;
                % さらに，チャネルがまだ残っていれば gap_gap_pts を足す
                if i_block < n_channels
                    y_start = y_start + gap_gap_pts;
                end
            end
        end
        delta_device_vec = delta_device(:);  % ベクトル化
        
        %% (C) 各チャネルの eta_k を定義
        %  2チャネルの時は eta1, eta2 を作っていましたが，
        %  nチャネルに拡張し，etaList{k} のように配列化します．
        
        etaList = cell(n_channels, 1);
        % チャネルの中心 y 座標をあらかじめ計算しておく
        ny_list = zeros(1, n_channels);
        
        y_temp = pos_src + spc_pts;  % 最初のブロックが始まる手前
        for k_ = 1 : n_channels
            % ブロック1 (Lpts) 後の gap 範囲中央にチャネルを置く
            y_temp = y_temp + Lpts;           % 最初のブロックを越える
            ny_  = floor(y_temp + gap_pts/2); % gap の中央
            ny_list(k_) = ny_;
            
            % 次のチャネルに進むために
            % gap を足し，もしまだチャネル残りがあるなら gap_gap も足す
            y_temp = y_temp + gap_pts;
            if k_ < n_channels
                y_temp = y_temp + gap_gap_pts;
            end
        end
        
        % 実際に eta_k(x,y) を作成
        for k_ = 1 : n_channels
            eta_k = zeros(Nx, Ny);
            
            % y = ny_list(k_) に対して，式:  exp(2*pi*1i * dlx*(0:Nx-1)/lambda0/beta)
            % を与える (2チャネルコードの eta1, eta2 相当)
            this_ny = ny_list(k_);
            eta_k(:, this_ny) = 1/Nx * exp(2*pi*1i * dlx*(0:Nx-1)/lambda0/beta);
            
            etaList{k_} = eta_k;  % セル配列に格納
        end
        
        display(ny_list(1))
        display(ny_list(2))
        
        
        % シミュレーション (全て真空) で基準となる E0 を計算
        [fields_empty, ~] = FDFD_TFSF(ones(Nx,Ny), MuR, RES, NPML, BC, lambda0, Pol, b, kinc);
        Ex0 = fields_empty.Ex;
        Ey0 = fields_empty.Ey;
        
        % 真空基準の E0 を，(nx, ny_list(1)) などから取得
        % (1チャネルコード,2チャネルコードでのやり方と同様)
        % ここではチャネル1の位置を参考にします (複数チャネルでもOK)
        ny_ref = ny_list(1);
        E0 = sqrt(abs(Ex0(nx, ny_ref))^2 + abs(Ey0(nx, ny_ref))^2);
        
        % 最適化に向け，初期の ER を設定 (starting に応じて)
        for ix = 1 : Nx
            for iy = 1 : Ny
                if delta_device(ix, iy) == 1
                    if starting == 1
                        ER(ix, iy) = rand*(eps - 1) + 1;
                    elseif starting == 2
                        ER(ix, iy) = eps/2 + 0.5;
                    else
                        % 0 -> vacuumのまま(=1)
                    end
                end
            end
        end
        
        % iteration記録用
        Gs     = zeros(N,1);
        E_maxs = zeros(N,1);
        phis   = zeros(N,1);
        AVM_prev = zeros(Nx, Ny);
        
        if display_plots
            figure('Name', sprintf('n-channels = %d, gap = %d nm', n_channels, gap_nm));
        end
        
        disp(['Start optimization for n_channels=', num2str(n_channels), ...
            ', min_G_Emax=', num2str(min_G_Emax)]);
        upd = textprogressbar(N);  % プログレスバー
        
        %% (D) メインの反復ループ
        for iter = 1 : N
            upd(iter);
            
            % 構造ありでFDFD_TFSF
            [fields, extra] = FDFD_TFSF(ER, MuR, RES, NPML, BC, lambda0, Pol, b, kinc);
            Ex = fields.Ex / E0;
            Ey = fields.Ey / E0;
            
            % チャネルごとの g_k を計算し，合計 g を得る
            g_sum = 0;
            for k_ = 1 : n_channels
                g_k = sum(sum( etaList{k_} .* Ex ));
                g_sum = g_sum + g_k;
            end
            G = real(g_sum);
            
            % フェーズ (参考)
            phis(iter) = angle(g_sum);
            
            % 数値微分演算子
            DEY = extra.derivatives.DEY;
            DEX = extra.derivatives.DEX;
            
            % ER_vec
            ER_vec = ER(:);
            chi = delta_device .* (ER - 1);
            
            Ox = -1i * lambda0 /(2*pi*c0) * spdiags(1./ER_vec,0,Nx*Ny,Nx*Ny) * DEY;
            Oy =  1i * lambda0 /(2*pi*c0) * spdiags(1./ER_vec,0,Nx*Ny,Nx*Ny) * DEX;
            
            % E_max 用の E_abs
            if in_material
                E_abs = (chi/(eps-1)) .* sqrt( abs(Ex).^2 + abs(Ey).^2 );
            else
                E_abs = delta_device .* sqrt( abs(Ex).^2 + abs(Ey).^2 );
            end
            E_max_val = max(E_abs(:));
            E_maxs(iter) = E_max_val;
            
            % min_G_Emax による b_aj 分岐 (2チャネルコード参照)
            % ここでは簡単化して「Gのみ最大化」を例示 (min_G_Emax = 0)
            % 実際は G/E_max や G/Sa を組み込むなどの拡張が可能
            b_aj_x = zeros(Nx*Ny,1);
            b_aj_y = zeros(Nx*Ny,1);
            
            % 例: b_aj = - sum_k eta_k_aj (min_G_Emax=0の場合)
            for k_ = 1 : n_channels
                eta_k_vec = etaList{k_}(:);
                eta_k_aj  = [eta_k_vec; zeros(Nx*Ny,1)];
                b_aj_xk   = eta_k_aj(1:Nx*Ny);
                b_aj_yk   = eta_k_aj(Nx*Ny+1:end);
                
                b_aj_x = b_aj_x - b_aj_xk;
                b_aj_y = b_aj_y - b_aj_yk;
            end
            
            % adjointソース b_aj を Ox, Oy で作用させて2次元にreshape
            b_aj_2d = reshape(Ox*b_aj_x + Oy*b_aj_y, [Nx, Ny]);
            b_aj_2d(isnan(b_aj_2d)) = 0;
            
            % adjointシミュレーション
            AF = extra.AF;  % システム行列の factor
            [fields_aj, ~] = FDFD_fast(ER, MuR, RES, NPML, BC, lambda0, Pol, b_aj_2d, AF);
            
            x_aj = fields_aj.x / E0;
            Ex_aj = reshape(x_aj(1:Nx*Ny), [Nx, Ny]);
            Ey_aj = reshape(x_aj(Nx*Ny+1:end), [Nx, Ny]);
            
            % 感度分布 AVM
            AVM = -real( (Ex.*Ex_aj + Ey.*Ey_aj) .* delta_device );
            
            % 誘電率の更新 + モメンタム項
            ER = ER + alpha*AVM + alpha*gamma*AVM_prev;
            AVM_prev = AVM;
            
            % 上限下限クリップ
            ER(ER < 1)  = 1;
            ER(ER > eps)= eps;
            
            % ベスト更新
            if abs(g_sum) > G_best_local  % absを取るかrealを取るかは好みに応じる
                G_best_local = abs(g_sum);
                ER_best = ER;
            end
            
            % 適宜可視化
            if display_plots && mod(iter, skip) == 0
                clf;
                subplot(2,2,1);
                imagesc(repmat( real(ER), 5, 1 ), [1, eps]);  % 縦方向に5回複写して見やすく
                colormap(flipud(gray));
                colorbar(); axis image;
                title('relative permittivity');
                
                subplot(2,2,2);
                plot(Gs(1:iter),'k');
                hold on; plot(iter, G, 'ro');
                xlabel('iteration'); ylabel('G');
                title('Acceleration Gradient');
                grid on;
                
                subplot(2,2,3);
                plot(E_maxs(1:iter),'b');
                hold on; plot(iter, E_max_val, 'ro');
                xlabel('iteration'); ylabel('E_{max}');
                title('Max Electric Field');
                grid on;
                
                subplot(2,2,4);
                plot(phis(1:iter),'k'); hold on;
                plot(iter, phis(iter),'ro');
                xlabel('iteration'); ylabel('\phi');
                title('acceleration phase');
                grid on;
                
                drawnow;
            end
            Gs(iter) = G;
        end
        
        %% (E) 最終構造のバイナリ化や後処理
        eps_avg = (eps+1)/2;
        ER_best(ER_best < eps_avg) = 1;
        ER_best(ER_best >= eps_avg) = eps;
        
        [fields_best, ~] = FDFD_TFSF(ER_best, MuR, RES, NPML, BC, lambda0, Pol, b, kinc);
        Ex_best = fields_best.Ex / E0;
        Ey_best = fields_best.Ey / E0;
        
        % チャネル全体での加速勾配 g_best を計算
        g_sum_best = 0;
        for k_ = 1 : n_channels
            g_k_best = sum(sum( etaList{k_} .* Ex_best ));
            g_sum_best = g_sum_best + g_k_best;
        end
        G_best_local_abs = abs(g_sum_best);
        
        % E_max
        E_abs_best = delta_device .* sqrt(abs(Ex_best).^2 + abs(Ey_best).^2);
        E_max_best = max(E_abs_best(:));
        
        % 保存用
        G_best_values = [G_best_values; G_best_local_abs];
        
        % 結果をテキスト出力
        timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
        fname = sprintf('%s/final_gap_%d_nm_nch_%d_%s.txt', ...
            output_folder_name, gap_nm, n_channels, timestamp);
        fid = fopen(fname, 'w');
        fprintf(fid, 'G_best (abs) = %f\n', G_best_local_abs);
        fprintf(fid, 'E_max = %f\n', E_max_best);
        fprintf(fid, 'Nx = %d, Ny = %d\n', Nx, Ny);
        fprintf(fid, 'n_channels = %d\n', n_channels);
        fclose(fid);
        fprintf('File saved: %s\n', fname);
        
        % ベスト構造の可視化
        if display_plots
            figure('Name','Best Structure','Visible','on');
        else
            figure('Name','Best Structure','Visible','off');
        end
        imagesc(repmat(real(ER_best),5,1), [1, eps]);
        colormap(flipud(gray)); axis image; colorbar();
        title(sprintf('Best Structure (n\_channels=%d, gap=%d nm)', n_channels, gap_nm));
        figNameBest = sprintf('%s/Best_nch_%d_gap_%d_nm_%s.png', ...
            output_folder_name, n_channels, gap_nm, timestamp);
        saveas(gcf, figNameBest);
    end
end

%% (F) gap_nm_values を変えた場合などのプロット例
figure; plot(gap_nm_values, G_best_values, '-o');
xlabel('gap (nm)'); ylabel('G\_best'); grid on;
title(sprintf('n=%d channels: G\\_best vs gap', n_channels));
saveas(gcf, sprintf('%s/G_best_vs_gap_nch%d_%s.png', ...
    output_folder_name, n_channels, datestr(now,'yyyy-mm-dd_HHMMSS')));
