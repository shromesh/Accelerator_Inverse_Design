%% スクリプトの概要
%  既存コード "final_acceleration_gradients_gap_XXX_gapgap_YYY_タイムスタンプ.txt" で出力された
%  テキストファイルを読み込み，
%    1) (abs(g1)+abs(g2)) vs gap for each gap_gap
%    2) (abs(g1)+abs(g2))*gap vs gap for each gap_gap
%  の2種類のグラフをまとめて作成します．

clear; close all; clc;

%% ========================================================================
%  1. 解析対象の gap, gap_gap のリストを定義
%     （既存の最適化実行時と同じパラメータを設定してください）
% =========================================================================
gap_nm_values     = 200:20:1000;  % gap の候補 (nm)
gap_gap_nm_values = [300, 500, 700, 900];  % gap_gap の候補 (nm)

% 出力ファイルが置いてあるフォルダ
output_folder_name = 'result/double_channel_gap_step_20_gapgap_step_200_Jan19';


%% ========================================================================
%  2. データ格納用の配列を作成 (2次元: gap vs gap_gap)
%     - (abs(g1)+abs(g2))
%     - (abs(g1)+abs(g2)) * gap
% =========================================================================
ngap     = length(gap_nm_values);
ngapgap  = length(gap_gap_nm_values);

% 下記2つを格納するための変数を用意しておく
G_abs_sums_2D          = zeros(ngap, ngapgap);  % => (abs(g1)+abs(g2))
G_abs_sums_times_gap_2D = zeros(ngap, ngapgap); % => (abs(g1)+abs(g2)) * gap


%% ========================================================================
%  3. テキストファイルを読み込み，必要な値をパース(解析)して配列に格納
% =========================================================================
for iGap = 1:ngap
    gap_nm = gap_nm_values(iGap);
    
    for jGapGap = 1:ngapgap
        gap_gap_nm = gap_gap_nm_values(jGapGap);
        
        %------------------------------------------------------------------
        % 3-1. ファイル名のパターンを作成
        %      => "final_acceleration_gradients_gap_xxx_gapgap_yyy_*.txt" を想定
        %         (ワイルドカード * でタイムスタンプ部分を吸収）
        %------------------------------------------------------------------
        file_pattern = sprintf('final_acceleration_gradients_gap_%d_gapgap_%d_*.txt', ...
            gap_nm, gap_gap_nm);
        
        % フォルダ内を検索
        file_list = dir(fullfile(output_folder_name, file_pattern));
        if isempty(file_list)
            % 見つからなければスキップ (あるいはWarningなど)
            fprintf('Warning: %s が見つかりませんでした。\n', file_pattern);
            continue;
        end
        
        % 同じ gap, gap_gap で複数ファイルある場合は，とりあえず最初の1つを読む
        target_file = fullfile(output_folder_name, file_list(1).name);
        
        %------------------------------------------------------------------
        % 3-2. ファイルを読み込み，テキスト解析
        %      既存コードの出力例:
        %
        %       G_best (abs): XXXXXX
        %       G1_best (abs): XXXXXX
        %       G2_best (abs): XXXXXX
        %       ...
        %       (abs(g1)+abs(g2))*gap: XXXXXX
        %       abs(g1+g2)*gap: XXXXXX
        %       ...
        %
        %      ここでは (abs(g1)+abs(g2)) と (abs(g1)+abs(g2))*gap が欲しい
        %      既存コードの出力行:
        %         "G1_best (abs): %f"
        %         "G2_best (abs): %f"
        %         "(abs(g1)+abs(g2))*gap: %f"
        %
        %      を正規表現などで取り出す例を示す．
        %------------------------------------------------------------------
        file_text = fileread(target_file);
        
        % 正規表現で "G1_best (abs): 数値" の部分を抜き出す
        tokens_G1 = regexp(file_text, 'G1_best \(abs\): ([0-9e\+\-\.]+)', 'tokens', 'once');
        tokens_G2 = regexp(file_text, 'G2_best \(abs\): ([0-9e\+\-\.]+)', 'tokens', 'once');
        tokens_sum_gap = regexp(file_text, '\(abs\(g1\)\+abs\(g2\)\)\*gap: ([0-9e\+\-\.]+)', 'tokens', 'once');
        
        if isempty(tokens_G1) || isempty(tokens_G2) || isempty(tokens_sum_gap)
            % ファイルの書式が想定と違う場合はスキップやエラーなど
            fprintf('Warning: %s から必要な情報を取得できませんでした。\n', target_file);
            continue;
        end
        
        % 数値に変換
        val_G1_abs = str2double(tokens_G1{1});
        val_G2_abs = str2double(tokens_G2{1});
        val_sum_times_gap = str2double(tokens_sum_gap{1});
        
        % (abs(g1)+abs(g2)) を計算
        val_sum = val_G1_abs + val_G2_abs;
        
        %------------------------------------------------------------------
        % 3-3. 結果を 2次元配列に格納
        %------------------------------------------------------------------
        G_abs_sums_2D(iGap, jGapGap)          = val_sum;
        G_abs_sums_times_gap_2D(iGap, jGapGap) = val_sum_times_gap;  % 既に gap 済みの値
    end
end


%% ========================================================================
%  4. プロット: (abs(g1)+abs(g2)) vs gap for each gap_gap
% ========================================================================
figure('Name','(abs(g1)+abs(g2)) vs gap for each gap\_gap');
hold on; grid on;

for jGapGap = 1:ngapgap
    % gap_nm_values を x 軸にして，
    % G_abs_sums_2D(:, jGapGap) を y 値として折れ線プロット
    plot(gap_nm_values, G_abs_sums_2D(:, jGapGap), '-o', ...
        'DisplayName', sprintf('gap\\_gap = %d nm', gap_gap_nm_values(jGapGap)));
end

xlabel('gap (nm)');
ylabel('(abs(g1)+abs(g2))');
title('(abs(g1)+abs(g2)) vs gap for each gap\_gap');
legend('show');


%% ========================================================================
%  5. プロット: (abs(g1)+abs(g2))*gap vs gap for each gap_gap
%    ※ こちらが「新たに出力したかった」グラフ
% ========================================================================
figure('Name','(abs(g1)+abs(g2))*gap vs gap for each gap\_gap');
hold on; grid on;

for jGapGap = 1:ngapgap
    % gap_nm_values を x 軸にして，
    % G_abs_sums_times_gap_2D(:, jGapGap) を y 値として折れ線プロット
    plot(gap_nm_values, G_abs_sums_times_gap_2D(:, jGapGap), '-o', ...
        'DisplayName', sprintf('gap\\_gap = %d nm', gap_gap_nm_values(jGapGap)));
end

xlabel('gap (nm)');
ylabel('(abs(g1)+abs(g2)) \times gap');
title('(abs(g1)+abs(g2)) \times gap vs gap for each gap\_gap');
legend('show');

