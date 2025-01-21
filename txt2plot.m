%% スクリプトの概要
%  既存の "final_acceleration_gradients_gap_XXX_gapgap_YYY_..." テキストファイルを読み込み、
%    1) (abs(g1)+abs(g2)) vs gap for each gap_gap
%    2) (abs(g1)+abs(g2))*gap vs gap for each gap_gap
%  の2つのグラフを作成して、画面に表示したうえで PNG 保存も行います。

clear; close all; clc;

%% 1. gap, gap_gap のリスト
gap_nm_values     = 200:20:1000;  % gap の候補 (nm)
gap_gap_nm_values = [300, 500, 700, 900];  % gap_gap の候補 (nm)

% テキストファイルが置かれているフォルダ
output_folder_name = 'result/double_channel_gap_step_20_gapgap_step_200_Jan19';

%% 2. データ格納用の配列 (gap vs gap_gap)
ngap    = length(gap_nm_values);
ngapgap = length(gap_gap_nm_values);

G_abs_sums_2D           = zeros(ngap, ngapgap);  % (abs(g1)+abs(g2))
G_abs_sums_times_gap_2D = zeros(ngap, ngapgap);  % (abs(g1)+abs(g2))*gap

%% 3. テキストファイルを読み込んで、(abs(g1)+abs(g2)) と (abs(g1)+abs(g2))*gap を取得
for iGap = 1:ngap
    gap_nm = gap_nm_values(iGap);
    
    for jGapGap = 1:ngapgap
        gap_gap_nm = gap_gap_nm_values(jGapGap);
        
        %---------------------------------------------------------
        % 3-1. ファイル名パターン
        %      final_acceleration_gradients_gap_xxx_gapgap_yyy_*.txt
        %---------------------------------------------------------
        file_pattern = sprintf('final_acceleration_gradients_gap_%d_gapgap_%d_*.txt', ...
            gap_nm, gap_gap_nm);
        file_list = dir(fullfile(output_folder_name, file_pattern));
        
        if isempty(file_list)
            fprintf('Warning: %s が見つかりません。\n', file_pattern);
            continue;
        end
        
        % 同じ条件で複数ファイルあれば最初の1つを使う
        target_file = fullfile(output_folder_name, file_list(1).name);
        
        %---------------------------------------------------------
        % 3-2. テキスト読み込み・正規表現
        %---------------------------------------------------------
        file_text = fileread(target_file);
        
        % G1_best (abs): xxx
        tokens_G1 = regexp(file_text, 'G1_best \(abs\): ([0-9e\+\-\.]+)', 'tokens', 'once');
        % G2_best (abs): xxx
        tokens_G2 = regexp(file_text, 'G2_best \(abs\): ([0-9e\+\-\.]+)', 'tokens', 'once');
        % (abs(g1)+abs(g2))*gap: xxx
        tokens_sum_gap = regexp(file_text, '\(abs\(g1\)\+abs\(g2\)\)\*gap: ([0-9e\+\-\.]+)', 'tokens', 'once');
        
        if isempty(tokens_G1) || isempty(tokens_G2) || isempty(tokens_sum_gap)
            fprintf('Warning: %s から必要な情報を取得できません。\n', target_file);
            continue;
        end
        
        % 数値に変換
        val_G1_abs = str2double(tokens_G1{1});
        val_G2_abs = str2double(tokens_G2{1});
        val_sum_times_gap = str2double(tokens_sum_gap{1});
        
        % (abs(g1)+abs(g2))
        val_sum = val_G1_abs + val_G2_abs;
        
        %---------------------------------------------------------
        % 3-3. 2次元配列に格納
        %---------------------------------------------------------
        G_abs_sums_2D(iGap, jGapGap)           = val_sum;
        G_abs_sums_times_gap_2D(iGap, jGapGap) = val_sum_times_gap;  % 既に gap をかけた値
    end
end

%% 4. グラフ1: (abs(g1)+abs(g2)) vs gap for each gap_gap
figure('Name','(abs(g1)+abs(g2)) vs gap for each gap\_gap');
hold on; grid on;

for jGapGap = 1:ngapgap
    plot(gap_nm_values, G_abs_sums_2D(:, jGapGap), '-o', ...
        'DisplayName', sprintf('gap\\_gap = %d nm', gap_gap_nm_values(jGapGap)));
end

xlabel('gap (nm)');
ylabel('(abs(g1)+abs(g2))');
title('(abs(g1)+abs(g2)) vs gap for each gap\_gap');
legend('show');

% 保存 (例: PNG 形式)
timestamp_str = datestr(now,'yyyy-mm-dd_HHMMSS');  % 一意の名前にしたい場合
save_filename1 = fullfile(output_folder_name, ...
    sprintf('plot_abs_g1_plus_g2_vs_gap_for_each_gapgap_%s.png', timestamp_str));
saveas(gcf, save_filename1);
fprintf('Saved figure: %s\n', save_filename1);

%% 5. グラフ2: (abs(g1)+abs(g2))*gap vs gap for each gap_gap
figure('Name','(abs(g1)+abs(g2))*gap vs gap for each gap\_gap');
hold on; grid on;

for jGapGap = 1:ngapgap
    plot(gap_nm_values, G_abs_sums_times_gap_2D(:, jGapGap), '-o', ...
        'DisplayName', sprintf('gap\\_gap = %d nm', gap_gap_nm_values(jGapGap)));
end

xlabel('gap (nm)');
ylabel('(abs(g1)+abs(g2)) * gap');
title('(abs(g1)+abs(g2)) * gap vs gap for each gap\_gap');
legend('show');

% 保存 (例: PNG 形式)
save_filename2 = fullfile(output_folder_name, ...
    sprintf('plot_abs_g1_plus_g2_times_gap_vs_gap_for_each_gapgap_%s.png', timestamp_str));
saveas(gcf, save_filename2);
fprintf('Saved figure: %s\n', save_filename2);
