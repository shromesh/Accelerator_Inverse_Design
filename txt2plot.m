clear; close all; clc;

%% 1. gap, gap_gap のリスト
gap_nm_values = 40:40:1000;
gap_gap_nm_values = 100:200:1000;

% テキストファイルが置かれているフォルダ
input_folder_name = 'result/double_channel_step_40_gapgap_step_200_jan23_parallel_grids_100_L04';
% output_folder_name = 'result/double_channel_abs_g1_g2_gap_gap_jan24';
output_folder_name = 'result/double_channel_title_gap1_gap2_jan25';


%% 2. データ格納用の配列 (gap vs gap_gap)
ngap    = length(gap_nm_values);
ngapgap = length(gap_gap_nm_values);

% abs(g1), abs(g2), abs(g1)+abs(g2), (abs(g1)+abs(g2))*gap を格納
G1_abs_2D = zeros(ngap, ngapgap);  % abs(g1)
G2_abs_2D = zeros(ngap, ngapgap);  % abs(g2)

G_abs_sums_2D           = zeros(ngap, ngapgap);  % (abs(g1) + abs(g2))
G_abs_sums_times_gap_2D = zeros(ngap, ngapgap);  % (abs(g1) + abs(g2)) * gap

% g1 * gap, g2 * gap を格納
G1_times_gap_2D = zeros(ngap, ngapgap);  % g1 * gap
G2_times_gap_2D = zeros(ngap, ngapgap);  % g2 * gap

%% 3. テキストファイルを読み込んで、各種値を取得
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
        file_list = dir(fullfile(input_folder_name, file_pattern));
        
        if isempty(file_list)
            fprintf('Warning: %s が見つかりません。\n', file_pattern);
            continue;
        end
        
        % 同じ条件で複数ファイルあれば最初の1つを使う
        target_file = fullfile(input_folder_name, file_list(1).name);
        
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
        
        % g1 * gap, g2 * gap
        val_G1_times_gap = val_G1_abs * gap_nm;
        val_G2_times_gap = val_G2_abs * gap_nm;
        
        %---------------------------------------------------------
        % 3-3. 2次元配列に格納
        %---------------------------------------------------------
        G1_abs_2D(iGap, jGapGap)                = val_G1_abs;
        G2_abs_2D(iGap, jGapGap)                = val_G2_abs;
        G_abs_sums_2D(iGap, jGapGap)            = val_sum;
        G_abs_sums_times_gap_2D(iGap, jGapGap)  = val_sum_times_gap;
        G1_times_gap_2D(iGap, jGapGap)          = val_G1_times_gap;
        G2_times_gap_2D(iGap, jGapGap)          = val_G2_times_gap;
    end
end

%% タイムスタンプ文字列 (ファイル名に付加して一意化)
timestamp_str = datestr(now,'yyyy-mm-dd_HHMMSS');

% %% 4. グラフ1: abs(g1) * gap vs gap for each gap_gap
figure('Name','abs(g1) * gap vs gap for each gap\_gap');
hold on; grid on;

for jGapGap = 1:ngapgap
    plot(gap_nm_values, G1_abs_2D(:, jGapGap) .* gap_nm_values', '-o', ...
        'DisplayName', sprintf('gap\\_gap = %d nm', gap_gap_nm_values(jGapGap)));
end

xlabel('gap (nm)');
ylabel('abs(g1) * gap');
title('abs(g1) * gap vs gap for each gap\_gap');
legend('show');

% 保存
save_filename1 = fullfile(output_folder_name, ...
    sprintf('plot_abs_g1_times_gap_vs_gap_for_each_gapgap_%s.png', timestamp_str));
saveas(gcf, save_filename1);
fprintf('Saved figure: %s\n', save_filename1);

%% 5. グラフ2: abs(g2) * gap vs gap for each gap_gap
figure('Name','abs(g2) * gap vs gap for each gap\_gap');
hold on; grid on;

for jGapGap = 1:ngapgap
    plot(gap_nm_values, G2_abs_2D(:, jGapGap) .* gap_nm_values', '-o', ...
        'DisplayName', sprintf('gap\\_gap = %d nm', gap_gap_nm_values(jGapGap)));
end

xlabel('gap (nm)');
ylabel('abs(g2) * gap');
title('abs(g2) * gap vs gap for each gap\_gap');
legend('show');

% 保存
save_filename2 = fullfile(output_folder_name, ...
    sprintf('plot_abs_g2_times_gap_vs_gap_for_each_gapgap_%s.png', timestamp_str));
saveas(gcf, save_filename2);
fprintf('Saved figure: %s\n', save_filename2);

%% 6. グラフ3: abs(g1) * gap の gap_gap の中での最大値 vs gap
figure('Name','abs(g1) * gap with best gap\_gap vs gap');
hold on; grid on;

G1_abs_times_gap_max_over_gapgap = zeros(ngap,1);
best_gapgap_indices_for_G1 = zeros(ngap, 1); % Store the best gap_gap index

for iGap = 1:ngap
    [G1_abs_times_gap_max_over_gapgap(iGap), best_gapgap_indices_for_G1(iGap)] = ...
        max(G1_abs_2D(iGap,:) .* gap_nm_values(iGap));
    
end

plot(gap_nm_values, G1_abs_times_gap_max_over_gapgap, '-o', 'DisplayName', 'Max (abs(g1) * gap)');

xlabel('gap (nm)');
ylabel('abs(g1) * gap with best gap\_gap');
title('abs(g1) * gap with best gap\_gap vs gap');

% 保存
save_filename3 = fullfile(output_folder_name, ...
    sprintf('plot_max_abs_g1_times_gap_vs_gap_%s.png', timestamp_str));
saveas(gcf, save_filename3);
fprintf('Saved figure: %s\n', save_filename3);

% 追加: abs(g1) with best gap_gap vs gap
figure('Name','abs(g1) with best gap\_gap vs gap');
hold on; grid on;

G1_abs_max_over_gapgap = zeros(ngap,1);

for iGap = 1:ngap
    G1_abs_max_over_gapgap(iGap) = G1_abs_2D(iGap, best_gapgap_indices_for_G1(iGap));
end

plot(gap_nm_values, G1_abs_max_over_gapgap, '-o', 'DisplayName', 'Max (abs(g1))');

xlabel('gap (nm)');
ylabel('abs(g1) with best gap\_gap');
title('abs(g1) with best gap\_gap vs gap');

% 保存
save_filename4 = fullfile(output_folder_name, ...
    sprintf('plot_max_abs_g1_vs_gap_%s.png', timestamp_str));
saveas(gcf, save_filename4);
fprintf('Saved figure: %s\n', save_filename4);

%% 7. グラフ4: abs(g2) * gap の gap_gap の中での最大値 vs gap
figure('Name','abs(g2) * gap with best gap\_gap vs gap');
hold on; grid on;

G2_abs_times_gap_max_over_gapgap = zeros(ngap,1);
best_gapgap_indices_for_G2 = zeros(ngap,1);

for iGap = 1:ngap
    [G2_abs_times_gap_max_over_gapgap(iGap), best_gapgap_indices_for_G2(iGap)] = ...
        max(G2_abs_2D(iGap,:) .* gap_nm_values(iGap));
end


plot(gap_nm_values, G2_abs_times_gap_max_over_gapgap, '-o', 'DisplayName', 'Max (abs(g2) * gap)');


xlabel('gap (nm)');
ylabel('abs(g2) * gap with best gap\_gap');
title('abs(g2) * gap with best gap\_gap vs gap');

% 保存
save_filename5 = fullfile(output_folder_name, ...
    sprintf('plot_max_abs_g2_times_gap_vs_gap_%s.png', timestamp_str));
saveas(gcf, save_filename5);
fprintf('Saved figure: %s\n', save_filename5);

% 追加: abs(g2) with best gap_gap vs gap
figure('Name','abs(g2) with best gap\_gap vs gap');
hold on; grid on;

G2_abs_max_over_gapgap = zeros(ngap,1);

for iGap = 1:ngap
    G2_abs_max_over_gapgap(iGap) = G2_abs_2D(iGap, best_gapgap_indices_for_G2(iGap));
end

plot(gap_nm_values, G2_abs_max_over_gapgap, '-o', 'DisplayName', 'Max (abs(g2))');

xlabel('gap (nm)');
ylabel('abs(g2) with best gap\_gap');
title('abs(g2) with best gap\_gap vs gap');

% 保存
save_filename6 = fullfile(output_folder_name, ...
    sprintf('plot_max_abs_g2_vs_gap_%s.png', timestamp_str));
saveas(gcf, save_filename6);
fprintf('Saved figure: %s\n', save_filename6);

%% 8. グラフ5: abs(g1) * gap + abs(g2) * gap vs gap for each gap_gap
figure('Name','abs(g1) * gap + abs(g2) * gap vs gap for each gap\_gap');
hold on; grid on;

for jGapGap = 1:ngapgap
    plot(gap_nm_values, G_abs_sums_times_gap_2D(:, jGapGap), '-o', ...
        'DisplayName', sprintf('gap\\_gap = %d nm', gap_gap_nm_values(jGapGap)));
end

xlabel('gap (nm)');
ylabel('abs(g1) * gap + abs(g2) * gap');
title('abs(g1) * gap + abs(g2) * gap vs gap for each gap\_gap');
legend('show');

% 保存
save_filename7 = fullfile(output_folder_name, ...
    sprintf('plot_abs_g1_times_gap_plus_abs_g2_times_gap_vs_gap_for_each_gapgap_%s.png', timestamp_str));
saveas(gcf, save_filename7);
fprintf('Saved figure: %s\n', save_filename7);

%% 9. グラフ6: abs(g1) * gap + abs(g2) * gap の gap_gap の中での最大値 vs gap
figure('Name','abs(g1) * gap + abs(g2) * gap with best gap\_gap vs gap');
hold on; grid on;

G_abs_sums_times_gap_max_over_gapgap = zeros(ngap,1);
best_gapgap_indices_for_sums = zeros(ngap,1);


for iGap = 1:ngap
    [G_abs_sums_times_gap_max_over_gapgap(iGap), best_gapgap_indices_for_sums(iGap)] = ...
        max(G_abs_sums_times_gap_2D(iGap,:));
end

plot(gap_nm_values, G_abs_sums_times_gap_max_over_gapgap, '-o', 'DisplayName', 'Max ((abs(g1) + abs(g2)) * gap)');


xlabel('gap (nm)');
ylabel('abs(g1) * gap + abs(g2) * gap with best gap\_gap');
title('abs(g1) * gap + abs(g2) * gap with best gap\_gap vs gap');

% 保存
save_filename8 = fullfile(output_folder_name, ...
    sprintf('plot_max_abs_g1_times_gap_plus_abs_g2_times_gap_vs_gap_%s.png', timestamp_str));
saveas(gcf, save_filename8);
fprintf('Saved figure: %s\n', save_filename8);

% // 追加: g1 * gap vs gap for each gap_gap
figure('Name','g1 * gap vs gap for each gap\_gap');
hold on; grid on;

for jGapGap = 1:ngapgap
    plot(gap_nm_values, G1_times_gap_2D(:, jGapGap), '-o', ...
        'DisplayName', sprintf('gap\\_gap = %d nm', gap_gap_nm_values(jGapGap)));
end

xlabel('gap (nm)');
ylabel('g1 * gap');
title('g1 * gap vs gap for each gap\_gap');
legend('show');

% 保存
save_filename9 = fullfile(output_folder_name, ...
    sprintf('plot_g1_times_gap_vs_gap_for_each_gapgap_%s.png', timestamp_str));
saveas(gcf, save_filename9);
fprintf('Saved figure: %s\n', save_filename9);

% // 追加: g2 * gap vs gap for each gap_gap
figure('Name','g2 * gap vs gap for each gap\_gap');
hold on; grid on;

for jGapGap = 1:ngapgap
    plot(gap_nm_values, G2_times_gap_2D(:, jGapGap), '-o', ...
        'DisplayName', sprintf('gap\\_gap = %d nm', gap_gap_nm_values(jGapGap)));
end

xlabel('gap (nm)');
ylabel('g2 * gap');
title('g2 * gap vs gap for each gap\_gap');
legend('show');

% 保存
save_filename10 = fullfile(output_folder_name, ...
    sprintf('plot_g2_times_gap_vs_gap_for_each_gapgap_%s.png', timestamp_str));
saveas(gcf, save_filename10);
fprintf('Saved figure: %s\n', save_filename10);
