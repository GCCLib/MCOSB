% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function PlotOSBSatSTD(Data, PRN)
    % Threshold for detecting jumps
    Thres = 2.5;
    % osb(:,x,:) x is code type
    [satGC1C, osbGC1C] = non_zero_rows(squeeze(Data(:,  1, :)), PRN, 'GPS');
    [satGC1W, osbGC1W] = non_zero_rows(squeeze(Data(:,  3, :)), PRN, 'GPS');
    [satGC2C, osbGC2C] = non_zero_rows(squeeze(Data(:, 14, :)), PRN, 'GPS');
    [satGC2W, osbGC2W] = non_zero_rows(squeeze(Data(:, 20, :)), PRN, 'GPS');
    [satGC2L, osbGC2L] = non_zero_rows(squeeze(Data(:, 17, :)), PRN, 'GPS');
    [satGC2S, osbGC2S] = non_zero_rows(squeeze(Data(:, 16, :)), PRN, 'GPS');
    [satGC2X, osbGC2X] = non_zero_rows(squeeze(Data(:, 18, :)), PRN, 'GPS');
    [satGC5Q, osbGC5Q] = non_zero_rows(squeeze(Data(:, 25, :)), PRN, 'GPS');
    [satGC5X, osbGC5X] = non_zero_rows(squeeze(Data(:, 26, :)), PRN, 'GPS');

    [satEC1C, osbEC1C] = non_zero_rows(squeeze(Data(:,  1, :)), PRN, 'GAL');
    [satEC1X, osbEC1X] = non_zero_rows(squeeze(Data(:, 12, :)), PRN, 'GAL');
    [satEC6C, osbEC6C] = non_zero_rows(squeeze(Data(:, 32, :)), PRN, 'GAL');
    [satEC6X, osbEC6X] = non_zero_rows(squeeze(Data(:, 33, :)), PRN, 'GAL');
    [satEC5Q, osbEC5Q] = non_zero_rows(squeeze(Data(:, 25, :)), PRN, 'GAL');
    [satEC5X, osbEC5X] = non_zero_rows(squeeze(Data(:, 26, :)), PRN, 'GAL');
    [satEC7Q, osbEC7Q] = non_zero_rows(squeeze(Data(:, 28, :)), PRN, 'GAL');
    [satEC7X, osbEC7X] = non_zero_rows(squeeze(Data(:, 29, :)), PRN, 'GAL');
    [satEC8Q, osbEC8Q] = non_zero_rows(squeeze(Data(:, 38, :)), PRN, 'GAL');
    [satEC8X, osbEC8X] = non_zero_rows(squeeze(Data(:, 39, :)), PRN, 'GAL');

    [satCC2I, osbCC2I] = non_zero_rows(squeeze(Data(:, 40, :)), PRN, 'BDS');
    [satCC7I, osbCC7I] = non_zero_rows(squeeze(Data(:, 27, :)), PRN, 'BDS');
    [satCC6I, osbCC6I] = non_zero_rows(squeeze(Data(:, 42, :)), PRN, 'BDS');
    [satCC1X, osbCC1X] = non_zero_rows(squeeze(Data(:, 12, :)), PRN, 'BDS');
    [satCC1P, osbCC1P] = non_zero_rows(squeeze(Data(:,  2, :)), PRN, 'BDS');
    [satCC5X, osbCC5X] = non_zero_rows(squeeze(Data(:, 26, :)), PRN, 'BDS');
    [satCC5P, osbCC5P] = non_zero_rows(squeeze(Data(:, 58, :)), PRN, 'BDS');
    [satCC7Z, osbCC7Z] = non_zero_rows(squeeze(Data(:, 63, :)), PRN, 'BDS');
    [satCC8X, osbCC8X] = non_zero_rows(squeeze(Data(:, 39, :)), PRN, 'BDS');
      
    % GPS
    Gcode = {'C1C', 'C1W', 'C2C', 'C2W', 'C2L', 'C2S', 'C2X', 'C5Q', 'C5X'};
    Gosb  = {osbGC1C, osbGC1W, osbGC2C, osbGC2W, osbGC2L, osbGC2S, osbGC2X, osbGC5Q, osbGC5X};
    Gsat  = {satGC1C, satGC1W, satGC2C, satGC2W, satGC2L, satGC2S, satGC2X, satGC5Q, satGC5X}; 
    Gsats = [];
    for i = 1:length(Gsat)
        Gsats = [Gsats; Gsat{i}];
    end
    Gsats = unique(Gsats); Gsats = sort(Gsats);
    Gstd = NaN(length(Gcode), length(Gsats));
    % Calculate STD for each code type
    for idx = 1:length(Gcode)
        osb_data = Gosb{idx};
        sat_data = Gsat{idx}; 
        % Calculate STD for each satellite
        for sat = 1:length(sat_data)
            col = Gsats == sat_data(sat); 
            data = osb_data(sat, :);
            if ~isempty(data)
                % delete zeros
                data = data(:, ~all(data==0, 1));
                data = filloutliers(data, NaN, 'mean'); data = rmmissing(data);
                data = filloutliers(data, NaN, 'mean'); data = rmmissing(data);
                data = filloutliers(data, NaN, 'mean'); data = rmmissing(data);
                % jump detection with continuous condition
                diff_data = [0 diff(data)];
                % find all points where diff_data exceeds threshold
                exceed_idx = find(abs(diff_data) > Thres);
                if isempty(exceed_idx)
                    std_val = std(data);
                    Gstd(idx, col) = std_val;
                    continue;
                end
                [~,jj] = size(data); count = 1; jcount = 1; i = 1; jump_idx = [];
                while i <= (length(exceed_idx))
                    if jj < exceed_idx(i) + 20
                        data(exceed_idx(i) : end) = [];
                        break;
                    end
                    if abs(data(exceed_idx(i)+jcount)-data(exceed_idx(i)-3)) > Thres
                        count = count + 1;
                    end
                    jcount = jcount + 1;
                    if jcount > 10
                        i = i + 1; jcount = 1; count = 1;
                        continue;
                    end
                    if count > 7
                        jump_idx = [jump_idx, exceed_idx(i)];
                        i = i+1; jcount = 1; count = 1; continue;
                    end
                end                
                if isempty(jump_idx)
                    std_val = std(data);
                else            
                    i = 0;
                    while i <= (length(jump_idx)-1)
                        i = i + 1;
                        if length(jump_idx) > 1 && i < length(jump_idx) && (jump_idx(i+1)- jump_idx(i)) < 20
                            data(jump_idx(i):jump_idx(i+1)) = 0;
                            continue;
                        end
                    end
                    % apply bias correction for each jump
                    for j = 1:length(jump_idx)
                        bias = diff_data(jump_idx(j));
                        data(jump_idx(j):end) = data(jump_idx(j):end) - bias;
                    end
                    data(:, all(data == 0, 1)) = [];
                    std_val = std(data);
                end
                Gstd(idx, col) = std_val;
            end
        end
    end  
    Gstd(isnan(Gstd)) = 0;
    figure('Position', [100, 100, 1000, 600]);
    cmap = jet(256);
    cmap_with_white = [1 1 1; cmap];
    colormap(cmap_with_white);
    Gdisplay = Gstd;
    min_val = min(Gstd(:), [], 'omitnan');
    max_val = max(Gstd(:), [], 'omitnan');
    if isnan(min_val) || isnan(max_val)
        min_val = 0;
        max_val = 1;
    end
    Gdisplay(isnan(Gdisplay)) = min_val - 0.1 * (max_val - min_val);
    imagesc(Gdisplay);
    colorbar;
    set(gca, 'YTick', 1:length(Gcode), 'YTickLabel', Gcode);
    set(gca, 'XTick', 1:length(Gsats), 'XTickLabel', arrayfun(@num2str, Gsats, 'UniformOutput', false));
    xlabel('PRN', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('OSB type', 'FontSize', 12, 'FontWeight', 'bold');
    title('GPS STD', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca,'FontSize',12);
    % add std value
    for i = 1:size(Gstd, 1)
        for j = 1:size(Gstd, 2)
            if (Gstd(i, j)) ~= 0
                text(j, i, sprintf('%.1f', Gstd(i, j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', 'k');
            end
        end
    end
    clim([0, 1]);

    % Galileo
    Ecode = {'C1C', 'C1X', 'C6C', 'C6X', 'C5Q', 'C5X', 'C7Q', 'C7X', 'C8Q', 'C8X'};
    Eosb = {osbEC1C, osbEC1X, osbEC6C, osbEC6X, osbEC5Q, osbEC5X, osbEC7Q, osbEC7X, osbEC8Q, osbEC8X};
    Esat = {satEC1C, satEC1X, satEC6C, satEC6X, satEC5Q, satEC5X, satEC7Q, satEC7X, satEC8Q, satEC8X};
    Esats = [];
    for i = 1:length(Esat)
        Esats = [Esats; Esat{i}];
    end
    Esats = unique(Esats); Esats = sort(Esats);
    Estd = NaN(length(Ecode), length(Esats));
    % Calculate STD for each code type
    for idx = 1:length(Ecode)
        osb_data = Eosb{idx};
        sat_data = Esat{idx}; 
        % Calculate STD for each satellite
        for sat = 1:length(sat_data)
            col = Esats == sat_data(sat); 
            data = osb_data(sat, :);
            if ~isempty(data)
                % delete zeros
                data = data(:, ~all(data==0, 1));
                data = filloutliers(data, NaN, 'mean'); data = rmmissing(data);
                data = filloutliers(data, NaN, 'mean'); data = rmmissing(data);
                data = filloutliers(data, NaN, 'mean'); data = rmmissing(data);
                % jump detection with continuous condition
                diff_data = [0 diff(data)];
                % find all points where diff_data exceeds threshold
                exceed_idx = find(abs(diff_data) > Thres);
                if isempty(exceed_idx)
                    std_val = std(data);
                    Estd(idx, col) = std_val;
                    continue;
                end
                [~,jj] = size(data); count = 1; jcount = 1; i = 1; jump_idx = [];
                while i <= (length(exceed_idx))
                    if jj < exceed_idx(i) + 20
                        data(exceed_idx(i) : end) = [];
                        break;
                    end
                    if abs(data(exceed_idx(i)+jcount)-data(exceed_idx(i)-3)) > Thres
                        count = count + 1;
                    end
                    jcount = jcount + 1;
                    if jcount > 10
                        i = i + 1; jcount = 1; count = 1;
                        continue;
                    end
                    if count > 7
                        jump_idx = [jump_idx, exceed_idx(i)];
                        i = i+1; jcount = 1; count = 1; continue;
                    end
                end                
                if isempty(jump_idx)
                    std_val = std(data);
                else            
                    i = 0;
                    while i <= (length(jump_idx)-1)
                        i = i + 1;
                        if length(jump_idx) > 1 && i < length(jump_idx) && (jump_idx(i+1)- jump_idx(i)) < 20
                            data(jump_idx(i):jump_idx(i+1)) = 0;
                            continue;
                        end
                    end
                    % apply bias correction for each jump
                    for j = 1:length(jump_idx)
                        bias = diff_data(jump_idx(j));
                        data(jump_idx(j):end) = data(jump_idx(j):end) - bias;
                    end
                    data(:, all(data == 0, 1)) = [];
                    std_val = std(data);
                end
                Estd(idx, col) = std_val;
            end
        end
    end  
    Estd(isnan(Estd)) = 0;
    figure('Position', [100, 100, 1000, 600]);
    cmap = jet(256);
    cmap_with_white = [1 1 1; cmap];
    colormap(cmap_with_white);
    Edisplay = Estd;
    min_val = min(Estd(:), [], 'omitnan');
    max_val = max(Estd(:), [], 'omitnan');
    if isnan(min_val) || isnan(max_val)
        min_val = 0;
        max_val = 1;
    end
    Edisplay(isnan(Edisplay)) = min_val - 0.1 * (max_val - min_val);
    imagesc(Edisplay);
    colorbar;
    set(gca, 'YTick', 1:length(Ecode), 'YTickLabel', Ecode);
    set(gca, 'XTick', 1:length(Esats), 'XTickLabel', arrayfun(@num2str, Esats, 'UniformOutput', false));
    xlabel('PRN', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('OSB type', 'FontSize', 12, 'FontWeight', 'bold');
    title('GAL STD', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca,'FontSize',12);
    % add std value
    for i = 1:size(Estd, 1)
        for j = 1:size(Estd, 2)
            if (Estd(i, j)) ~= 0
                text(j, i, sprintf('%.1f', Estd(i, j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', 'k');
            end
        end
    end
    clim([0, 1]);

    % BDS 
    Ccode = {'C2I', 'C7I', 'C6I', 'C1X', 'C1P', 'C5X', 'C5P', 'C7Z', 'C8X'};
    Cosb = {osbCC2I, osbCC7I, osbCC6I, osbCC1X, osbCC1P, osbCC5X, osbCC5P, osbCC7Z, osbCC8X};
    Csat = {satCC2I, satCC7I, satCC6I, satCC1X, satCC1P, satCC5X, satCC5P, satCC7Z, satCC8X};
    Csats = [];
    for i = 1:length(Csat)
        Csats = [Csats; Csat{i}];
    end
    Csats = unique(Csats); Csats = sort(Csats);
    Cstd = NaN(length(Ccode), length(Csats));
    % Calculate STD for each code type
    for idx = 1:length(Ccode)
        osb_data = Cosb{idx};
        sat_data = Csat{idx}; 
        % Calculate STD for each satellite
        for sat = 1:length(sat_data)
            col = Csats == sat_data(sat); 
            data = osb_data(sat, :);
            if ~isempty(data)
                % delete zeros
                data = data(:, ~all(data==0, 1));
                data = filloutliers(data, NaN, 'mean'); data = rmmissing(data);
                data = filloutliers(data, NaN, 'mean'); data = rmmissing(data);
                data = filloutliers(data, NaN, 'mean'); data = rmmissing(data);
                % jump detection with continuous condition
                diff_data = [0 diff(data)];
                % find all points where diff_data exceeds threshold
                exceed_idx = find(abs(diff_data) > Thres);
                if isempty(exceed_idx)
                    std_val = std(data);
                    Cstd(idx, col) = std_val;
                    continue;
                end
                [~,jj] = size(data); count = 1; jcount = 1; i = 1; jump_idx = [];
                while i <= (length(exceed_idx))
                    if jj < exceed_idx(i) + 20
                        data(exceed_idx(i) : end) = [];
                        break;
                    end
                    if abs(data(exceed_idx(i)+jcount)-data(exceed_idx(i)-3)) > Thres
                        count = count + 1;
                    end
                    jcount = jcount + 1;
                    if jcount > 10
                        i = i + 1; jcount = 1; count = 1;
                        continue;
                    end
                    if count > 7
                        jump_idx = [jump_idx, exceed_idx(i)];
                        i = i+1; jcount = 1; count = 1; continue;
                    end
                end                
                if isempty(jump_idx)
                    std_val = std(data);
                else            
                    i = 0;
                    while i <= (length(jump_idx)-1)
                        i = i + 1;
                        if length(jump_idx) > 1 && i < length(jump_idx) && (jump_idx(i+1)- jump_idx(i)) < 20
                            data(jump_idx(i):jump_idx(i+1)) = 0;
                            continue;
                        end
                    end
                    % apply bias correction for each jump
                    for j = 1:length(jump_idx)
                        bias = diff_data(jump_idx(j));
                        data(jump_idx(j):end) = data(jump_idx(j):end) - bias;
                    end
                    data(:, all(data == 0, 1)) = [];
                    std_val = std(data);
                end
                Cstd(idx, col) = std_val;
            end
        end
    end  
    Cstd(isnan(Cstd)) = 0;
    figure('Position', [100, 100, 1000, 600]);
    cmap = jet(256);
    cmap_with_white = [1 1 1; cmap];
    colormap(cmap_with_white);
    Cdisplay = Cstd;
    min_val = min(Cstd(:), [], 'omitnan');
    max_val = max(Cstd(:), [], 'omitnan');
    if isnan(min_val) || isnan(max_val)
        min_val = 0;
        max_val = 1;
    end
    Cdisplay(isnan(Cdisplay)) = min_val - 0.1 * (max_val - min_val);
    imagesc(Cdisplay);
    colorbar;
    set(gca, 'YTick', 1:length(Ccode), 'YTickLabel', Ccode);
    set(gca, 'XTick', 1:length(Csats), 'XTickLabel', arrayfun(@num2str, Csats, 'UniformOutput', false));
    xlabel('PRN', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('OSB type', 'FontSize', 12, 'FontWeight', 'bold');
    title('BDS STD', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca,'FontSize',12);
    % add std value
    for i = 1:size(Cstd, 1)
        for j = 1:size(Cstd, 2)
            if (Cstd(i, j)) ~= 0
                text(j, i, sprintf('%.2f', Cstd(i, j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 8, 'Color', 'k');
            end
        end
    end
    clim([0, 1]);
end
function [sat, new_matrix] = non_zero_rows(matrix, prn, sys)
    if isempty(matrix)
        sat = [];
        new_matrix = [];
    else
        if contains(sys, 'GPS')
            matrix(33:end, :) = 0;
        end
        if contains(sys, 'GAL')
            matrix(1:32, :)   = 0;
            matrix(64:end, :) = 0;
        end
        if contains(sys, 'BDS')
            matrix(1:63, :)   = 0;
        end    
        % find non_zero_rows
        ind = any(matrix, 2);
        sat = prn(ind);
        sat = sat';
        new_matrix = matrix(ind, :);
    end
end