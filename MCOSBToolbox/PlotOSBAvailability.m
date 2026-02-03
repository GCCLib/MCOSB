% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function PlotOSBAvailability(Data, PRN)
    % osb(:,x,:) x is code type
    [sat_C1W, osb_C1W] = non_zero_rows(squeeze(Data(:,  3, :)), PRN, 'GPS');
    [sat_C1C, osb_C1C] = non_zero_rows(squeeze(Data(:,  1, :)), PRN, 'GAL');
    [sat_C2I, osb_C2I] = non_zero_rows(squeeze(Data(:, 40, :)), PRN, 'BDS');
    % 0 or 1
    A = logical(osb_C1W); A = double(A); 
    B = logical(osb_C1C); B = double(B); 
    C = logical(osb_C2I); C = double(C);
    Code_availability = [A; B; C];
    SatList = [sat_C1W; sat_C1C; sat_C2I];
    SatNum = length(SatList);
    [~, days] = size(osb_C1W);
    
    figure('Position', [100, 100, 1200, 600]);
    ax1 = axes('Position', [0.1, 0.15, 0.75, 0.75]);
    imagesc(1:SatNum, 1:days, Code_availability');
    
    % red or blue
    colormap(ax1, [1 0.3 0.3; 0.3 0.8 0.3]);
    xlabel('PRN', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('DOY', 'FontSize', 16, 'FontWeight', 'bold');
    title('OSB availability', 'FontSize', 16, 'FontWeight', 'bold');

    set(ax1, 'XTick', 1:SatNum);
    tickLabels = cell(1, SatNum);
    for i = 1:SatNum
        if mod(i, 2) == 1
            tickLabels{i} = num2str(SatList(i));
        else
            tickLabels{i} = '';
        end
    end
    set(ax1, 'XTickLabel', tickLabels);
    xtickangle(ax1, 45);

    if days > 360
        yticks = 1:30:days;
        if yticks(end) ~= days
            yticks = [yticks, days];
        end
        set(ax1, 'YTick', yticks);
    else
        yticks = 1:10:days;
        if yticks(end) ~= days
            yticks = [yticks, days];
        end
        set(ax1, 'YTick', yticks);        
    end  
    % color bar
    c = colorbar('Position', [0.87, 0.15, 0.02, 0.75]);
    c.Label.String = 'Availability';
    c.Ticks = '';
    % c.TickLabels = {'missing', 'normal'};
    c.Label.FontSize = 16;
    
    grid(ax1, 'on');
    set(ax1, 'FontSize', 16, 'FontWeight', 'normal');
    set(gcf, 'Color', 'white');
    
    if ~isempty(sat_C1W) && ~isempty(sat_C1C)
        gps_end = length(sat_C1W);
        x_line = gps_end + 0.5;
        y_limits = ylim(ax1);
        line(ax1, [x_line x_line], y_limits, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
    end 
    if ~isempty(sat_C1C) && ~isempty(sat_C2I)
        gal_end = length(sat_C1W) + length(sat_C1C);
        x_line = gal_end + 0.5;
        y_limits = ylim(ax1);
        line(ax1, [x_line x_line], y_limits, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
    end
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