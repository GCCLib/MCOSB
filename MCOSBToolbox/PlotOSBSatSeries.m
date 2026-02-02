% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function PlotOSBSatSeries(Data, PRN, colors)
    % osb(:,x,:) x is code type
    [satGC1C, osbGC1C] = non_zero_rows(squeeze(Data(:,   1, :)), PRN, 'GPS');
    [satGC1W, osbGC1W] = non_zero_rows(squeeze(Data(:,   3, :)), PRN, 'GPS');
    [satGC2C, osbGC2C] = non_zero_rows(squeeze(Data(:,  14, :)), PRN, 'GPS');
    [satGC2W, osbGC2W] = non_zero_rows(squeeze(Data(:,  20, :)), PRN, 'GPS');
    [satGC2L, osbGC2L] = non_zero_rows(squeeze(Data(:,  17, :)), PRN, 'GPS');
    [satGC2S, osbGC2S] = non_zero_rows(squeeze(Data(:,  16, :)), PRN, 'GPS');
    [satGC2X, osbGC2X] = non_zero_rows(squeeze(Data(:,  18, :)), PRN, 'GPS');
    [satGC5Q, osbGC5Q] = non_zero_rows(squeeze(Data(:,  25, :)), PRN, 'GPS');
    [satGC5X, osbGC5X] = non_zero_rows(squeeze(Data(:,  26, :)), PRN, 'GPS');

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
    if ~isempty(osbGC1C)
        figure;
        temp_data = osbGC1C;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satGC1C);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GPS C1C', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end
    if ~isempty(osbGC1W)
        figure;
        temp_data = osbGC1W;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satGC1W);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GPS C1W', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end
    if ~isempty(osbGC2C)
        figure;
        temp_data = osbGC2C;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satGC2C);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GPS C2C', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end
    if ~isempty(osbGC2W)
        figure;
        temp_data = osbGC2W;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satGC2W);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GPS C2W', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end
    if ~isempty(osbGC2L)
        figure;
        temp_data = osbGC2L;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satGC2L);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GPS C2L', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end
    if ~isempty(osbGC2S)
        figure;
        temp_data = osbGC2S;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satGC2S);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GPS C2S', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end    
    if ~isempty(osbGC2X)
        figure;
        temp_data = osbGC2X;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satGC2X);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GPS C2X', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end
    if ~isempty(osbGC5Q)
        figure;
        temp_data = osbGC5Q;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satGC5Q);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GPS C5Q', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end    
    if ~isempty(osbGC5X)
        figure;
        temp_data = osbGC5X;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satGC5X);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GPS C5X', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end 
    
    % Galileo
    if ~isempty(osbEC1C)
        figure;
        temp_data = osbEC1C;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satEC1C);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GAL C1C', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end
    if ~isempty(osbEC1X)
        figure;
        temp_data = osbEC1X;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satEC1X);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GAL C1X', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end
    if ~isempty(osbEC6C)
        figure;
        temp_data = osbEC6C;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satEC6C);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GAL C6C', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end    
    if ~isempty(osbEC6X)
        figure;
        temp_data = osbEC6X;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satEC6X);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GAL C6X', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end
    if ~isempty(osbEC5Q)
        figure;
        temp_data = osbEC5Q;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satEC5Q);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GAL C5Q', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end    
    if ~isempty(osbEC5X)
        figure;
        temp_data = osbEC5X;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satEC5X);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GAL C5X', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end          
    if ~isempty(osbEC7Q)
        figure;
        temp_data = osbEC7Q;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satEC7Q);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GAL C7Q', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end  
    if ~isempty(osbEC7X)
        figure;
        temp_data = osbEC7X;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satEC7X);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GAL C7X', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end 
    if ~isempty(osbEC8Q)
        figure;
        temp_data = osbEC8Q;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satEC8Q);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GAL C8Q', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end
    if ~isempty(osbEC8X)
        figure;
        temp_data = osbEC8X;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satEC8X);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('GAL C8X', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end    
    % BDS
    if ~isempty(osbCC2I)
        figure;
        temp_data = osbCC2I;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satCC2I);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('BDS C2I', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end
    if ~isempty(osbCC7I)
        figure;
        temp_data = osbCC7I;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satCC7I);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('BDS C7I', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end    
    if ~isempty(osbCC6I)
        figure;
        temp_data = osbCC6I;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satCC6I);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('BDS C6I', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end      
    if ~isempty(osbCC1X)
        figure;
        temp_data = osbCC1X;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satCC1X);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('BDS C1X', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end     
    if ~isempty(osbCC1P)
        figure;
        temp_data = osbCC1P;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satCC1P);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('BDS C1P', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end    
    if ~isempty(osbCC5X)
        figure;
        temp_data = osbCC5X;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satCC5X);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('BDS C5X', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end
    if ~isempty(osbCC5P)
        figure;
        temp_data = osbCC5P;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satCC5P);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('BDS C5P', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end    
    if ~isempty(osbCC7Z)
        figure;
        temp_data = osbCC7Z;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satCC7Z);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('BDS C7Z', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
    end    
    if ~isempty(osbCC8X)
        figure;
        temp_data = osbCC8X;
        temp_data(temp_data == 0) = NaN;
        [~, days] = size(temp_data);
        [sats, ~] = size(satCC8X);
        count = 0;
        for i=1:sats
            count = count + 1;
            scatter(1:days,temp_data(i,:),[],colors(count,:),'filled');
            hold on;box on;grid on;
        end
        xlabel('DOY', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('OSB [ns]', 'FontSize', 12, 'FontWeight', 'bold');
        title('BDS C8X', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca,'FontSize',12);
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
