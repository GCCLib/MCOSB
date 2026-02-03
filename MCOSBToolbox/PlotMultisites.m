% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [] = PlotMultisites(lat1,lat2,lon1,lon2,Sites_InfoGPS,Sites_InfoGAL,Sites_InfoBDS,StationGPS,StationGAL,StationBDS)
    StationGPS = mat2cell(StationGPS, ones(length(StationGPS), 1));
    StationGAL = mat2cell(StationGAL, ones(length(StationGAL), 1));
    StationBDS = mat2cell(StationBDS, ones(length(StationBDS), 1));
    
    tol_G=size(Sites_InfoGPS.name,2);
    inf_G=[];
    for j=1:tol_G
        if Sites_InfoGPS.coor(j,1) == 0
            continue;
        end
        if sum(ismember(StationGPS,Sites_InfoGPS.name(j)))>0
            [latitude, longitude]=XYZtoBLH(Sites_InfoGPS.coor(j,1),Sites_InfoGPS.coor(j,2),Sites_InfoGPS.coor(j,3));
            latitude=latitude/pi*180;
            longitude=longitude/pi*180;
            inf_G=[inf_G;Sites_InfoGPS.name(j), latitude, longitude];
        end
    end
    
    tol_E=size(Sites_InfoGAL.name,2);
    inf_E=[];
    for j=1:tol_E
        if Sites_InfoGAL.coor(j,1) == 0
            continue;
        end
        if sum(ismember(StationGAL,Sites_InfoGAL.name(j)))>0
            [latitude, longitude]=XYZtoBLH(Sites_InfoGAL.coor(j,1),Sites_InfoGAL.coor(j,2),Sites_InfoGAL.coor(j,3));
            latitude=latitude/pi*180;
            longitude=longitude/pi*180;
            inf_E=[inf_E;Sites_InfoGAL.name(j), latitude, longitude];
        end
    end
    
    tol_C=size(Sites_InfoBDS.name,2);
    inf_C=[];
    for j=1:tol_C
        if Sites_InfoBDS.coor(j,1) == 0
            continue;
        end
        if sum(ismember(StationBDS,Sites_InfoBDS.name(j)))>0
            [latitude, longitude]=XYZtoBLH(Sites_InfoBDS.coor(j,1),Sites_InfoBDS.coor(j,2),Sites_InfoBDS.coor(j,3));
            latitude=latitude/pi*180;
            longitude=longitude/pi*180;
            inf_C=[inf_C;Sites_InfoBDS.name(j), latitude, longitude];
        end
    end  
    m_proj('Miller Cylindrical','longitudes',[lon1 lon2], 'latitudes',[lat1,lat2]);
    m_coast('linewidth',1,'color','k');
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','none');
    m_grid('box','fancy','tickdir','out','linestyle','--','gridcolor','k','ytick',[-90:90:90],'FontSize',20);
    hold on;
    
    C_num=size(inf_C,1);
    for i=1:C_num        
        lat=cell2mat(inf_C(i,2));
        lon=cell2mat(inf_C(i,3));
        h3(1)=m_scatter(lon, lat,90,'Marker','o','MarkerEdgeColor',[1.00,0.00,0.00],'MarkerFaceColor',[1.00,0.00,0.00], 'LineWidth', 1); % plot bds sites
        hold on;        
    end
    disp('------> BDS has been plotted');
    
    E_num=size(inf_E,1);
    for i=1:E_num       
        lat=cell2mat(inf_E(i,2));
        lon=cell2mat(inf_E(i,3));
        h2(1)=m_scatter(lon, lat,60,'Marker','o','MarkerEdgeColor',[0.00,1.00,0.00],'MarkerFaceColor',[0.00,1.00,0.00], 'LineWidth', 1); % plot galileo sites
        hold on;        
    end
    disp('------> GAL has been plotted');
    
    G_num=size(inf_G,1);
    for i=1:G_num       
        lat=cell2mat(inf_G(i,2));
        lon=cell2mat(inf_G(i,3));
        h1(1)=m_scatter(lon, lat,40,'Marker','o','MarkerEdgeColor',[0.00,0.00,1.00],'MarkerFaceColor',[0.00,0.00,1.00], 'LineWidth', 1);  % plot gps sites
        hold on;
    end
    disp('------> GPS has been plotted');
    legend([h1(1),h2(1),h3(1)],'GPS','GAL','BDS','location','southwest','FontSize',20,'Box','off');
end