% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [StationGPS,StationGAL,StationBDS,Sites_InfoGPS,Sites_InfoGAL,Sites_InfoBDS] = GetStation(MCOSBPath,DOY)
    StationGPS = [];
    Path = [MCOSBPath,'\MCB\GPS\',num2str(DOY,'%03d'),'\GPSC5Q'];
    list=dir([Path '/*.mat']);
    for i = 1:length(list)
        StationGPS = [StationGPS;list(i).name(1:4)];
    end
    Path = [MCOSBPath,'\MCB\GPS\',num2str(DOY,'%03d'),'\GPSC5X'];
    list=dir([Path '/*.mat']);
    for i = 1:length(list)
        StationGPS = [StationGPS;list(i).name(1:4)];
    end 
    Path = [MCOSBPath,'\STA\GPS\',num2str(DOY,'%03d')];
    list=dir([Path '/*.mat']);
    FullFilePath = fullfile(Path, list(1).name);
    Sites_InfoGPS = load(FullFilePath);

    StationGAL = [];
    Path = [MCOSBPath,'\MCB\GAL\',num2str(DOY,'%03d'),'\GALC1C'];
    list=dir([Path '/*.mat']);
    for i = 1:length(list)
        StationGAL = [StationGAL;list(i).name(1:4)];
    end
    Path = [MCOSBPath,'\MCB\GAL\',num2str(DOY,'%03d'),'\GALC1X'];
    list=dir([Path '/*.mat']);
    for i = 1:length(list)
        StationGAL = [StationGAL;list(i).name(1:4)];
    end 
    Path = [MCOSBPath,'\STA\GAL\',num2str(DOY,'%03d')];
    list=dir([Path '/*.mat']);
    FullFilePath = fullfile(Path, list(1).name);
    Sites_InfoGAL = load(FullFilePath);

    StationBDS = [];
    Path = [MCOSBPath,'\MCB\BDS3\',num2str(DOY,'%03d'),'\BDSC1P'];
    list=dir([Path '/*.mat']);
    for i = 1:length(list)
        StationBDS = [StationBDS;list(i).name(1:4)];
    end
    Path = [MCOSBPath,'\MCB\BDS3\',num2str(DOY,'%03d'),'\BDSC1X'];
    list=dir([Path '/*.mat']);
    for i = 1:length(list)
        StationBDS = [StationBDS;list(i).name(1:4)];
    end    
    Path = [MCOSBPath,'\STA\BDS3\',num2str(DOY,'%03d')];
    list=dir([Path '/*.mat']);
    FullFilePath = fullfile(Path, list(1).name);
    Sites_InfoBDS = load(FullFilePath);    
end