% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
clc; clear;
%------------------------------------------------------------------------------------------------------
clc;clear;
disp('Step four: analyze multi-GNSS and multi-frequency code OSB characteristics !')
Path = 'E:\MCOSB';
% e.g., CAS, ECU
Center = 'CAS';
% Plot station distribution, 1 for on and 0 for off, only for the processing based on the MCOSB software
Staion = 0;
% Plot which day
DOY = 1;
%------------------------------------------------------------------------------------------------------
SatPRN = ["G01","G02","G03","G04","G05","G06","G07","G08","G09","G10","G11",...    % 1,2,3,4,5,6,7,8,9,10,11
          "G12","G13","G14","G15","G16","G17","G18","G19","G20","G21","G22",...    % 12,13,14,15,16,17,18,19,20,21,22
          "G23","G24","G25","G26","G27","G28","G29","G30","G31","G32",...          % 23,24,25,26,27,28,29,30,31,32
          "E01","E02","E03","E04","E05","E06","E07","E08","E09","E10","E11",...    % 33,34,35,36,37,38,39,40,41,42,43
          "E12","E13","E14","E15","E16","E18","E19","E21","E22","E23","E24",...    % 44,45,46,47,48,49,50,51,52,53,54
          "E25","E26","E27","E29","E30","E31","E33","E34","E36",...                % 55,56,57,58,59,60,61,62,63
          "C01","C02","C03","C04","C05","C06","C07","C08","C09","C10","C11",...    % 64,65,66,67,68,69,70,71,72,73,74
          "C12","C13","C14","C16","C19","C20","C21","C22","C23","C24","C25",...    % 75,76,77,78,79,80,81,82,83,84,85
          "C26","C27","C28","C29","C30","C32","C33","C34","C35","C36","C37",...    % 86,87,88,89,90,91,92,93,94,95,96
          "C38","C39","C40","C41","C42","C43","C44","C45","C46"];                  % 97,98,99,100,101,102,103,104,105
%------------------------------------------------------------------------------------------------------
Signal = ["1C","1P","1W","1Y","1M","1N","1S","1L","1E","1A",...% /*  1-10 */
          "1B","1X","1Z","2C","2D","2S","2L","2X","2P","2W",...% /* 11-20 */
          "2Y","2M","2N","5I","5Q","5X","7I","7Q","7X","6A",...% /* 21-30 */
          "6B","6C","6X","6Z","6S","6L","8L","8Q","8X","2I",...% /* 31-40 */
          "2Q","6I","6Q","3I","3Q","3X","1I","1Q","5A","5B",...% /* 41-50 */
          "5C","9A","9B","9C","9X","1D","5D","5P","5Z","6E",...% /* 51-60 */
          "7D","7P","7Z","8D","8P","4A","4B","4X" ];           % /* 61-68 */
%------------------------------------------------------------------------------------------------------
addpath('MCOSBToolbox'); load('colors.mat');
OSBPath = [Path '\5_FilesAnalysis\' Center '\'];
% Plot used GNSS stations or not
if (Staion)
    addpath('MCOSBToolbox/m_map');
    % Set latitude and longitude range
    lat2=-90;lat1=90;lon1=-180;lon2=180;
    [StationGPS,StationGAL,StationBDS,Sites_InfoGPS,Sites_InfoGAL,Sites_InfoBDS] = GetStation(Path,DOY);
    PlotMultisites(lat2,lat1,lon1,lon2,Sites_InfoGPS.STAInformation,Sites_InfoGAL.STAInformation,Sites_InfoBDS.STAInformation,StationGPS,StationGAL,StationBDS);
end
% AC
if contains(Center,'ECU')
    [Data] = ReadOSBFiles(OSBPath,SatPRN,Signal);
    [OSB,SATG,IDG,SATE,IDE,SATC,IDC] = GetCodeOSB(Data,SatPRN,Signal);
    % Plot Availability of code OSB
    PlotOSBAvailability(OSB,SatPRN);
    % Plot time series of code OSB for each satellite
    PlotOSBSatSeries(OSB,SatPRN,colors);
    % Plot value of code OSB for each satellite
    PlotOSBSatValue(OSB,SatPRN);
    % Plot STD of code OSB for each satellite
    PlotOSBSatSTD(OSB,SatPRN);
end
if contains(Center,'CAS')
    [Data] = ReadOSBFiles(OSBPath,SatPRN,Signal);
    [OSB,SATG,IDG,SATE,IDE,SATC,IDC] = GetCodeOSB(Data,SatPRN,Signal);
    % Plot Availability of code OSB
    PlotOSBAvailability(OSB,SatPRN);
    % Plot time series of code OSB for each satellite
    PlotOSBSatSeries(OSB,SatPRN,colors);
    % Plot value of code OSB for each satellite
    PlotOSBSatValue(OSB,SatPRN);
    % Plot STD of code OSB for each satellite
    PlotOSBSatSTD(OSB,SatPRN);
end
disp('Step four: completing !')