% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
clc;clear;
% Extracting Satellite Plus Receiver (SPR) observables of dual-frequecny GF and triple-frequecny MCB
disp('Step two: extract MCB observables !')
Path = 'E:\MCOSB';
% Satellite cut-off elevation (degree) for OBS
Lim = 15;
% Year, DOY of start, DOY of end
Year = 2025; DOYStart = 1; DOYEnd = 5;
% Processing system. Here, you can select GPS, BDS3, or GAL
% It is noted that there are two pairs of code types for BDS3 or Galileo here
%------------------------------------------------------------------------------------------------------
% SYS = 'GPS';
% SYS = 'BDS3';
% SYS = 'GAL';
SYS = 'ALL';
%------------------------------------------------------------------------------------------------------
addpath('MCOSBToolbox'); load('SATMark.mat');
% Loop
for DOYth = DOYStart:DOYEnd
    % It is noted that the SP3 files should use the 5-min sampling rate
    disp('read SP3 files !');
    SP3Path = [Path,'\2_FilesSP3'];
    if strcmp(SYS, 'GPS') || strcmp(SYS, 'ALL')
        sateG = ReadSP3GPS(SP3Path,DOYth,SATMark,Year);
    end
    if strcmp(SYS, 'BDS3') || strcmp(SYS, 'ALL')
        sateC = ReadSP3BDS(SP3Path,DOYth,SATMark,Year);
    end
    if strcmp(SYS, 'GAL') || strcmp(SYS, 'ALL')
        sateE = ReadSP3GAL(SP3Path,DOYth,SATMark,Year);
    end
    % Read ION for calculating the ionospheric error
    disp('read ION files !'); 
    IONPath = [Path,'\3_FilesION'];
    ION = ReadION(IONPath,DOYth);
    % Calculate the SPR MCB value for each satellite
    disp('Extract SPR MCB observales !')
    if strcmp(SYS, 'GPS') || strcmp(SYS, 'ALL')
        MCBPath      = 'MCB\GPS';
        OBSPath      = [Path,'\OBS\GPS\'];
        STAPath      = [Path,'\STA\GPS\'];
        StrDOYth     = num2str(DOYth,'%03d');
        STAPathDOYth = [STAPath,StrDOYth,'\'];
        load([STAPathDOYth,'STAInformation',StrDOYth,'.mat']);
        GetSPRMCBGPS(OBSPath,MCBPath,sateG,STAInformation,Lim*pi/180,SATMark,ION,Year,StrDOYth);
    end
    if strcmp(SYS, 'BDS3') || strcmp(SYS, 'ALL')
        MCBPath      = 'MCB\BDS3';
        OBSPath      = [Path,'\OBS\BDS3\'];
        STAPath      = [Path,'\STA\BDS3\'];
        StrDOYth     = num2str(DOYth,'%03d');
        STAPathDOYth = [STAPath,StrDOYth,'\'];
        load([STAPathDOYth,'STAInformation',StrDOYth,'.mat']);
        B1C = 'BDSC1X'; B1I = 'BDSC2I'; B3I = 'BDSC6I'; B2B = 'BDSC7Z'; B2A = 'BDSC5X';
        GetSPRMCBBDS(OBSPath,MCBPath,sateC,STAInformation,Lim*pi/180,SATMark,ION,B1C,B1I,B2A,B3I,B2B,Year,StrDOYth);
        B1C = 'BDSC1P'; B1I = 'BDSC2I'; B3I = 'BDSC6I'; B2B = 'BDSC7D'; B2A = 'BDSC5P';
        GetSPRMCBBDS(OBSPath,MCBPath,sateC,STAInformation,Lim*pi/180,SATMark,ION,B1C,B1I,B2A,B3I,B2B,Year,StrDOYth);
    end
    if strcmp(SYS, 'GAL') || strcmp(SYS, 'ALL')
        MCBPath      = 'MCB\GAL';
        OBSPath      = [Path,'\OBS\GAL\'];
        STAPath      = [Path,'\STA\GAL\'];
        StrDOYth     = num2str(DOYth,'%03d');
        STAPathDOYth = [STAPath,StrDOYth,'\'];
        load([STAPathDOYth,'STAInformation',StrDOYth,'.mat']);
        E1='GALC1C'; E6='GALC6C'; E5B='GALC7Q'; E5='GALC8Q'; E5A='GALC5Q';
        GetSPRMCBGAL(OBSPath,MCBPath,sateE,STAInformation,Lim*pi/180,SATMark,ION,E1,E6,E5B,E5,E5A,Year,StrDOYth);
        E1='GALC1X'; E6='GALC6X'; E5B='GALC7X'; E5='GALC8X'; E5A='GALC5X';
        GetSPRMCBGAL(OBSPath,MCBPath,sateE,STAInformation,Lim*pi/180,SATMark,ION,E1,E6,E5B,E5,E5A,Year,StrDOYth);
    end
end
disp('Step two: completing !')