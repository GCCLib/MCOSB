% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
clc;clear;
disp('Step one: read GNSS rinex files !')
% Filepath of OBS
Path='E:\MCOSB';
% Processing system. Here, you can select GPS, BDS3, GAL, or ALL
% SYS = 'GPS';
% SYS = 'BDS3';
% SYS = 'GAL';
SYS = 'ALL';
%------------------------------------------------------------------------------------------------------
addpath('MCOSBToolbox');
OBSPath = [Path,'\1_FilesOBS\'];
list=dir(OBSPath);len=length(list);
% Loop each day
for i = 3:len
    % Output filepath of OBS
    if strcmp(SYS, 'GPS') || strcmp(SYS, 'ALL')
        OutputPathG='OBS\GPS';
    end
    if strcmp(SYS, 'BDS3') || strcmp(SYS, 'ALL')
        OutputPathC='OBS\BDS3';
    end
    if strcmp(SYS, 'GAL') || strcmp(SYS, 'ALL')
        OutputPathE='OBS\GAL';
    end
    % Output filepath of station information
    if strcmp(SYS, 'GPS') || strcmp(SYS, 'ALL')
        STAPathG='STA\GPS';
    end
    if strcmp(SYS, 'BDS3') || strcmp(SYS, 'ALL')
        STAPathC='STA\BDS3';
    end
    if strcmp(SYS, 'GAL') || strcmp(SYS, 'ALL')
        STAPathE='STA\GAL';
    end
    % Filepath
    StrDOY  = num2str(list(i).name,'%03d');
    PathDOY = [OBSPath,StrDOY,'\'];
    if strcmp(SYS, 'GPS') || strcmp(SYS, 'ALL')
        ReadRinexGPS(PathDOY,OutputPathG,STAPathG);
    end
    if strcmp(SYS, 'BDS3') || strcmp(SYS, 'ALL')
        ReadRinexBDS(PathDOY,OutputPathC,STAPathC);
    end
    if strcmp(SYS, 'GAL') || strcmp(SYS, 'ALL')
        ReadRinexGAL(PathDOY,OutputPathE,STAPathE);
    end
end
disp('Step one: completing !')