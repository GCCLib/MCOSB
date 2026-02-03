% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
clc;clear;
disp('Step three: estimate multi-GNSS and multi-frequency code OSB products !')
Path = 'E:\MCOSB';
% Year, DOY of start, DOY of end
Year = 2025; DOYStart = 1; DOYEnd = 5;
%------------------------------------------------------------------------------------------------------
% Processing system. Here, you can select GPS, BDS3, GAL, or ALL
% SYS = 'GPS';
% SYS = 'BDS3';
% SYS = 'GAL';
SYS = 'ALL';
%------------------------------------------------------------------------------------------------------
addpath('MCOSBToolbox');
NumFre = 5;
%------------------------------------------------------------------------------------------------------
if strcmp(SYS, 'GPS') || strcmp(SYS, 'ALL')
    GPSSatOSBC5Q = [];
    GPSSatOSBC5X = [];
    % Loop each day
    for DOYTh = DOYStart : DOYEnd
        %------------------------------------------------------------------------------------------------------
        StrDOY = num2str(DOYTh,'%03d');
        InPath = ['MCB\GPS\',StrDOY,'\','GPSC5Q\'];
        list_gps=dir([InPath '/*.mat']);
        % Number of used receivers
        NumRec=length(list_gps);
        % Check the number of each satellite's observations
        [GSatTypeNum,GsatC1W,GsatC2W,GsatC5Q,GsatC1C,GsatC2L,GPRNC1C,GPRNC1W,GPRNC2W,GPRNC2L,GPRNC5Q] = GetSatTypeNumC5Q(NumRec,InPath,list_gps);
        [GRecTypeNum] = GetRecTypeNumC5Q(InPath);
        % Total number of estimated parameters
        NumPar    = sum(GSatTypeNum) + sum(GRecTypeNum);
        NumRecPar = sum(GRecTypeNum);
        NBB = zeros(NumPar,NumPar);
        W   = zeros(NumPar,1);
        %Loop each receiver
        for j = 1:NumRec
            load([InPath '/' list_gps(j).name],'-mat');
            %Delete zero colum and get combination value again
            [C1WC2W,C1WC2WC5Q,C1WC1C,C2WC2L] = GetCombinationC5Q(GPSSPRMCB,GsatC1W,GsatC5Q,GsatC1C,GsatC2L);
            %Get the obs equation
            [sN,sl,P]=GetMatrixC5Q(GSatTypeNum,GRecTypeNum,C1WC2W,C1WC2WC5Q,C1WC1C,C2WC2L,j,NumPar,NumRecPar);
            NBB=NBB+sN'*P*sN;
            W=W+sN'*P*sl;
            clear GPSSPRMCB C1WC2W C1WC2WC5Q C1WC1C C2WC2L;
            disp(['------------> [',num2str(DOYTh),'] [',num2str(j),'/',num2str(NumRec),']  ',' % GPSC5Q has been solved!!']);
        end
        [GPSOSB,SatC1C,SatC1W,SatC2W,SatC2L] = EstimateOSBGPSC5Q(GSatTypeNum,GRecTypeNum,NumPar,NumRecPar,NumRec,NBB,W);
        [GPSSatOSB] = GetSatOSBGPSC5Q(SatC1C,SatC1W,SatC2W,SatC2L,GPRNC1C,GPRNC1W,GPRNC2W,GPRNC2L,GPRNC5Q,GPSOSB);
        GPSSatOSBC5Q(DOYTh-DOYStart+1).OSB=GPSSatOSB;
        clear NBB W GSatTypeNum GRecTypeNum GPSSatOSB list_gps;
        %------------------------------------------------------------------------------------------------------
        InPath = ['MCB\GPS\',StrDOY,'\','GPSC5X\'];
        list_gps=dir([InPath '/*.mat']);
        % Number of used receivers
        NumRec=length(list_gps);
        % Check the number of each satellite's observations
        [GSatTypeNum,GsatC1W,GsatC2W,GsatC5X,GsatC1C,GsatC2X,GPRNC1C,GPRNC1W,GPRNC2W,GPRNC2X,GPRNC5X] = GetSatTypeNumC5X(NumRec,InPath,list_gps);
        [GRecTypeNum] = GetRecTypeNumC5X(InPath); 
        % Total number of estimated parameters
        NumPar    = sum(GSatTypeNum) + sum(GRecTypeNum);
        NumRecPar = sum(GRecTypeNum);
        NBB = zeros(NumPar,NumPar);
        W   = zeros(NumPar,1);
        % Loop each receiver
        for j = 1:NumRec
            load([InPath '/' list_gps(j).name],'-mat');
            %Delete zero colum and get combination value again
            [C1WC2W,C1WC2WC5X,C1WC1C,C2WC2X] = GetCombinationC5X(GPSSPRMCB,GsatC1W,GsatC5X,GsatC1C,GsatC2X);
            %Get the obs equation
            [sN,sl,P]=GetMatrixC5X(GSatTypeNum,GRecTypeNum,C1WC2W,C1WC2WC5X,C1WC1C,C2WC2X,j,NumPar,NumRecPar);
            NBB=NBB+sN'*P*sN;
            W=W+sN'*P*sl;
            clear GPSSPRMCB C1WC2W C1WC2WC5X C1WC1C C2WC2X;
            disp(['------------> [',num2str(DOYTh),'] [',num2str(j),'/',num2str(NumRec),']  ',' % GPSC5X has been solved!']);
        end
        [GPSOSB,SatC1C,SatC1W,SatC2W,SatC2X] = EstimateOSBGPSC5X(GSatTypeNum,GRecTypeNum,NumPar,NumRecPar,NumRec,NBB,W);
        [GPSSatOSB] = GetSatOSBGPSC5X(SatC1C,SatC1W,SatC2W,SatC2X,GPRNC1C,GPRNC1W,GPRNC2W,GPRNC2X,GPRNC5X,GPSOSB);
        GPSSatOSBC5X(DOYTh-DOYStart+1).OSB=GPSSatOSB;
        clear NBB W GSatTypeNum GRecTypeNum GPSSatOSB list_gps;
    end
    if exist('OSB\GPS','dir')==0
        mkdir('OSB\GPS');
    end
    save(['OSB\GPS','\GPSSatOSBC5Q.mat'], 'GPSSatOSBC5Q', '-mat');
    save(['OSB\GPS','\GPSSatOSBC5X.mat'], 'GPSSatOSBC5X', '-mat');    
end
%------------------------------------------------------------------------------------------------------
if strcmp(SYS, 'BDS3') || strcmp(SYS, 'ALL')
    BDSSatOSBC1P = [];
    BDSSatOSBC1X = [];
    % Loop each day
    for DOYTh = DOYStart:DOYEnd
        StrDOY = num2str(DOYTh,'%03d');
        InPath = ['MCB\BDS3\',StrDOY,'\','BDSC1P\'];
        list_bds=dir([InPath '/*.mat']);
        % Number of used receivers
        NumRec=length(list_bds);      
        % Check the number of each satellite's observations
        [BSatTypeNum,Bsat,BPRN] = GetSatTypeNumBDS(NumRec,InPath,list_bds);
        NumPar = (BSatTypeNum + NumRec) * NumFre;
        NBB    = zeros(NumPar,NumPar); 
        W      = zeros(NumPar,1);
        for j = 1:NumRec
            load([InPath '/' list_bds(j).name],'-mat');
            % Delete zero colum, i.e., without the observed satellites
            [TempBDSSPRMCB] = GetCombinationBDS(Bsat,BDSSPRMCB);
            % Get the obs equation
            [sN,sl,P]=GetMatrixBDS(TempBDSSPRMCB,NumRec,BSatTypeNum,j,NumPar,NumFre);
            NBB=NBB+sN'*P*sN;
            W=W+sN'*P*sl;
            clear TempBDSSPRMCB BDSSPRMCB;
            disp(['------------> [',num2str(DOYTh),'] [',num2str(j),'/',num2str(NumRec),']  ',' % BDSC1P has been solved!']);
        end
        [BDSOSB] = EstimateOSBBDS(BSatTypeNum,NumPar,NumRec,NBB,W,NumFre);
        [BDSSatOSB] = GetSatOSBBDSC1P(BSatTypeNum,BPRN,BDSOSB);
        BDSSatOSBC1P(DOYTh-DOYStart+1).OSB=BDSSatOSB;
        clear NBB W BDSSatOSB list_bds;
        %------------------------------------------------------------------------------------------------------
        InPath = ['MCB\BDS3\',StrDOY,'\','BDSC1X\'];
        list_bds=dir([InPath '/*.mat']);
        % Number of used receivers
        NumRec=length(list_bds);      
        % Check the number of each satellite's observations
        [BSatTypeNum,Bsat,BPRN] = GetSatTypeNumBDS(NumRec,InPath,list_bds);
        NumPar = (BSatTypeNum + NumRec) * NumFre;
        NBB = zeros(NumPar,NumPar);
        W   = zeros(NumPar,1);
        for j = 1:NumRec
            load([InPath '/' list_bds(j).name],'-mat');
            [TempBDSSPRMCB] = GetCombinationBDS(Bsat,BDSSPRMCB);
            % Get the obs equation
            [sN,sl,P]=GetMatrixBDS(TempBDSSPRMCB,NumRec,BSatTypeNum,j,NumPar,NumFre);
            NBB=NBB+sN'*P*sN;
            W=W+sN'*P*sl;
            clear TempBDSSPRMCB BDSSPRMCB;
            disp(['------------> [',num2str(DOYTh),'] [',num2str(j),'/',num2str(NumRec),']  ',' % BDSC1X has been solved!']);
        end
        [BDSOSB]    = EstimateOSBBDS(BSatTypeNum,NumPar,NumRec,NBB,W,NumFre);
        [BDSSatOSB] = GetSatOSBBDSC1X(BSatTypeNum,BPRN,BDSOSB);
        BDSSatOSBC1X(DOYTh-DOYStart+1).OSB=BDSSatOSB;
        clear NBB W BDSSatOSB list_bds;
    end
    if exist('OSB\BDS3','dir')==0
        mkdir('OSB\BDS3');
    end
    save(['OSB\BDS3','\BDSSatOSBC1P.mat'], 'BDSSatOSBC1P', '-mat');
    save(['OSB\BDS3','\BDSSatOSBC1X.mat'], 'BDSSatOSBC1X', '-mat');    
end
%------------------------------------------------------------------------------------------------------
if strcmp(SYS, 'GAL') || strcmp(SYS, 'ALL')
    GALSatOSBC1C = [];
    GALSatOSBC1X = [];
    % Loop each day
    for DOYTh = DOYStart:DOYEnd
        StrDOY = num2str(DOYTh,'%03d');
        InPath = ['MCB\GAL\',StrDOY,'\','GALC1C\'];
        list_gal=dir([InPath '/*.mat']);
        % Number of used receivers
        NumRec=length(list_gal);      
        % Check the number of each satellite's observations
        [ESatTypeNum,Esat,EPRN] = GetSatTypeNumGAL(NumRec,InPath,list_gal);
        NumPar = (ESatTypeNum + NumRec) * NumFre;
        NBB = zeros(NumPar,NumPar);
        W   = zeros(NumPar,1);
        for j = 1:NumRec
            load([InPath '/' list_gal(j).name],'-mat');
            % Delete zero colum, i.e., without the observed satellites
            [TempGALSPRMCB] = GetCombinationGAL(Esat,GALSPRMCB);
            % Get the obs equation
            [sN,sl,P]=GetMatrixGAL(TempGALSPRMCB,NumRec,ESatTypeNum,j,NumPar,NumFre);
            NBB=NBB+sN'*P*sN;
            W=W+sN'*P*sl;
            clear TempGALSPRMCB GALSPRMCB;
            disp(['------------> [',num2str(DOYTh),'] [',num2str(j),'/',num2str(NumRec),']  ',' % GALC1C has been solved!']);
        end
        [GALOSB] = EstimateOSBGAL(ESatTypeNum,NumPar,NumRec,NBB,W,NumFre);
        [GALSatOSB] = GetSatOSBGALC1C(ESatTypeNum,EPRN,GALOSB);
        GALSatOSBC1C(DOYTh-DOYStart+1).OSB=GALSatOSB;
        clear NBB W GALSatOSB list_gal;
        %------------------------------------------------------------------------------------------------------
        InPath = ['MCB\GAL\',StrDOY,'\','GALC1X\'];
        list_gal=dir([InPath '/*.mat']);
        % Number of used receivers
        NumRec=length(list_gal);      
        % Check the number of each satellite's observations
        [ESatTypeNum,Esat,EPRN] = GetSatTypeNumGAL(NumRec,InPath,list_gal);
        NumPar = (ESatTypeNum + NumRec) * NumFre;
        NBB = zeros(NumPar,NumPar);
        W   = zeros(NumPar,1);
        for j = 1:NumRec
            load([InPath '/' list_gal(j).name],'-mat');
            % Delete zero colum, i.e., without the observed satellites
            [TempGALSPRMCB] = GetCombinationGAL(Esat,GALSPRMCB);
            % Get the obs equation
            [sN,sl,P]=GetMatrixGAL(TempGALSPRMCB,NumRec,ESatTypeNum,j,NumPar,NumFre);
            NBB=NBB+sN'*P*sN;
            W=W+sN'*P*sl;
            clear TempGALSPRMCB GALSPRMCB;
            disp(['------------> [',num2str(DOYTh),'] [',num2str(j),'/',num2str(NumRec),']  ',' % GALC1X has been solved!']);
        end
        [GALOSB] = EstimateOSBGAL(ESatTypeNum,NumPar,NumRec,NBB,W,NumFre);
        [GALSatOSB] = GetSatOSBGALC1X(ESatTypeNum,EPRN,GALOSB);
        GALSatOSBC1X(DOYTh-DOYStart+1).OSB=GALSatOSB;
        clear NBB W GALSatOSB list_gal;
    end
    if exist('OSB\GAL','dir')==0
        mkdir('OSB\GAL');
    end
    save(['OSB\GAL','\GALSatOSBC1C.mat'], 'GALSatOSBC1C', '-mat');
    save(['OSB\GAL','\GALSatOSBC1X.mat'], 'GALSatOSBC1X', '-mat');
end
%------------------------------------------------------------------------------------------------------
% Generate code OSB SINEX files
% Loop
for DOYTh = DOYStart:DOYEnd
    %SINEX-file head
    StrDOY = num2str(DOYTh,'%03d');
    fid_read  = fopen([Path '\4_FilesOSB\FileHead.BIA'], 'r');
    % Depending on the users
    FileName = ['\4_FilesOSB\ECU0MGXFIN_' num2str(Year), StrDOY, '0000_01D_01D_OSB.BIA'];
    fid_write = fopen([Path FileName], 'w');
    if fid_read == -1 || fid_write == -1
        if fid_read  ~= -1,  fclose(fid_read);  end
        if fid_write ~= -1,  fclose(fid_write); end
        error('error: no file');
    end
    while ~feof(fid_read)
        line = fgetl(fid_read);
        fprintf(fid_write, '%s\n', line);
    end 
    %SINEX-file body
    index = DOYTh - DOYStart + 1;
    fid_write = WriteOSBGPS(GPSSatOSBC5Q,GPSSatOSBC5X,fid_write,index,DOYTh,Year);
    fid_write = WriteOSBGAL(GALSatOSBC1C,GALSatOSBC1X,fid_write,index,DOYTh,Year);
    fid_write = WriteOSBBDS(BDSSatOSBC1P,BDSSatOSBC1X,fid_write,index,DOYTh,Year);
    %Close file
    fclose(fid_read); fclose(fid_write);
    disp(['>----- ' num2str(Year), num2str(DOYTh), ' % SINEX file has been written！']);
end
disp('Step three: completing !')