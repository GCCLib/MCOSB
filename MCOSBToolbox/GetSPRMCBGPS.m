% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [] = GetSPRMCBGPS(OBSPath,MCBPath,sate,STAInformation,lim,SATMark,ION,year,current_doy)
    % Station position
    Coor=STAInformation.coor;
    % Station name
    stations=STAInformation.name;
    % load saved OBS
    list_obs=dir([OBSPath,'/',current_doy,'/*.mat']);
    len=length(list_obs);
    DOY = str2double(current_doy);
    jd=doy2jd_GT(year,DOY);
    [~,sow,~] = jd2gps_GT(jd);
    % Loop each station
    for i=1:len
        load([OBSPath,'/',current_doy,'/',list_obs(i).name],'-mat');
        site=list_obs(i).name(1:4);
        if contains(site,'ACRG')
            debug = 1;
        end
        doy=list_obs(i).name(5:7);
        % Which station
        index=find(strcmpi(site,stations), 1);
        sx=Coor(index,1);
        sy=Coor(index,2);
        sz=Coor(index,3);
        % Exclude sites without receivers coordinates
        if isempty(sx)|| isempty(sy)|| isempty(sz)
            continue;
        end
        if sx==0 || sy==0 || sz==0
            continue;
        end
        % cut obs
        obs=cutobs(sate,sx,sy,sz,obs,lim,SATMark);
        [GPSSPRMCB]=GPS_SPRMCB(sate,sx,sy,sz,obs,lim,ION,sow);
        % Quality control, when nsat smaller than 5
        no_zero = ~all(GPSSPRMCB == 0, 1);      
        [~,col] = size(GPSSPRMCB(:, no_zero));
        if col < 4
            continue;
        end
        clear no_zero;
        % Different OBS types for GPS
        %C5X
        if sum(GPSSPRMCB(3,:))~=0
            if exist([MCBPath '/' doy '/GPSC5X/'],'dir')==0
                mkdir([MCBPath '/' doy '/GPSC5X/']);
            end
            filenameSPRMCB = [MCBPath '/' doy '/GPSC5X/' site doy 'MCB.mat'];
            save(filenameSPRMCB,'GPSSPRMCB','-mat');
        %C5Q
        elseif sum(GPSSPRMCB(2,:))~=0
            if exist([MCBPath '/' doy '/GPSC5Q/'],'dir')==0
                mkdir([MCBPath '/' doy '/GPSC5Q/']);
            end
            filenameSPRMCB = [MCBPath '/' doy '/GPSC5Q/' site doy 'MCB.mat'];
            save(filenameSPRMCB,'GPSSPRMCB','-mat');
        end
        disp(['-------- [',num2str(DOY),'] [ ',num2str(i),'/',num2str(len),' ] GPSSPRMCB is processed !']);
    end
end
%% ----------------subfunction-----------------
function obs=cutobs(sate,sx,sy,sz,obs,lim,sate_mark)
% get the observations under specified cut-off angle
% INPUT: 
%     sate: precise coordinates of the satellites
%     sx: X coordinate of the station
%     sy: Y coordinate of the station
%     sz: Z coordinate of the station
%     obs: original observation structs
%     lim: cut-off angle
%     sate_mark: satellite status identification
% OUTPUT：
%     obs: updated observation structs
fields=fieldnames(obs);
% cut GPS obs
if ~isnan(find(strcmp(fields, 'GPSC1C' )))
    gpsline=size(obs.GPSC1C,1);
    gpssate=size(sate_mark.gps,2);
    if gpsline<2880
        obs.GPSC1C(gpsline+1:2880,:)=0;
    end
    gpsl=size(obs.GPSC1C,2);
    if gpsl<gpssate
        obs.GPSC1C(:,gpsl+1:gpssate)=0;
    end    
    gpsdelete=find(sate_mark.gps==0);
    gpsnum=size(sate.gpsx,2);    
    % if gpsnum<size(obs.GPSC1C,2)
    %     if ~isnan(gpsdelete)
    %         k=length(gpsdelete);
    %         for i=k:-1:1
    %             obs.GPSC1C(:,gpsdelete(i))=[];
    %         end
    %     end
    %     obs.GPSC1C=obs.GPSC1C(:,1:gpsnum);
    % end
    obs.GPSC1C(isnan(obs.GPSC1C))=0;    
    gpsx=sate.gpsx;gpsy=sate.gpsy;gpsz=sate.gpsz;
    for i=1:gpsnum
        for j=1:2880
            if obs.GPSC1C(j,i)==0
                obs.GPSC1C(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,gpsx(j,i),gpsy(j,i),gpsz(j,i));
            if el<lim
                obs.GPSC1C(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'GPSC1W' )))
    gpsline=size(obs.GPSC1W,1);
    gpssate=size(sate_mark.gps,2);
    if gpsline<2880
        obs.GPSC1W(gpsline+1:2880,:)=0;
    end
    gpsl=size(obs.GPSC1W,2);
    if gpsl<gpssate
        obs.GPSC1W(:,gpsl+1:gpssate)=0;
    end    
    gpsdelete=find(sate_mark.gps==0);
    gpsnum=size(sate.gpsx,2);    
    % if gpsnum<size(obs.GPSC1W,2)
    %     if ~isnan(gpsdelete)
    %         k=length(gpsdelete);
    %         for i=k:-1:1
    %             obs.GPSC1W(:,gpsdelete(i))=[];
    %         end
    %     end
    %     obs.GPSC1W=obs.GPSC1W(:,1:gpsnum);
    % end
    obs.GPSC1W(isnan(obs.GPSC1W))=0;    
    gpsx=sate.gpsx;gpsy=sate.gpsy;gpsz=sate.gpsz;
    for i=1:gpsnum
        for j=1:2880
            if obs.GPSC1W(j,i)==0
                obs.GPSC1W(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,gpsx(j,i),gpsy(j,i),gpsz(j,i));
            if el<lim
                obs.GPSC1W(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'GPSC2W' )))
    gpsline=size(obs.GPSC2W,1);
    gpssate=size(sate_mark.gps,2);
    if gpsline<2880
        obs.GPSC2W(gpsline+1:2880,:)=0;
    end
    gpsl=size(obs.GPSC2W,2);
    if gpsl<gpssate
        obs.GPSC2W(:,gpsl+1:gpssate)=0;
    end    
    gpsdelete=find(sate_mark.gps==0);
    gpsnum=size(sate.gpsx,2);    
    % if gpsnum<size(obs.GPSC2W,2)
    %     if ~isnan(gpsdelete)
    %         k=length(gpsdelete);
    %         for i=k:-1:1
    %             obs.GPSC2W(:,gpsdelete(i))=[];
    %         end
    %     end
    %     obs.GPSC2W=obs.GPSC2W(:,1:gpsnum);
    % end
    obs.GPSC2W(isnan(obs.GPSC2W))=0;    
    gpsx=sate.gpsx;gpsy=sate.gpsy;gpsz=sate.gpsz;
    for i=1:gpsnum
        for j=1:2880
            if obs.GPSC2W(j,i)==0
                obs.GPSC2W(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,gpsx(j,i),gpsy(j,i),gpsz(j,i));
            if el<lim
                obs.GPSC2W(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'GPSC2L' )))
    gpsline=size(obs.GPSC2L,1);
    gpssate=size(sate_mark.gps,2);
    if gpsline<2880
        obs.GPSC2L(gpsline+1:2880,:)=0;
    end
    gpsl=size(obs.GPSC2L,2);
    if gpsl<gpssate
        obs.GPSC2L(:,gpsl+1:gpssate)=0;
    end    
    gpsdelete=find(sate_mark.gps==0);
    gpsnum=size(sate.gpsx,2);    
    % if gpsnum<size(obs.GPSC2L,2)
    %     if ~isnan(gpsdelete)
    %         k=length(gpsdelete);
    %         for i=k:-1:1
    %             obs.GPSC2L(:,gpsdelete(i))=[];
    %         end
    %     end
    %     obs.GPSC2L=obs.GPSC2L(:,1:gpsnum);
    % end
    obs.GPSC2L(isnan(obs.GPSC2L))=0;    
    gpsx=sate.gpsx;gpsy=sate.gpsy;gpsz=sate.gpsz;
    for i=1:gpsnum
        for j=1:2880
            if obs.GPSC2L(j,i)==0
                obs.GPSC2L(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,gpsx(j,i),gpsy(j,i),gpsz(j,i));
            if el<lim
                obs.GPSC2L(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'GPSC2X' )))
    gpsline=size(obs.GPSC2X,1);
    gpssate=size(sate_mark.gps,2);
    if gpsline<2880
        obs.GPSC2X(gpsline+1:2880,:)=0;
    end
    gpsl=size(obs.GPSC2X,2);
    if gpsl<gpssate
        obs.GPSC2X(:,gpsl+1:gpssate)=0;
    end    
    gpsdelete=find(sate_mark.gps==0);
    gpsnum=size(sate.gpsx,2);    
    % if gpsnum<size(obs.GPSC2X,2)
    %     if ~isnan(gpsdelete)
    %         k=length(gpsdelete);
    %         for i=k:-1:1
    %             obs.GPSC2X(:,gpsdelete(i))=[];
    %         end
    %     end
    %     obs.GPSC2X=obs.GPSC2X(:,1:gpsnum);
    % end
    obs.GPSC2X(isnan(obs.GPSC2X))=0;    
    gpsx=sate.gpsx;gpsy=sate.gpsy;gpsz=sate.gpsz;
    for i=1:gpsnum
        for j=1:2880
            if obs.GPSC2X(j,i)==0
                obs.GPSC2X(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,gpsx(j,i),gpsy(j,i),gpsz(j,i));
            if el<lim
                obs.GPSC2X(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'GPSC5Q' )))
    gpsline=size(obs.GPSC5Q,1);
    gpssate=size(sate_mark.gps,2);
    if gpsline<2880
        obs.GPSC5Q(gpsline+1:2880,:)=0;
    end
    gpsl=size(obs.GPSC5Q,2);
    if gpsl<gpssate
        obs.GPSC5Q(:,gpsl+1:gpssate)=0;
    end    
    gpsdelete=find(sate_mark.gps==0);
    gpsnum=size(sate.gpsx,2);    
    % if gpsnum<size(obs.GPSC5Q,2)
    %     if ~isnan(gpsdelete)
    %         k=length(gpsdelete);
    %         for i=k:-1:1
    %             obs.GPSC5Q(:,gpsdelete(i))=[];
    %         end
    %     end
    %     obs.GPSC5Q=obs.GPSC5Q(:,1:gpsnum);
    % end
    obs.GPSC5Q(isnan(obs.GPSC5Q))=0;    
    gpsx=sate.gpsx;gpsy=sate.gpsy;gpsz=sate.gpsz;
    for i=1:gpsnum
        for j=1:2880
            if obs.GPSC5Q(j,i)==0
                obs.GPSC5Q(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,gpsx(j,i),gpsy(j,i),gpsz(j,i));
            if el<lim
                obs.GPSC5Q(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'GPSC5X' )))
    gpsline=size(obs.GPSC5X,1);
    gpssate=size(sate_mark.gps,2);
    if gpsline<2880
        obs.GPSC5X(gpsline+1:2880,:)=0;
    end
    gpsl=size(obs.GPSC5X,2);
    if gpsl<gpssate
        obs.GPSC5X(:,gpsl+1:gpssate)=0;
    end    
    gpsdelete=find(sate_mark.gps==0);
    gpsnum=size(sate.gpsx,2);    
    % if gpsnum<size(obs.GPSC5X,2)
    %     if ~isnan(gpsdelete)
    %         k=length(gpsdelete);
    %         for i=k:-1:1
    %             obs.GPSC5X(:,gpsdelete(i))=[];
    %         end
    %     end
    %     obs.GPSC5X=obs.GPSC5X(:,1:gpsnum);
    % end
    obs.GPSC5X(isnan(obs.GPSC5X))=0;    
    gpsx=sate.gpsx;gpsy=sate.gpsy;gpsz=sate.gpsz;
    for i=1:gpsnum
        for j=1:2880
            if obs.GPSC5X(j,i)==0
                obs.GPSC5X(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,gpsx(j,i),gpsy(j,i),gpsz(j,i));
            if el<lim
                obs.GPSC5X(j,i)=0;
                continue;
            end
        end
    end
end
end
%% ---------------------------subfunction-------------------------
function [GPSSPRMCB]=GPS_SPRMCB(sate,sx,sy,sz,obs,lim,ION,sow)
    size2=size(obs.GPSC1W,2);
    % C1WC2W/C1WC2WC5Q/C1WC2WC5X/C1WC1C/C2WC2L/C2WC2X
    GPSSPRMCB=[];
    SPRC1WC2W     = zeros(2880,32);
    SPRC1WC2WC5Q  = zeros(2880,32);
    SPRC1WC2WC5X  = zeros(2880,32);
    SPRC1WC1C     = zeros(2880,32);
    SPRC2WC2L     = zeros(2880,32);
    SPRC2WC2X     = zeros(2880,32);
    fields=fieldnames(obs);
    resultString = strjoin(fields);
    % Constants
    CLIGHT = 299792458;   
    FL1  = 1.57542E9;
    FL2  = 1.22760E9;
    FL5 = 1.17645E9;
    WL1 = CLIGHT/FL1; 
    WL2 = CLIGHT/FL2;
    WL5 = CLIGHT/FL5;
    A32 = SQR(WL5) - SQR(WL2);
    A13 = SQR(WL1) - SQR(WL5);
    A21 = SQR(WL2) - SQR(WL1);    
    [B,L,H]=XYZtoBLH(sx,sy,sz);
    gpsx=sate.gpsx;gpsy=sate.gpsy;gpsz=sate.gpsz;
    % i is PRN number
    for i=1:size2 
        for j=1:2880
            OBSC1W = obs.GPSC1W(j, i);
            if OBSC1W == 0
                continue;
            end
            [E,A]=Get_EA(sx,sy,sz,gpsx(j,i),gpsy(j,i),gpsz(j,i));
            if E<lim
                continue;
            end
            % Calculate the ionospheric delay
            OBSC1W = obs.GPSC1W(j, i);
            OBSC2W = obs.GPSC2W(j, i);
            CodeDistance = [OBSC1W;OBSC2W];
            CodeDistance(CodeDistance(:,1)==0,:)=[];
            tau = mean(CodeDistance)/CLIGHT;
            Ttr = sow+(30*(j-1))-tau;
            % GIM
            elev = E*180/pi;
            mf = 'COSZ';
            az = A*180/pi;
            pos_WGS84.lat=B;
            pos_WGS84.lon=L;
            pos_WGS84.h  =H;
            radius = 6371;
            hgt = [450 450 0];
            interpol = 'Consecutive Rotated Maps';
            [mappingf, Lat_IPP, Lon_IPP] = iono_mf(elev, mf, pos_WGS84, az, radius, hgt);
            vtec = iono_gims(Lat_IPP, Lon_IPP, Ttr, ION.ionex, interpol); % interpolate VTEC
            iono(1) = mappingf * 40.3e16/FL1^2 * vtec; %delta_iono [m]
            iono(2) = mappingf * 40.3e16/FL2^2 * vtec;      
            % Calculate SPRMCB
            if contains(resultString,'GPSC1C'),OBSC1C = obs.GPSC1C(j, i); else, OBSC1C = 0; end
            if contains(resultString,'GPSC1W'),OBSC1W = obs.GPSC1W(j, i); else, OBSC1W = 0; end
            if contains(resultString,'GPSC2W'),OBSC2W = obs.GPSC2W(j, i); else, OBSC2W = 0; end
            if contains(resultString,'GPSC2L'),OBSC2L = obs.GPSC2L(j, i); else, OBSC2L = 0; end
            if contains(resultString,'GPSC2X'),OBSC2X = obs.GPSC2X(j, i); else, OBSC2X = 0; end
            if contains(resultString,'GPSC5Q'),OBSC5Q = obs.GPSC5Q(j, i); else, OBSC5Q = 0; end
            if contains(resultString,'GPSC5X'),OBSC5X = obs.GPSC5X(j, i); else, OBSC5X = 0; end    
            % Reference OBS types
            SPRC1WC2W(j,i)     = (OBSC1W - OBSC2W) - (iono(1) - iono(2));
            % GFIF
            SPRC1WC2WC5Q(j,i)  = A32 * OBSC1W + A13 * OBSC2W + A21 * OBSC5Q;
            SPRC1WC2WC5X(j,i)  = A32 * OBSC1W + A13 * OBSC2W + A21 * OBSC5X;
            SPRC1WC1C(j,i)     = OBSC1W - OBSC1C;
            SPRC2WC2L(j,i)     = OBSC2W - OBSC2L;
            SPRC2WC2X(j,i)     = OBSC2W - OBSC2X;
            % Quality control
            if OBSC1W == 0 || OBSC2W == 0 || vtec == 0
                SPRC1WC2W(j,i)=0.0;
            end
            if OBSC1W == 0 || OBSC2W == 0 || OBSC5Q == 0
                SPRC1WC2WC5Q(j,i)=0.0;
            end
            if OBSC1W == 0 || OBSC2W == 0 || OBSC5X == 0
                SPRC1WC2WC5X(j,i)=0.0;
            end
            if OBSC1W == 0 || OBSC1C == 0 
                SPRC1WC1C(j,i)=0.0;
            end
            if OBSC2W == 0 || OBSC2L == 0 
                SPRC2WC2L(j,i)=0.0;
            end
            if OBSC2W == 0 || OBSC2X == 0 
                SPRC2WC2X(j,i)=0.0;
            end      
        end
    end  
    % Quality control
    for ii = 1:size2
        [C1WC2W,STDC1WC2W]         = RemoveOutlier(SPRC1WC2W,ii);
        [C1WC2WC5Q,STDC1WC2WC5Q]   = RemoveOutlier(SPRC1WC2WC5Q,ii);
        [C1WC2WC5X,STDC1WC2WC5X]   = RemoveOutlier(SPRC1WC2WC5X,ii);
        [C1WC1C,STDC1WC1C]         = RemoveOutlier(SPRC1WC1C,ii);
        [C2WC2L,STDC2WC2L]         = RemoveOutlier(SPRC2WC2L,ii);
        [C2WC2X,STDC2WC2X]         = RemoveOutlier(SPRC2WC2X,ii);
        % No reference OBS types
        if C1WC2W == 0
            C1WC2WC5Q = 0; C1WC2WC5X = 0; C1WC1C = 0; C2WC2L = 0; C2WC2X = 0;
            STDC1WC2W=0; STDC1WC2WC5Q=0; STDC1WC2WC5X=0; STDC1WC1C=0; STDC2WC2L=0; STDC2WC2X=0;
        end
        % Store
        SPRSAT=[C1WC2W;C1WC2WC5Q;C1WC2WC5X;C1WC1C;C2WC2L;C2WC2X;STDC1WC2W;STDC1WC2WC5Q;STDC1WC2WC5X;STDC1WC1C;STDC2WC2L;STDC2WC2X];
        GPSSPRMCB=[GPSSPRMCB SPRSAT];
        clear SPRSAT;
    end
end
%% ---------------------------subfunction-------------------------
function [output,STDoutput] = RemoveOutlier(SPR,flag)
    TempSPR = SPR;
    TempSPR(TempSPR(:, flag)==0,:)=[];
    a = 2.0;
    STD = std(TempSPR(:,flag)) * a;
    [row,~] = size(TempSPR);
    % quality control
    if row < 120 || std(TempSPR(:,flag)) > 0.6
        output    = 0;
        STDoutput = 0;
    else
        outputArg2 = mean(TempSPR(:,flag));
        SPR(SPR(:,flag)==0,:)=[];
        SPR(SPR(:,flag)>outputArg2+STD,:)=[];
        SPR(SPR(:,flag)<outputArg2-STD,:)=[];  
        output    = mean(SPR(:,flag));
        STDoutput = std(SPR(:,flag));
    end         
end
% SQR
function [outputArg1] = SQR(inputArg1)
    outputArg1 = inputArg1 * inputArg1;
end