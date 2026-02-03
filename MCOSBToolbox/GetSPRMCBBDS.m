% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [] = GetSPRMCBBDS(OBSPath,MCBPath,sate,STAInformation,lim,SATMark,ION,B1C,B1I,B2A,B3I,B2B,year,current_doy)
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
        fields=fieldnames(obs);
        resultString = strjoin(fields);
        if contains(resultString,B1C)&&contains(resultString,B1I)&&contains(resultString,B2A)&&contains(resultString,B3I)&&contains(resultString,B2B) 
            [BDSSPRMCB]=BDS_SPRMCB(sate,sx,sy,sz,obs,B1C,B1I,B2A,B3I,B2B,lim,ION,sow);
            % Quality control, when nsat smaller than 5
            no_zero = ~all(BDSSPRMCB == 0, 1);          
            [~,col] = size(BDSSPRMCB(:, no_zero));
            if col < 4
                continue;
            end
            clear no_zero;
            % Different OBS types for BDS3
            if exist([MCBPath '/' doy '/' B1C '/'],'dir')==0
                mkdir([MCBPath '/' doy '/'  B1C '/']);
            end
            filenameSPR = [MCBPath '/' doy '/' B1C '/' site doy 'MCB.mat'];
            save(filenameSPR,'BDSSPRMCB','-mat');
        end
        disp(['-------- [',num2str(DOY),'] [ ',num2str(i),'/',num2str(len),' ] BDS3SPRMCB is processed !']);
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
% cut BDS obs
if ~isnan(find(strcmp(fields, 'BDSC2I' )))
    bdsline=size(obs.BDSC2I,1);
    bdssate=size(sate_mark.bds,2);
    if bdsline<2880
        obs.BDSC2I(bdsline+1:2880,:)=0;
    end
    bdsl=size(obs.BDSC2I,2);
    if bdsl<bdssate
        obs.BDSC2I(:,bdsl+1:bdssate)=0;
    end    
    bdsdelete=find(sate_mark.bds==0);
    bdsnum=size(sate.bdsx,2);    
    % if bdsnum<size(obs.BDSC2I,2)
    %     if ~isnan(bdsdelete)
    %         k=length(bdsdelete);
    %         for i=k:-1:1
    %             obs.BDSC2I(:,bdsdelete(i))=[];
    %         end
    %     end
    %     obs.BDSC2I=obs.BDSC2I(:,1:bdsnum);
    % end
    obs.BDSC2I(isnan(obs.BDSC2I))=0;    
    bdsx=sate.bdsx;bdsy=sate.bdsy;bdsz=sate.bdsz;
    for i=1:bdsnum
        for j=1:2880
            if obs.BDSC2I(j,i)==0
                obs.BDSC2I(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,bdsx(j,i),bdsy(j,i),bdsz(j,i));
            if el<lim
                obs.BDSC2I(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'BDSC1P' )))
    bdsline=size(obs.BDSC1P,1);
    bdssate=size(sate_mark.bds,2);
    if bdsline<2880
        obs.BDSC1P(bdsline+1:2880,:)=0;
    end
    bdsl=size(obs.BDSC1P,2);
    if bdsl<bdssate
        obs.BDSC1P(:,bdsl+1:bdssate)=0;
    end    
    bdsdelete=find(sate_mark.bds==0);
    bdsnum=size(sate.bdsx,2);    
    % if bdsnum<size(obs.BDSC1P,2)
    %     if ~isnan(bdsdelete)
    %         k=length(bdsdelete);
    %         for i=k:-1:1
    %             obs.BDSC1P(:,bdsdelete(i))=[];
    %         end
    %     end
    %     obs.BDSC1P=obs.BDSC1P(:,1:bdsnum);
    % end
    obs.BDSC1P(isnan(obs.BDSC1P))=0;    
    bdsx=sate.bdsx;bdsy=sate.bdsy;bdsz=sate.bdsz;
    for i=1:bdsnum
        for j=1:2880
            if obs.BDSC1P(j,i)==0
                obs.BDSC1P(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,bdsx(j,i),bdsy(j,i),bdsz(j,i));
            if el<lim
                obs.BDSC1P(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'BDSC1X' )))
    bdsline=size(obs.BDSC1X,1);
    bdssate=size(sate_mark.bds,2);
    if bdsline<2880
        obs.BDSC1X(bdsline+1:2880,:)=0;
    end
    bdsl=size(obs.BDSC1X,2);
    if bdsl<bdssate
        obs.BDSC1X(:,bdsl+1:bdssate)=0;
    end    
    bdsdelete=find(sate_mark.bds==0);
    bdsnum=size(sate.bdsx,2);    
    % if bdsnum<size(obs.BDSC1X,2)
    %     if ~isnan(bdsdelete)
    %         k=length(bdsdelete);
    %         for i=k:-1:1
    %             obs.BDSC1X(:,bdsdelete(i))=[];
    %         end
    %     end
    %     obs.BDSC1X=obs.BDSC1X(:,1:bdsnum);
    % end
    obs.BDSC1X(isnan(obs.BDSC1X))=0;    
    bdsx=sate.bdsx;bdsy=sate.bdsy;bdsz=sate.bdsz;
    for i=1:bdsnum
        for j=1:2880
            if obs.BDSC1X(j,i)==0
                obs.BDSC1X(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,bdsx(j,i),bdsy(j,i),bdsz(j,i));
            if el<lim
                obs.BDSC1X(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'BDSC5P' )))
    bdsline=size(obs.BDSC5P,1);
    bdssate=size(sate_mark.bds,2);
    if bdsline<2880
        obs.BDSC5P(bdsline+1:2880,:)=0;
    end
    bdsl=size(obs.BDSC5P,2);
    if bdsl<bdssate
        obs.BDSC5P(:,bdsl+1:bdssate)=0;
    end    
    bdsdelete=find(sate_mark.bds==0);
    bdsnum=size(sate.bdsx,2);    
    % if bdsnum<size(obs.BDSC5P,2)
    %     if ~isnan(bdsdelete)
    %         k=length(bdsdelete);
    %         for i=k:-1:1
    %             obs.BDSC5P(:,bdsdelete(i))=[];
    %         end
    %     end
    %     obs.BDSC5P=obs.BDSC5P(:,1:bdsnum);
    % end
    obs.BDSC5P(isnan(obs.BDSC5P))=0;    
    bdsx=sate.bdsx;bdsy=sate.bdsy;bdsz=sate.bdsz;
    for i=1:bdsnum
        for j=1:2880
            if obs.BDSC5P(j,i)==0
                obs.BDSC5P(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,bdsx(j,i),bdsy(j,i),bdsz(j,i));
            if el<lim
                obs.BDSC5P(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'BDSC5X' )))
    bdsline=size(obs.BDSC5X,1);
    bdssate=size(sate_mark.bds,2);
    if bdsline<2880
        obs.BDSC5X(bdsline+1:2880,:)=0;
    end
    bdsl=size(obs.BDSC5X,2);
    if bdsl<bdssate
        obs.BDSC5X(:,bdsl+1:bdssate)=0;
    end    
    bdsdelete=find(sate_mark.bds==0);
    bdsnum=size(sate.bdsx,2);    
    % if bdsnum<size(obs.BDSC5X,2)
    %     if ~isnan(bdsdelete)
    %         k=length(bdsdelete);
    %         for i=k:-1:1
    %             obs.BDSC5X(:,bdsdelete(i))=[];
    %         end
    %     end
    %     obs.BDSC5X=obs.BDSC5X(:,1:bdsnum);
    % end
    obs.BDSC5X(isnan(obs.BDSC5X))=0;    
    bdsx=sate.bdsx;bdsy=sate.bdsy;bdsz=sate.bdsz;
    for i=1:bdsnum
        for j=1:2880
            if obs.BDSC5X(j,i)==0
                obs.BDSC5X(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,bdsx(j,i),bdsy(j,i),bdsz(j,i));
            if el<lim
                obs.BDSC5X(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'BDSC6I' )))
    bdsline=size(obs.BDSC6I,1);
    bdssate=size(sate_mark.bds,2);
    if bdsline<2880
        obs.BDSC6I(bdsline+1:2880,:)=0;
    end
    bdsl=size(obs.BDSC6I,2);
    if bdsl<bdssate
        obs.BDSC6I(:,bdsl+1:bdssate)=0;
    end    
    bdsdelete=find(sate_mark.bds==0);
    bdsnum=size(sate.bdsx,2);    
    % if bdsnum<size(obs.BDSC6I,2)
    %     if ~isnan(bdsdelete)
    %         k=length(bdsdelete);
    %         for i=k:-1:1
    %             obs.BDSC6I(:,bdsdelete(i))=[];
    %         end
    %     end
    %     obs.BDSC6I=obs.BDSC6I(:,1:bdsnum);
    % end
    obs.BDSC6I(isnan(obs.BDSC6I))=0;    
    bdsx=sate.bdsx;bdsy=sate.bdsy;bdsz=sate.bdsz;
    for i=1:bdsnum
        for j=1:2880
            if obs.BDSC6I(j,i)==0
                obs.BDSC6I(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,bdsx(j,i),bdsy(j,i),bdsz(j,i));
            if el<lim
                obs.BDSC6I(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'BDSC7D' )))
    bdsline=size(obs.BDSC7D,1);
    bdssate=size(sate_mark.bds,2);
    if bdsline<2880
        obs.BDSC7D(bdsline+1:2880,:)=0;
    end
    bdsl=size(obs.BDSC7D,2);
    if bdsl<bdssate
        obs.BDSC7D(:,bdsl+1:bdssate)=0;
    end    
    bdsdelete=find(sate_mark.bds==0);
    bdsnum=size(sate.bdsx,2);    
    % if bdsnum<size(obs.BDSC7D,2)
    %     if ~isnan(bdsdelete)
    %         k=length(bdsdelete);
    %         for i=k:-1:1
    %             obs.BDSC7D(:,bdsdelete(i))=[];
    %         end
    %     end
    %     obs.BDSC7D=obs.BDSC7D(:,1:bdsnum);
    % end
    obs.BDSC7D(isnan(obs.BDSC7D))=0;    
    bdsx=sate.bdsx;bdsy=sate.bdsy;bdsz=sate.bdsz;
    for i=1:bdsnum
        for j=1:2880
            if obs.BDSC7D(j,i)==0
                obs.BDSC7D(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,bdsx(j,i),bdsy(j,i),bdsz(j,i));
            if el<lim
                obs.BDSC7D(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'BDSC7Z' )))
    bdsline=size(obs.BDSC7Z,1);
    bdssate=size(sate_mark.bds,2);
    if bdsline<2880
        obs.BDSC7Z(bdsline+1:2880,:)=0;
    end
    bdsl=size(obs.BDSC7Z,2);
    if bdsl<bdssate
        obs.BDSC7Z(:,bdsl+1:bdssate)=0;
    end    
    bdsdelete=find(sate_mark.bds==0);
    bdsnum=size(sate.bdsx,2);    
    % if bdsnum<size(obs.BDSC7Z,2)
    %     if ~isnan(bdsdelete)
    %         k=length(bdsdelete);
    %         for i=k:-1:1
    %             obs.BDSC7Z(:,bdsdelete(i))=[];
    %         end
    %     end
    %     obs.BDSC7Z=obs.BDSC7Z(:,1:bdsnum);
    % end
    obs.BDSC7Z(isnan(obs.BDSC7Z))=0;    
    bdsx=sate.bdsx;bdsy=sate.bdsy;bdsz=sate.bdsz;
    for i=1:bdsnum
        for j=1:2880
            if obs.BDSC7Z(j,i)==0
                obs.BDSC7Z(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,bdsx(j,i),bdsy(j,i),bdsz(j,i));
            if el<lim
                obs.BDSC7Z(j,i)=0;
                continue;
            end
        end
    end
end
end
%% ---------------------------subfunction-------------------------
function [ALLBDSSPRMCB]=BDS_SPRMCB(sate,sx,sy,sz,obs,B1C,B1I,B2A,B3I,B2B,lim,IONO,sow)
    size2=46;
    BDSSPR123 = zeros(2880,46);
    BDSSPR234 = zeros(2880,46);
    BDSSPR235 = zeros(2880,46);
    BDSSPR23  = zeros(2880,46);
    % Constants
    CLIGHT = 299792458;   
    FB1C = 1.57542E9;
    FB1I = 1.561098E9;
    FB3I = 1.26852E9;
    FB2B = 1.20714E9;
    FB2A = 1.17645E9;
    % 1
    WLB1C = CLIGHT/FB1C; 
    % 2
    WLB1I = CLIGHT/FB1I;
    % 3
    WLB3I = CLIGHT/FB3I; 
    % 4
    WLB2B = CLIGHT/FB2B;
    % 5
    WLB2A = CLIGHT/FB2A;
    A32 = SQR(WLB3I) - SQR(WLB1I);
    A13 = SQR(WLB1C) - SQR(WLB3I);
    A21 = SQR(WLB1I) - SQR(WLB1C);
    A43 = SQR(WLB2B) - SQR(WLB3I);
    A53 = SQR(WLB2A) - SQR(WLB3I);
    A24 = SQR(WLB1I) - SQR(WLB2B);
    A25 = SQR(WLB1I) - SQR(WLB2A);  
    [B,L,H]=XYZtoBLH(sx,sy,sz);
    bdsx=sate.bdsx;bdsy=sate.bdsy;bdsz=sate.bdsz;
    %i is PRN number
    for i=1:size2 
        % Exclude the BDS2 satellites
        if i < 19
            continue;
        end
        for j=1:2880
            OBSB1I = obs.(B1I)(j, i);
            if OBSB1I == 0
                continue;
            end     
            [E,A]=Get_EA(sx,sy,sz,bdsx(j,i),bdsy(j,i),bdsz(j,i));
            if E<lim
                continue;
            end
            % Calculate the ionospheric delay
            OBSB1I = obs.(B1I)(j, i);
            OBSB3I = obs.(B3I)(j, i);
            CodeDistance = [OBSB1I;OBSB3I];
            CodeDistance(CodeDistance(:,1)==0,:)=[];
            tau = mean(CodeDistance)/CLIGHT;
            Ttr = sow + (30*(j-1)) - tau;
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
            % get value of mapping-function and IPP
            [mappingf, Lat_IPP, Lon_IPP] = iono_mf(elev, mf, pos_WGS84, az, radius, hgt);
            vtec = iono_gims(Lat_IPP, Lon_IPP, Ttr,IONO.ionex,interpol); % interpolate VTEC
            % GF combination B1I/B3I(C2I/C6I)
            iono(1) = mappingf * 40.3e16/FB1I^2* vtec;
            iono(2) = mappingf * 40.3e16/FB3I^2* vtec;
            % Calculate SPRMCB
            OBSB1C = obs.(B1C)(j, i);
            OBSB1I = obs.(B1I)(j, i);
            OBSB3I = obs.(B3I)(j, i);
            OBSB2B = obs.(B2B)(j, i);
            OBSB2A = obs.(B2A)(j, i);  
            % GFIF
            BDSSPR123(j,i) = A32 * OBSB1C + A13 * OBSB1I + A21 * OBSB3I + 0.0 * OBSB2B + 0.0 * OBSB2A;
            BDSSPR234(j,i) = 0.0 * OBSB1C + A43 * OBSB1I + A24 * OBSB3I + A32 * OBSB2B + 0.0 * OBSB2A;
            BDSSPR235(j,i) = 0.0 * OBSB1C + A53 * OBSB1I + A25 * OBSB3I + 0.0 * OBSB2B + A32 * OBSB2A;
            % GF
            BDSSPR23(j,i) = (OBSB1I - OBSB3I) - (iono(1) - iono(2));
            % Quality control
            if OBSB1C == 0 || OBSB1I == 0 || OBSB3I==0 || OBSB2B==0 || OBSB2A==0
                BDSSPR123(j,i)=0.0; BDSSPR234(j,i)=0.0;BDSSPR235(j,i)=0.0;
            end
            if OBSB1I == 0 || OBSB3I == 0 || vtec == 0
                BDSSPR23(j,i)=0.0;
            end      
        end
    end
    % Quality control
    for ii = 1:size2
        SATBDSSPRMCB=[];
        SATSTDBDSSPRMCB=[];
        % Exclude the BDS2 satellites
        if ii < 19
            continue;
        end
        [BDSSPRMCB,STDBDSSPRMCB] = RemoveOutlier(BDSSPR123,ii);
        SATBDSSPRMCB=[SATBDSSPRMCB;BDSSPRMCB];
        SATSTDBDSSPRMCB=[SATSTDBDSSPRMCB;STDBDSSPRMCB];
        [BDSSPRMCB,STDBDSSPRMCB] = RemoveOutlier(BDSSPR234,ii);
        SATBDSSPRMCB=[SATBDSSPRMCB;BDSSPRMCB];
        SATSTDBDSSPRMCB=[SATSTDBDSSPRMCB;STDBDSSPRMCB];
        [BDSSPRMCB,STDBDSSPRMCB] = RemoveOutlier(BDSSPR235,ii);
        SATBDSSPRMCB=[SATBDSSPRMCB;BDSSPRMCB];
        SATSTDBDSSPRMCB=[SATSTDBDSSPRMCB;STDBDSSPRMCB];

        [outBDSSPR23,STDBDSSPR23] = RemoveOutlierGF(BDSSPR23,ii); 
        % No reference GF value
        if outBDSSPR23 == 0
            SATBDSSPRMCB(:,1)    = 0;
            SATSTDBDSSPRMCB(:,1) = 0;
            STDBDSSPR23          = 0;
        end
        ALLBDSSPRMCB(:,ii-18) = [SATBDSSPRMCB;outBDSSPR23;SATSTDBDSSPRMCB;STDBDSSPR23];
        clear SATBDSSPRMCB SATSTDBDSSPRMCB;
    end
end
%% ---------------------------subfunction-------------------------
function [output,STDoutput] = RemoveOutlier(SPR,flag)
    TempSPR = SPR;
    TempSPR(TempSPR(:, flag)==0,:)=[];
    a = 2.0;
    STD = std(TempSPR(:,flag)) * a;

    [row,~] = size(TempSPR);
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
%% ---------------------------subfunction-------------------------
function [output,STDout] = RemoveOutlierGF(SPR12,flag)
    TempSPR = SPR12;
    TempSPR(TempSPR(:, flag)==0,:)=[];
    a = 2.0;
    STD = std(TempSPR(:,flag))*a;

    [row,~] = size(TempSPR);
    if row < 120 || std(TempSPR(:,flag)) > 0.6
        output = 0;
        STDout=0;
    else
        outputArg2 = mean(TempSPR(:,flag));
        SPR12(SPR12(:,flag)==0,:)=[];
        SPR12(SPR12(:,flag)>outputArg2+STD,:)=[];
        SPR12(SPR12(:,flag)<outputArg2-STD,:)=[];
        output  = mean(SPR12(:,flag));
        STDout = std(SPR12(:,flag));
    end         
end
% SQR
function [outputArg1] = SQR(inputArg1)
    outputArg1 = inputArg1 * inputArg1;
end