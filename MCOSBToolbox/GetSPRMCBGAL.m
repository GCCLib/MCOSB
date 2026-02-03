% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [] = GetSPRMCBGAL(OBSPath,MCBPath,sate,STAInformation,lim,SATMark,ION,E1,E6,E5B,E5,E5A,year,current_doy)
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
        if contains(resultString,E1)&&contains(resultString,E6)&&contains(resultString,E5B)&&contains(resultString,E5)&&contains(resultString,E5A) 
            [GALSPRMCB]=GAL_SPRMCB(sate,sx,sy,sz,obs,E1,E6,E5B,E5,E5A,lim,ION,sow);
            % Quality control, when nsat smaller than 5
            no_zero = ~all(GALSPRMCB == 0, 1);          
            [~,col] = size(GALSPRMCB(:, no_zero));
            if col < 4
                continue;
            end
            clear no_zero;
            % Different OBS types for GAL
            if exist([MCBPath '/' doy '/' E1 '/'],'dir')==0
                mkdir([MCBPath '/' doy '/'  E1 '/']);
            end
            filenameSPR = [MCBPath '/' doy '/' E1 '/' site doy 'MCB.mat'];
            save(filenameSPR,'GALSPRMCB','-mat');
        end
        disp(['-------- [',num2str(DOY),'] [ ',num2str(i),'/',num2str(len),' ] GALSPRMCB is processed !']);
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
% cut GAL obs
if ~isnan(find(strcmp(fields, 'GALC1C' )))
    galline=size(obs.GALC1C,1);
    galsate=size(sate_mark.gal,2);
    if galline<2880
        obs.GALC1C(galline+1:2880,:)=0;
    end
    gall=size(obs.GALC1C,2);
    if gall<galsate
        obs.GALC1C(:,gall+1:galsate)=0;
    end    
    galdelete=find(sate_mark.gal==0);
    galnum=size(sate.galx,2);    
    % if galnum<size(obs.GALC1C,2)
    %     if ~isnan(galdelete)
    %         k=length(galdelete);
    %         for i=k:-1:1
    %             obs.GALC1C(:,galdelete(i))=[];
    %         end
    %     end
    %     obs.GALC1C=obs.GALC1C(:,1:galnum);
    % end
    obs.GALC1C(isnan(obs.GALC1C))=0;    
    galx=sate.galx;galy=sate.galy;galz=sate.galz;
    for i=1:galnum
        for j=1:2880
            if obs.GALC1C(j,i)==0
                obs.GALC1C(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,galx(j,i),galy(j,i),galz(j,i));
            if el<lim
                obs.GALC1C(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'GALC1X' )))
    galline=size(obs.GALC1X,1);
    galsate=size(sate_mark.gal,2);
    if galline<2880
        obs.GALC1X(galline+1:2880,:)=0;
    end
    gall=size(obs.GALC1X,2);
    if gall<galsate
        obs.GALC1X(:,gall+1:galsate)=0;
    end    
    galdelete=find(sate_mark.gal==0);
    galnum=size(sate.galx,2);    
    % if galnum<size(obs.GALC1X,2)
    %     if ~isnan(galdelete)
    %         k=length(galdelete);
    %         for i=k:-1:1
    %             obs.GALC1X(:,galdelete(i))=[];
    %         end
    %     end
    %     obs.GALC1X=obs.GALC1X(:,1:galnum);
    % end
    obs.GALC1X(isnan(obs.GALC1X))=0;    
    galx=sate.galx;galy=sate.galy;galz=sate.galz;
    for i=1:galnum
        for j=1:2880
            if obs.GALC1X(j,i)==0
                obs.GALC1X(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,galx(j,i),galy(j,i),galz(j,i));
            if el<lim
                obs.GALC1X(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'GALC6C' )))
    galline=size(obs.GALC6C,1);
    galsate=size(sate_mark.gal,2);
    if galline<2880
        obs.GALC6C(galline+1:2880,:)=0;
    end
    gall=size(obs.GALC6C,2);
    if gall<galsate
        obs.GALC6C(:,gall+1:galsate)=0;
    end    
    galdelete=find(sate_mark.gal==0);
    galnum=size(sate.galx,2);    
    % if galnum<size(obs.GALC6C,2)
    %     if ~isnan(galdelete)
    %         k=length(galdelete);
    %         for i=k:-1:1
    %             obs.GALC6C(:,galdelete(i))=[];
    %         end
    %     end
    %     obs.GALC6C=obs.GALC6C(:,1:galnum);
    % end
    obs.GALC6C(isnan(obs.GALC6C))=0;    
    galx=sate.galx;galy=sate.galy;galz=sate.galz;
    for i=1:galnum
        for j=1:2880
            if obs.GALC6C(j,i)==0
                obs.GALC6C(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,galx(j,i),galy(j,i),galz(j,i));
            if el<lim
                obs.GALC6C(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'GALC6X' )))
    galline=size(obs.GALC6X,1);
    galsate=size(sate_mark.gal,2);
    if galline<2880
        obs.GALC6X(galline+1:2880,:)=0;
    end
    gall=size(obs.GALC6X,2);
    if gall<galsate
        obs.GALC6X(:,gall+1:galsate)=0;
    end    
    galdelete=find(sate_mark.gal==0);
    galnum=size(sate.galx,2);    
    % if galnum<size(obs.GALC6X,2)
    %     if ~isnan(galdelete)
    %         k=length(galdelete);
    %         for i=k:-1:1
    %             obs.GALC6X(:,galdelete(i))=[];
    %         end
    %     end
    %     obs.GALC6X=obs.GALC6X(:,1:galnum);
    % end
    obs.GALC6X(isnan(obs.GALC6X))=0;    
    galx=sate.galx;galy=sate.galy;galz=sate.galz;
    for i=1:galnum
        for j=1:2880
            if obs.GALC6X(j,i)==0
                obs.GALC6X(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,galx(j,i),galy(j,i),galz(j,i));
            if el<lim
                obs.GALC6X(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'GALC7Q' )))
    galline=size(obs.GALC7Q,1);
    galsate=size(sate_mark.gal,2);
    if galline<2880
        obs.GALC7Q(galline+1:2880,:)=0;
    end
    gall=size(obs.GALC7Q,2);
    if gall<galsate
        obs.GALC7Q(:,gall+1:galsate)=0;
    end    
    galdelete=find(sate_mark.gal==0);
    galnum=size(sate.galx,2);    
    % if galnum<size(obs.GALC7Q,2)
    %     if ~isnan(galdelete)
    %         k=length(galdelete);
    %         for i=k:-1:1
    %             obs.GALC7Q(:,galdelete(i))=[];
    %         end
    %     end
    %     obs.GALC7Q=obs.GALC7Q(:,1:galnum);
    % end
    obs.GALC7Q(isnan(obs.GALC7Q))=0;    
    galx=sate.galx;galy=sate.galy;galz=sate.galz;
    for i=1:galnum
        for j=1:2880
            if obs.GALC7Q(j,i)==0
                obs.GALC7Q(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,galx(j,i),galy(j,i),galz(j,i));
            if el<lim
                obs.GALC7Q(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'GALC7X' )))
    galline=size(obs.GALC7X,1);
    galsate=size(sate_mark.gal,2);
    if galline<2880
        obs.GALC7X(galline+1:2880,:)=0;
    end
    gall=size(obs.GALC7X,2);
    if gall<galsate
        obs.GALC7X(:,gall+1:galsate)=0;
    end    
    galdelete=find(sate_mark.gal==0);
    galnum=size(sate.galx,2);    
    % if galnum<size(obs.GALC7X,2)
    %     if ~isnan(galdelete)
    %         k=length(galdelete);
    %         for i=k:-1:1
    %             obs.GALC7X(:,galdelete(i))=[];
    %         end
    %     end
    %     obs.GALC7X=obs.GALC7X(:,1:galnum);
    % end
    obs.GALC7X(isnan(obs.GALC7X))=0;    
    galx=sate.galx;galy=sate.galy;galz=sate.galz;
    for i=1:galnum
        for j=1:2880
            if obs.GALC7X(j,i)==0
                obs.GALC7X(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,galx(j,i),galy(j,i),galz(j,i));
            if el<lim
                obs.GALC7X(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'GALC8Q' )))
    galline=size(obs.GALC8Q,1);
    galsate=size(sate_mark.gal,2);
    if galline<2880
        obs.GALC8Q(galline+1:2880,:)=0;
    end
    gall=size(obs.GALC8Q,2);
    if gall<galsate
        obs.GALC8Q(:,gall+1:galsate)=0;
    end    
    galdelete=find(sate_mark.gal==0);
    galnum=size(sate.galx,2);    
    % if galnum<size(obs.GALC8Q,2)
    %     if ~isnan(galdelete)
    %         k=length(galdelete);
    %         for i=k:-1:1
    %             obs.GALC8Q(:,galdelete(i))=[];
    %         end
    %     end
    %     obs.GALC8Q=obs.GALC8Q(:,1:galnum);
    % end
    obs.GALC8Q(isnan(obs.GALC8Q))=0;    
    galx=sate.galx;galy=sate.galy;galz=sate.galz;
    for i=1:galnum
        for j=1:2880
            if obs.GALC8Q(j,i)==0
                obs.GALC8Q(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,galx(j,i),galy(j,i),galz(j,i));
            if el<lim
                obs.GALC8Q(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'GALC8X' )))
    galline=size(obs.GALC8X,1);
    galsate=size(sate_mark.gal,2);
    if galline<2880
        obs.GALC8X(galline+1:2880,:)=0;
    end
    gall=size(obs.GALC8X,2);
    if gall<galsate
        obs.GALC8X(:,gall+1:galsate)=0;
    end    
    galdelete=find(sate_mark.gal==0);
    galnum=size(sate.galx,2);    
    % if galnum<size(obs.GALC8X,2)
    %     if ~isnan(galdelete)
    %         k=length(galdelete);
    %         for i=k:-1:1
    %             obs.GALC8X(:,galdelete(i))=[];
    %         end
    %     end
    %     obs.GALC8X=obs.GALC8X(:,1:galnum);
    % end
    obs.GALC8X(isnan(obs.GALC8X))=0;    
    galx=sate.galx;galy=sate.galy;galz=sate.galz;
    for i=1:galnum
        for j=1:2880
            if obs.GALC8X(j,i)==0
                obs.GALC8X(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,galx(j,i),galy(j,i),galz(j,i));
            if el<lim
                obs.GALC8X(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'GALC5Q' )))
    galline=size(obs.GALC5Q,1);
    galsate=size(sate_mark.gal,2);
    if galline<2880
        obs.GALC5Q(galline+1:2880,:)=0;
    end
    gall=size(obs.GALC5Q,2);
    if gall<galsate
        obs.GALC5Q(:,gall+1:galsate)=0;
    end    
    galdelete=find(sate_mark.gal==0);
    galnum=size(sate.galx,2);    
    % if galnum<size(obs.GALC5Q,2)
    %     if ~isnan(galdelete)
    %         k=length(galdelete);
    %         for i=k:-1:1
    %             obs.GALC5Q(:,galdelete(i))=[];
    %         end
    %     end
    %     obs.GALC5Q=obs.GALC5Q(:,1:galnum);
    % end
    obs.GALC5Q(isnan(obs.GALC5Q))=0;    
    galx=sate.galx;galy=sate.galy;galz=sate.galz;
    for i=1:galnum
        for j=1:2880
            if obs.GALC5Q(j,i)==0
                obs.GALC5Q(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,galx(j,i),galy(j,i),galz(j,i));
            if el<lim
                obs.GALC5Q(j,i)=0;
                continue;
            end
        end
    end
end
if ~isnan(find(strcmp(fields, 'GALC5X' )))
    galline=size(obs.GALC5X,1);
    galsate=size(sate_mark.gal,2);
    if galline<2880
        obs.GALC5X(galline+1:2880,:)=0;
    end
    gall=size(obs.GALC5X,2);
    if gall<galsate
        obs.GALC5X(:,gall+1:galsate)=0;
    end    
    galdelete=find(sate_mark.gal==0);
    galnum=size(sate.galx,2);    
    % if galnum<size(obs.GALC5X,2)
    %     if ~isnan(galdelete)
    %         k=length(galdelete);
    %         for i=k:-1:1
    %             obs.GALC5X(:,galdelete(i))=[];
    %         end
    %     end
    %     obs.GALC5X=obs.GALC5X(:,1:galnum);
    % end
    obs.GALC5X(isnan(obs.GALC5X))=0;    
    galx=sate.galx;galy=sate.galy;galz=sate.galz;
    for i=1:galnum
        for j=1:2880
            if obs.GALC5X(j,i)==0
                obs.GALC5X(j,i)=0;
                continue;
            end
            [el,~]=Get_EA(sx,sy,sz,galx(j,i),galy(j,i),galz(j,i));
            if el<lim
                obs.GALC5X(j,i)=0;
                continue;
            end
        end
    end
end
end
%% ---------------------------subfunction-------------------------
function [ALLGALSPRMCB]=GAL_SPRMCB(sate,sx,sy,sz,obs,E1,E6,E5B,E5,E5A,lim,IONO,sow)
    size2=36;
    GALSPR125 = zeros(2880,36);
    GALSPR135 = zeros(2880,36);
    GALSPR145 = zeros(2880,36);
    GALSPR15  = zeros(2880,36);
    % Constants
    CLIGHT = 299792458;   
    FE1 = 1.57542E9;
    FE6 = 1.27875E9;
    FE5B = 1.20714E9;
    FE5 = 1.191795E9;
    FE5A = 1.17645E9;
    % 1
    WLE1 = CLIGHT/FE1; 
    % 2
    WLE6 = CLIGHT/FE6;
    % 3
    WLE5B = CLIGHT/FE5B; 
    % 4
    WLE5 = CLIGHT/FE5;
    % 5
    WLE5A = CLIGHT/FE5A;  
    A21 = SQR(WLE6) - SQR(WLE1);
    A52 = SQR(WLE5A) - SQR(WLE6);
    A15 = SQR(WLE1) - SQR(WLE5A);
    A31 = SQR(WLE5B) - SQR(WLE1);
    A53 = SQR(WLE5A) - SQR(WLE5B);
    A54 = SQR(WLE5A) - SQR(WLE5);
    A41 = SQR(WLE5) - SQR(WLE1);   
    [B,L,H]=XYZtoBLH(sx,sy,sz);
    galx=sate.galx;galy=sate.galy;galz=sate.galz;
    %i is PRN number
    for i=1:size2 
        for j=1:2880
            OBSE1 = obs.(E1)(j, i);
            if OBSE1 == 0
                continue;
            end     
            [E,A]=Get_EA(sx,sy,sz,galx(j,i),galy(j,i),galz(j,i));
            if E<lim
                continue;
            end
            % Calculate the ionospheric delay
            OBSE1 = obs.(E1)(j, i);
            OBSE5A = obs.(E5A)(j, i);
            CodeDistance = [OBSE1;OBSE5A];
            CodeDistance(CodeDistance(:,1)==0,:)=[];
            tau = mean(CodeDistance)/CLIGHT;
            Ttr = sow + (30*(j-1)) - tau;
            % GIM
            elev = E*180/pi;
            mf = 'COSZ';
            az = A*180/pi;
            pos_WGS84.lat=B;
            pos_WGS84.lon=L;
            pos_WGS84.h=H;
            radius = 6371;
            hgt = [450 450 0];
            interpol = 'Consecutive Rotated Maps';
            % get value of mapping-function and IPP
            [mappingf, Lat_IPP, Lon_IPP] = iono_mf(elev, mf, pos_WGS84, az, radius, hgt);
            vtec = iono_gims(Lat_IPP, Lon_IPP, Ttr, IONO.ionex, interpol); % interpolate VTEC
            % GF combination E1/E5A(C1X/C5X or C1C/C5Q)
            iono(1) = mappingf * 40.3e16/FE1^2 * vtec;
            iono(2) = mappingf * 40.3e16/FE5A^2* vtec;
            % Calculate SPRMCB
            OBSE1 = obs.(E1)(j, i);
            OBSE6 = obs.(E6)(j, i);
            OBSE5B = obs.(E5B)(j, i);
            OBSE5 = obs.(E5)(j, i);
            OBSE5A = obs.(E5A)(j, i);
            % GFIF
            GALSPR125(j,i) = A52 * OBSE1 + A15 * OBSE6 + 0.0 * OBSE5B + 0.0 * OBSE5 + A21 * OBSE5A;
            GALSPR135(j,i) = A53 * OBSE1 + 0.0 * OBSE6 + A15 * OBSE5B + 0.0 * OBSE5 + A31 * OBSE5A;
            GALSPR145(j,i) = A54 * OBSE1 + 0.0 * OBSE6 + 0.0 * OBSE5B + A15 * OBSE5 + A41 * OBSE5A;
            GALSPR15(j,i) = (OBSE1 - OBSE5A) - (iono(1) - iono(2));
            % Quality control
            if OBSE1 == 0 || OBSE6 == 0 || OBSE5B==0 || OBSE5==0 || OBSE5A==0
                GALSPR125(j,i)=0.0;GALSPR135(j,i)=0.0;GALSPR145(j,i)=0.0;
            end
            if OBSE1 == 0 || OBSE5A == 0 || vtec == 0
                GALSPR15(j,i)=0.0;
            end      
        end
    end
    % Quality control
    for ii = 1:size2
        SATGALSPRMCB=[];
        SATSTDGALSPRMCB=[]; 
    
        [GALSPRMCB,STDGALSPRMCB] = RemoveOutlier(GALSPR125,ii);
        SATGALSPRMCB=[SATGALSPRMCB;GALSPRMCB];
        SATSTDGALSPRMCB=[SATSTDGALSPRMCB;STDGALSPRMCB];
        [GALSPRMCB,STDGALSPRMCB] = RemoveOutlier(GALSPR135,ii);
        SATGALSPRMCB=[SATGALSPRMCB;GALSPRMCB];
        SATSTDGALSPRMCB=[SATSTDGALSPRMCB;STDGALSPRMCB];
        [GALSPRMCB,STDGALSPRMCB] = RemoveOutlier(GALSPR145,ii);
        SATGALSPRMCB=[SATGALSPRMCB;GALSPRMCB];
        SATSTDGALSPRMCB=[SATSTDGALSPRMCB;STDGALSPRMCB];
    
        [outGALSPR15,STDGALSPR15] = RemoveOutlierGF(GALSPR15,ii); 
        % No reference GF value
        if outGALSPR15 == 0
            SATGALSPRMCB(:,1) = 0;
            SATSTDGALSPRMCB(:,1) = 0;
            STDGALSPR15=0;
        end
        ALLGALSPRMCB(:,ii) = [SATGALSPRMCB;outGALSPR15;SATSTDGALSPRMCB;STDGALSPR15];
        clear SATGALSPRMCB SATSTDGALSPRMCB;
    end
end
%% ---------------------------subfunction-------------------------
function [output,STDoutput] = RemoveOutlier(SPR,flag)
    TempSPR = SPR;
    TempSPR(TempSPR(:, flag)==0,:)=[];
    a = 2.0;
    STD = std(TempSPR(:,flag))*a;

    [row,~] = size(TempSPR);
    if row < 120 || std(TempSPR(:,flag)) > 0.6
        output = 0;
        STDoutput=0;
    else
        outputArg2 = mean(TempSPR(:,flag));
        SPR(SPR(:,flag)==0,:)=[];
        SPR(SPR(:,flag)>outputArg2+STD,:)=[];
        SPR(SPR(:,flag)<outputArg2-STD,:)=[];
        output = mean(SPR(:,flag));  
        STDoutput=std(SPR(:,flag));
    end 
end
%% ---------------------------subfunction-------------------------
function [output,STDoutput] = RemoveOutlierGF(SPR12,flag)
    TempSPR = SPR12;
    TempSPR(TempSPR(:, flag)==0,:)=[];
    a = 2.0;
    STD = std(TempSPR(:,flag))*a;

    [row,~] = size(TempSPR);
    if row < 120 || std(TempSPR(:,flag)) > 0.6
        output = 0;
        STDoutput=0;
    else
        outputArg2 = mean(TempSPR(:,flag));
        SPR12(SPR12(:,flag)==0,:)=[];
        SPR12(SPR12(:,flag)>outputArg2+STD,:)=[];
        SPR12(SPR12(:,flag)<outputArg2-STD,:)=[];
        output = mean(SPR12(:,flag));
        STDoutput=std(SPR12(:,flag));
    end         
end
% SQR
function [outputArg1] = SQR(inputArg1)
    outputArg1 = inputArg1 * inputArg1;
end