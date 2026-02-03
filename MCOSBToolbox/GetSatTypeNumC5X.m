% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [GSatTypeNum,GsatC1W,GsatC2W,GsatC5X,GsatC1C,GsatC2X,GPRNC1C,GPRNC1W,GPRNC2W,GPRNC2X,GPRNC5X] = GetSatTypeNumC5X(NumRec,Path,list_gps)
    GPRNC1C=linspace(0,0,32);
    GPRNC1W=linspace(0,0,32);
    GPRNC2W=linspace(0,0,32);
    GPRNC2X=linspace(0,0,32);
    GPRNC5X=linspace(0,0,32);
    for i=1:NumRec
        load([Path '/' list_gps(i).name],'-mat');
        for j=1:32
            if GPSSPRMCB(1,j)~=0
                GPRNC1W(j)=GPRNC1W(j)+1;
                GPRNC2W(j)=GPRNC2W(j)+1;
            end
            if GPSSPRMCB(3,j)~=0
                GPRNC5X(j)=GPRNC5X(j)+1;
            end
            if GPSSPRMCB(4,j)~=0
                GPRNC1C(j)=GPRNC1C(j)+1;
            end 
            if GPSSPRMCB(6,j)~=0
                GPRNC2X(j)=GPRNC2X(j)+1;
            end          
        end
        clear GPSSPRMCB;
    end
    GsatC1W=find(GPRNC1W==0);
    GsatC2W=find(GPRNC2W==0);
    GsatC5X=find(GPRNC5X==0);
    GsatC1C=find(GPRNC1C==0);
    GsatC2X=find(GPRNC2X==0);
    if isempty(GsatC1W), NumSatC1W=32; else, NumSatC1W=32-length(GsatC1W); end  %the number of satellites
    if isempty(GsatC2W), NumSatC2W=32; else, NumSatC2W=32-length(GsatC2W); end  %the number of satellites
    if isempty(GsatC5X), NumSatC5X=32; else, NumSatC5X=32-length(GsatC5X); end  %the number of satellites
    if isempty(GsatC1C), NumSatC1C=32; else, NumSatC1C=32-length(GsatC1C); end  %the number of satellites
    if isempty(GsatC2X), NumSatC2X=32; else, NumSatC2X=32-length(GsatC2X); end  %the number of satellites
    GSatTypeNum = [NumSatC1C;NumSatC1W;NumSatC2W;NumSatC2X;NumSatC5X];
end