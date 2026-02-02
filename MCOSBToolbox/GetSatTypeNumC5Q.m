% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [GSatTypeNum,GsatC1W,GsatC2W,GsatC5Q,GsatC1C,GsatC2L,GPRNC1C,GPRNC1W,GPRNC2W,GPRNC2L,GPRNC5Q] = GetSatTypeNumC5Q(NumRec,Path,list_gps)
    GPRNC1C=linspace(0,0,32);
    GPRNC1W=linspace(0,0,32);
    GPRNC2W=linspace(0,0,32);
    GPRNC2L=linspace(0,0,32);
    GPRNC5Q=linspace(0,0,32);
    for i=1:NumRec
        load([Path '/' list_gps(i).name],'-mat');
        for j=1:32
            if GPSSPRMCB(1,j)~=0
                GPRNC1W(j)=GPRNC1W(j)+1;
                GPRNC2W(j)=GPRNC2W(j)+1;
            end
            if GPSSPRMCB(2,j)~=0
                GPRNC5Q(j)=GPRNC5Q(j)+1;
            end
            if GPSSPRMCB(4,j)~=0
                GPRNC1C(j)=GPRNC1C(j)+1;
            end 
            if GPSSPRMCB(5,j)~=0
                GPRNC2L(j)=GPRNC2L(j)+1;
            end          
        end
        clear GPSSPRMCB;
    end
    GsatC1W = find(GPRNC1W==0);
    GsatC2W = find(GPRNC2W==0);
    GsatC5Q = find(GPRNC5Q==0);
    GsatC1C = find(GPRNC1C==0);
    GsatC2L = find(GPRNC2L==0); 
    if isempty(GsatC1W), NumSatC1W=32; else, NumSatC1W=32-length(GsatC1W); end
    if isempty(GsatC2W), NumSatC2W=32; else, NumSatC2W=32-length(GsatC2W); end
    if isempty(GsatC5Q), NumSatC5Q=32; else, NumSatC5Q=32-length(GsatC5Q); end
    if isempty(GsatC1C), NumSatC1C=32; else, NumSatC1C=32-length(GsatC1C); end
    if isempty(GsatC2L), NumSatC2L=32; else, NumSatC2L=32-length(GsatC2L); end
    GSatTypeNum = [NumSatC1C;NumSatC1W;NumSatC2W;NumSatC2L;NumSatC5Q];
end