% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [ESatTypeNum,Esat,EPRN] = GetSatTypeNumGAL(NumRec,Path,list_gal)
    EPRN=linspace(0,0,36);
    for i=1:NumRec
        load([Path '/' list_gal(i).name],'-mat');
        for j=1:36
            if GALSPRMCB(4,j)~=0
                EPRN(j)=EPRN(j)+1;
            end
        end
        clear GALSPRMCB;
    end   
    Esat=find(EPRN==0);
    if isempty(Esat)
        ESatTypeNum=36;
    else
        % the number of satellites
        ESatTypeNum=36-length(Esat);
    end
end