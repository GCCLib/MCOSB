% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [C1WC2W,C1WC2WC5X,C1WC1C,C2WC2X] = GetCombinationC5X(GPSSPRMCB,GsatC1W,GsatC5X,GsatC1C,GsatC2X)
% Delete corresponding colum without observed satellites
    C1WC2W = [];
    if ~isempty(GsatC1W)
        for k = 1:32
            if sum(GsatC1W==k) > 0
               continue;
            else
               C1WC2W = [C1WC2W GPSSPRMCB(:,k)];
            end
        end
    else    
        %Full
        C1WC2W = GPSSPRMCB;
    end

    C1WC2WC5X = [];
    if ~isempty(GsatC5X)
        for k = 1:32
            if sum(GsatC5X==k) > 0
               continue;
            else
               C1WC2WC5X = [C1WC2WC5X GPSSPRMCB(:,k)];
            end
        end
    else      
        C1WC2WC5X = GPSSPRMCB;
    end

    C1WC1C = [];
    if ~isempty(GsatC1C)
        for k = 1:32
            if sum(GsatC1C==k) > 0
               continue;
            else
               C1WC1C = [C1WC1C GPSSPRMCB(:,k)];
            end
        end
    else     
        C1WC1C = GPSSPRMCB;
    end
    
    C2WC2X = [];
    if ~isempty(GsatC2X)
        for k = 1:32
            if sum(GsatC2X==k) > 0
               continue;
            else
               C2WC2X = [C2WC2X GPSSPRMCB(:,k)];
            end
        end
    else     
        C2WC2X = GPSSPRMCB;
    end
end