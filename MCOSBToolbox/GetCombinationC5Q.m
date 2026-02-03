% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [C1WC2W,C1WC2WC5Q,C1WC1C,C2WC2L] = GetCombinationC5Q(GPSSPRMCB,GsatC1W,GsatC5Q,GsatC1C,GsatC2L)
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

    C1WC2WC5Q = [];
    if ~isempty(GsatC5Q)
        for k = 1:32
            if sum(GsatC5Q==k) > 0
               continue;
            else
               C1WC2WC5Q = [C1WC2WC5Q GPSSPRMCB(:,k)];
            end
        end
    else
        C1WC2WC5Q = GPSSPRMCB;
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

    C2WC2L = [];
    if ~isempty(GsatC2L)
        for k = 1:32
            if sum(GsatC2L==k) > 0
               continue;
            else
               C2WC2L = [C2WC2L GPSSPRMCB(:,k)];
            end
        end
    else      
        C2WC2L = GPSSPRMCB;
    end
end