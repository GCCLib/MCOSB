% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [GPSSatOSB] = GetSatOSBGPSC5X(SatC1C,SatC1W,SatC2W,SatC2X,GPRNC1C,GPRNC1W,GPRNC2W,GPRNC2X,GPRNC5X,GPSOSB)
    % L1/L2/L5
    count = 0;
    for i=1:32
        if GPRNC1C(1,i)>0
            count = count + 1;
            GPSSatOSB(i).C1C = GPSOSB(count,1);
        else
            GPSSatOSB(i).C1C = 0;
        end
    end

    count = 0;
    for i=1:32
        if GPRNC1W(1,i)>0
            count = count + 1;
            GPSSatOSB(i).C1W = GPSOSB(count+SatC1C,1);
        else
            GPSSatOSB(i).C1W = 0;
        end
    end

    count = 0;
    for i=1:32
        if GPRNC2W(1,i)>0
            count = count + 1;
            GPSSatOSB(i).C2W = GPSOSB(count+SatC1C+SatC1W,1);
        else
            GPSSatOSB(i).C2W = 0;
        end
    end

    count = 0;
    for i=1:32
        if GPRNC2X(1,i)>0
            count = count + 1;
            GPSSatOSB(i).C2X = GPSOSB(count+SatC1C+SatC1W+SatC2W,1);
        else
            GPSSatOSB(i).C2X = 0;
        end
    end

    count = 0;
    for i=1:32
        if GPRNC5X(1,i)>0
            count = count + 1;
            GPSSatOSB(i).C5X = GPSOSB(count+SatC1C+SatC1W+SatC2W+SatC2X,1);
        else
            GPSSatOSB(i).C5X = 0;
        end
    end
end