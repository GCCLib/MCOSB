% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
function [GPSSatOSB] = GetSatOSBGPSC5Q(SatC1C,SatC1W,SatC2W,SatC2L,GPRNC1C,GPRNC1W,GPRNC2W,GPRNC2L,GPRNC5Q,GPSOSB)
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
        if GPRNC2L(1,i)>0
            count = count + 1;
            GPSSatOSB(i).C2L = GPSOSB(count+SatC1C+SatC1W+SatC2W,1);
        else
            GPSSatOSB(i).C2L = 0;
        end
    end

    count = 0;
    for i=1:32
        if GPRNC5Q(1,i)>0
            count = count + 1;
            GPSSatOSB(i).C5Q = GPSOSB(count+SatC1C+SatC1W+SatC2W+SatC2L,1);
        else
            GPSSatOSB(i).C5Q = 0;
        end
    end
end