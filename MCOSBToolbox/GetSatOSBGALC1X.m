% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%-----------------------------------------------------------------------------------------------------
function [GALSatOSB] = GetSatOSBGALC1X(NumSat,EPRN,GALOSB)
    count = 0;
    for i=1:36
        if EPRN(1,i)>0
            count = count + 1;
            OSB = GALOSB(count,1);
        else
            OSB = 0;
        end
        GALSatOSB(i).C1X = OSB;
    end

    count = 0;
    for i=1:36
        if EPRN(1,i)>0
            count = count + 1;
            OSB = GALOSB(count+NumSat,1);
        else
            OSB = 0;
        end
        GALSatOSB(i).C6X = OSB;
    end

    count = 0;
    for i=1:36
        if EPRN(1,i)>0
            count = count + 1;
            OSB = GALOSB(count+NumSat+NumSat,1);
        else
            OSB = 0;
        end
        GALSatOSB(i).C7X = OSB;
    end

    count = 0;
    for i=1:36
        if EPRN(1,i)>0
            count = count + 1;
            OSB = GALOSB(count+NumSat+NumSat+NumSat,1);
        else
            OSB = 0;
        end
        GALSatOSB(i).C8X = OSB;
    end

    count = 0;
    for i=1:36
        if EPRN(1,i)>0
            count = count + 1;
            OSB = GALOSB(count+NumSat+NumSat+NumSat+NumSat,1);
        else
            OSB = 0;
        end
        GALSatOSB(i).C5X = OSB;
    end
end