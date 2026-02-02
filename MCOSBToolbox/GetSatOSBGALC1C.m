% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [GALSatOSB] = GetSatOSBGALC1C(NumSat,EPRN,GALOSB)
    count = 0;
    for i=1:36
        if EPRN(1,i)>0
            count = count + 1;
            OSB = GALOSB(count,1);
        else
            OSB = 0;
        end
        GALSatOSB(i).C1C = OSB;
    end

    count = 0;
    for i=1:36
        if EPRN(1,i)>0
            count = count + 1;
            OSB = GALOSB(count+NumSat,1);
        else
            OSB = 0;
        end
        GALSatOSB(i).C6C = OSB;
    end

    count = 0;
    for i=1:36
        if EPRN(1,i)>0
            count = count + 1;
            OSB = GALOSB(count+NumSat+NumSat,1);
        else
            OSB = 0;
        end
        GALSatOSB(i).C7Q = OSB;
    end

    count = 0;
    for i=1:36
        if EPRN(1,i)>0
            count = count + 1;
            OSB = GALOSB(count+NumSat+NumSat+NumSat,1);
        else
            OSB = 0;
        end
        GALSatOSB(i).C8Q = OSB;
    end

    count = 0;
    for i=1:36
        if EPRN(1,i)>0
            count = count + 1;
            OSB = GALOSB(count+NumSat+NumSat+NumSat+NumSat,1);
        else
            OSB = 0;
        end
        GALSatOSB(i).C5Q = OSB;
    end
end