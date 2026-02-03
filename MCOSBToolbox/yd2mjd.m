function out_mjd = yd2mjd(year,doy)
%% convert days in year to year-month-day
% modified from: https://blog.csdn.net/qq_41696018/article/details/119647098
    mjd52 = 34012;
    yr = year;
    if(year<1900)
        if(year>60)
            yr = yr+1900;
        else
            yr = yr+2000;
        end
    end
    yr = yr - 1900;
    nyear = yr - 52;
    leap = (nyear+3)/4;
    nday = nyear*365;
    out_mjd = (mjd52-1)+nday+leap+doy;
    out_mjd = fix(out_mjd);
end