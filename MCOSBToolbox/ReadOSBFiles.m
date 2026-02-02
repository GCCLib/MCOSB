% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [data]=ReadOSBFiles(path,PRN,Signal)
    data = cell(105,68);
    if path == 0
        disp('No folder selected.');
        return;
    end
    fileListing = dir(fullfile(path,'*.BIA'));
    len = length(fileListing);
    i = 1;
    while i < len+1
       FileName = fullfile(path,fileListing(i).name);
       fid = fopen(FileName,'r');
         while ~feof(fid)
            tline = fgetl(fid);
            if contains(tline,'+BIAS/SOLUTION')
                % disp(tline)
                break
            end    
         end
         while ~feof(fid)   
             kline = fgetl(fid);
             if (contains(kline,'OSB'))
                if (kline(16:19)=="    ")
                      prn = kline(12:14);
                      signal = kline(27:28);
                      j = find(PRN == prn);
                      if isempty(j)
                          continue;
                      end
                      k = Signal == signal;
                      type = kline(26);
                      % code
                      if (contains(type,'C'))
                          data{j,k}{i,1} = str2double(kline(70:91)); %value
                      end
                 end
            end
        end
       fclose(fid);  
       i = i + 1;
    end
end