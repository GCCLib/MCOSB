% MCOSB: an open-source software for multi-GNSS and multi-frequency code OSB estimation and analysis
% Authored by Haijun Yuan (navyyuan@yeah.net)
% School of Surveying and Geoinformation Engineering, East China University of Technology (ECUT)
%------------------------------------------------------------------------------------------------------
function [osb,SATG,IDG,SATE,IDE,SATC,IDC] = GetCodeOSB(data,PRN,Signal)
    osb = zeros();
    satG = []; idG = [];
    satE = []; idE = [];
    satC = []; idC = [];    
    for j=1:68
         %GPS
         for i = 1:32    
             if(isempty(data{i,j}))
                data(i,j)={0};
             else 
                data{i,j}(cellfun('isempty',data{i,j})) = {0};
                satG  = [satG;PRN(i)];
                idG   = [idG;Signal(j)];
                len   = length(data{i,j});
                for k = 1:len
                    osb(i,j,k)=data{i,j}{k,:};
                end
             end
         end
         [SATG,~] = unique(satG,'sorted'); 
         [IDG,~]  = unique(idG,'sorted'); 
       %Galileo  
       for i = 33:63
             if(isempty(data{i,j}))
                data(i,j)={0};
             else 
                data{i,j}(cellfun('isempty',data{i,j})) = {0};
                satE = [satE;PRN(i)];
                idE  = [idE;Signal(j)];
                len  = length(data{i,j});
                for k = 1:len
                  osb(i,j,k)=data{i,j}{k,1};
                end
            end
       end
        [SATE,~] = unique(satE,'sorted'); 
        [IDE,~]  = unique(idE,'sorted'); 
       %BDS 
       for i = 64:105
          if(isempty(data{i,j}))
            data(i,j)={0};
          else 
            data{i,j}(cellfun('isempty',data{i,j})) = {0};
            satC = [satC;PRN(i)];
            idC = [idC;Signal(j)];
            len = length(data{i,j});
            for k = 1:len
              osb(i,j,k)=data{i,j}{k,1};
            end
          end 
       end
       [SATC,~] = unique(satC,'sorted'); 
       [IDC,~]  = unique(idC,'sorted'); 
    end
end