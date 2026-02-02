function [osb_figure] = zero2NaN(osb)
  if isempty(osb)
      osb_figure = [];
  else
      ind = osb == 0; 
      osb(ind) = NaN;
      osb_figure = osb';
  end
end