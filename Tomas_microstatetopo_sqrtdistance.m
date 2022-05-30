%%distance function based on spatial correlation between EEG maps

function D3 = Tomas_microstatetopo_sqrtdistance(XI,XJ)  

% n = size(XI,2);
D3 = sqrt(1-abs(corr(XI',XJ', 'Type','Pearson')));
end

