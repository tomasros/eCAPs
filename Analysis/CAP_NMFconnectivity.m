function [i] = CAP_NMFconnectivity(Cp,XON, xindp1)

    % Value of correlation below which 5% of all values lie (so gives an
    % estimate of a 'bad correlation for a trace belonging to a cluster'
%     CT = prctile(d,T);

    % cluster assignment in young
%     r = corr(Cp',XON);
%     [c,i] = max(r);

    % If the correlation value is too low (below threshold), the index is set
    % to a new non-existing group
%     i(c<CT) = size(Cp,1)+1;

% Number of subjects
n = size(xindp1,2);

for i = 1:n
    
for time=1:size(XON,2)
% time_course(:,time)=nnls(CP', XONn_filtered(:,time));
time_course(:,time)=tntnn(Cp', XON(:,time));
end

end


end