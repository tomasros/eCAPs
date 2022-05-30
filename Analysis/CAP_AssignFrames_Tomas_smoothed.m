function [i] = CAP_AssignFrames_Tomas_smoothed(Cp,XON,d,T)

    % Value of correlation below which 5% of all values lie (so gives an
    % estimate of a 'bad correlation for a trace belonging to a cluster'
    CT = prctile(d,T);

    % cluster assignment in young
    r = corr(Cp',XON);
    [c,i] = max(r);
    
    
    
     % option4: windowed smoothing via Poulsen microstate toolbox
% opts.b=2; %%smoothing window in samples (@40Hz 1 sample=25ms)
% opts.lambda=5; %smoothing weight
% i = MicroSmooth(XON, Cp', 'windowed' ,opts); %% 'reject segments' (default) or 'windowed'.
% i=i';

    
    

    % If the correlation value is too low (below threshold), the index is set
    % to a new non-existin,g group
    i(c<CT) = size(Cp,1)+1;
    

    
end