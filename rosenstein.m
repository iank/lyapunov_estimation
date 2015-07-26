function [di] = rosenstein(tss, T, name)
    %% Rosenstein et al / Lyapunov exponents from small data sets

    %% Estimate lag and mean period using the FFT
    % Lag value calculated from MI (see fraser / swinney, "independent
    % coordinates for strange attractors from mutual information" )

    % mean period
    Fs = 204800;

    %TODO: get this from FFT
    % dominant frequency 100 Hz
    mean_period = 1/100; % seconds
    mean_lag = floor(mean_period*Fs);

    %% Reconstruct attractor dynamics using method of delays

    N = length(tss(:,1));
    x1 = tss(1:(N-8*T-1));
    x2 = tss((1:N-8*T-1)+T);
    x3 = tss((1:N-8*T-1)+2*T);
    x4 = tss((1:N-8*T-1)+3*T);
    x5 = tss((1:N-8*T-1)+4*T);
    x6 = tss((1:N-8*T-1)+5*T);
    x7 = tss((1:N-8*T-1)+6*T);
    x8 = tss((1:N-8*T-1)+7*T);
    x9 = tss((1:N-8*T-1)+8*T);


    X = [x1 x2 x3 x4 x5 x6 x7 x8 x9];

    %% Find nearest neighbors, constrain temporal separation

    % TODO: can speed this up quite a bit using pdist() & then the temporal constraint.
    % - actually pdist uses too much memory, subsample X

    % How many i can we use?
    max_i = 40960;

    % Number of points to use
    npts = 10000;
    nn = [];
    samp = randperm(size(X,1) - max_i);
    samp = samp(1:npts);


    for k=samp  % for each point in X
        x = X(k,:);
        
        dists = sqrt(sum(abs(repmat(x, [size(X,1) 1]) - X).^2,2));
        % find set of points w/ temporal separation greater than 1 mean period
        unacceptable_idxs = ((abs((1:length(dists)) - k) < mean_lag) | (1:length(dists) > length(dists) - max_i));
        dists(unacceptable_idxs) = max(dists);
        nn_idx = find((dists == min(dists))');
        
        nn = [nn; nn_idx];
    end

    %%

    di = zeros(max_i, 1);
    for ii=1:max_i
        % Compute distance between nns
        dists = sqrt(sum(abs(X(samp+ii,:) - X(nn+ii,:)).^2,2));
        di(ii) = mean(dists);
    end

end
