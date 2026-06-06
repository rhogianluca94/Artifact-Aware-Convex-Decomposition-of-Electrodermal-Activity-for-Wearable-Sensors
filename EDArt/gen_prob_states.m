function [artifact_vec, pvec, states] = gen_prob_states(expDur, N, fs)

    N = N*fs;
    states = [];
    pvec = [];
    Ns = 7; % number of states
    p0 = 1/(Ns-3);
    prob_art = [.5 .5 0 0 0 0]; % probability of observing a given artifact
%     prob_art = [p0 p0 p0 0 p0 0]; % probability of observing a given artifact
    prob_plateau_pp = [.3 0 0 0 .7 0 0]; % probability to observe a decay after plateau
    prob_plateau_nn = [.3 0 0 0 0 0 .7]; % probability to observe a decay after plateau
    P = [0 prob_art; 1 zeros(1,(Ns-1)); 1 zeros(1,(Ns-1));...
         prob_plateau_pp; 1 zeros(1,(Ns-1));...
         prob_plateau_nn; 1 zeros(1,(Ns-1))]; % transition matrix
    expDurPlateau = 5;
    lambda = [expDur 0 0 expDurPlateau 0 expDurPlateau 0];
    x = 1;              % the initial state is the 0
    nt = 1;             % counter of states transitions
    T0 = 0;
    while T0 < N
        % determine time spent in state X (in samples)
        T = round(exprnd(lambda(x))) * fs;
        if x==1 && T<3*fs         % go on if the clean signal interval lasts less than 1s
            continue;
        end
        T0 = T + T0;
        
        if T0>=N 
            break; 
        end

        % determine state transition
        p = rand(1);
        x = find( p <= cumsum(P(x,:)), 1 );
        if( x > 1 )     % if we transitate to an artifact state, mark its occurrence in time and its type
            % here, generate the artifact and accumulate it
            pvec(nt,1) = T0;
            states(nt,1) = x;
            nt = nt + 1;

            artType = x-1;
            artifact = gen_artifact(T0, N, fs, artType);
            artifact_vec(:, nt) = artifact;
        end
        
    end

end