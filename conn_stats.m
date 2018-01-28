
% K6.
e = [[1 2]; [1 3]; [1 4]; [1 5]; [1 6]; [2 3]; [2 4]; [2 5]; [2 6]; [3 4]; [3 5]; [3 6]; [4 5]; [4 6]; [5 6]];
u = 1;
v = 6;

% e = [[1 2]; [1 3]; [2 5]; [3 4]; [2 4]; [3 5]; [4 6]; [5 6]; [2 3]; [1 6]];
% u = 1;
% v = 6;


% e = [[1 3]; [3 2]; [1 2]; [1 4]; [1 5]; [5 2]];
% u = 1;
% v = 2;

n = max(max(e));
m = size(e,1);

num_trials = 1000;

poolobj = parpool;
conn_probs = zeros(1,num_trials);
exp_corr2s = zeros(1,num_trials);
parfor trialnum = 1:num_trials
  %  trialnum
  if mod(trialnum,100) == 0
      trialnum
  end

    % Initialize the channels randomly.
    a = rand(1,m);
    b = rand(1,m); % Asymmetric case.
   %  b = 1-a; % Symmetric case. 
   
    % Calculate probability $u$ and $v$ are connected.
    % Each of the edges is on with probability equal to its SNR.
    snr = (1/2)*(a-b).^2 .* (1./(a+b) + 1./(2-a-b));
    
    % TODO: Can calculate this faster using monotonicity of 
    conn_prob = 0;
    not_conn_prob = 0;
    for i=0:(2^m-1)
%         if mod(i,1000) == 0
%             i
%             conn_prob
%         end
        estate = dec2bin(i,m);
        edgin = find(estate == '1');
        edgnotin = find(estate == '0');
        prob_estate = prod(snr(edgin))*prod(1-snr(edgnotin));
        
        G = graph(e(edgin,1), e(edgin,2),1,n);
        [~,d] = shortestpath(G,u,v,'Method','unweighted');
        isconn = ~isinf(d);
        if isconn
            conn_prob = conn_prob + prob_estate;
        else
            not_conn_prob = not_conn_prob + prob_estate;
        end
    end
    conn_probs(trialnum) = conn_prob;
    assert(abs(conn_prob + not_conn_prob - 1) < 1e-12);

    % Calculate expected correlation between $u$ and $v$, given edge
    % labels.
    % E[E[uv | edges]^2]
    exp_corr2 = 0;
    for i=0:(2^m-1)
%         if mod(i,1000) == 0
%             exp_corr2
%             i
%         end
        estate = dec2bin(i,m);
%         edgin = find(estate == '1');
%         edgnotin = find(estate == '0');
%         prob_estate = prod(snr(edgin))*prod(1-snr(edgnotin));
        
        exp_corr_contrib = 0;
        probestate = 0;
        for j=0:(2^n-1)
            vstate = dec2bin(j,n);
            p_contrib = 1;
            for k = 1:m
                if vstate(e(k,1)) == vstate(e(k,2))
                    if estate(k) == '0'
                        p_contrib = p_contrib * a(k);
                    else
                        p_contrib = p_contrib * (1-a(k));
                    end
                else
                    if estate(k) == '0'
                        p_contrib = p_contrib * b(k);
                    else
                        p_contrib = p_contrib * (1-b(k));
                    end
                end
            end
            if vstate(u) == vstate(v)
                sgn = 1;
            else
                sgn = -1;
            end
            corr_contrib = p_contrib * sgn / 2^n;
            probestate = probestate + p_contrib / 2^n;
            
            exp_corr_contrib = exp_corr_contrib + corr_contrib;
        end
        cond_corr_contrib = exp_corr_contrib / probestate;
        
        exp_corr2 = exp_corr2 + probestate * cond_corr_contrib^2;
    end
    exp_corr2s(trialnum) = exp_corr2;
    run_vals = [conn_prob exp_corr2]
end

save('collected_data.mat','exp_corr2s','conn_probs');
assert(all(exp_corr2s <= conn_probs));

% hold off
% scatter(exp_corr2s,conn_probs,'.')
% hold on
% ezplot(@(x) x, [0 1])
% disp('Expected correlation seems <= conn_prob.')
