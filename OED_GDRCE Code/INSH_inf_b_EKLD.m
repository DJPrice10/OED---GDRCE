% clear

myCluster = parcluster('local');
myCluster.NumWorkers = 4;
saveProfile(myCluster);
                           
rng(1)

totaltime = tic;

number_of_data=100000;

% % % % % % % % % % % % % % % % % % % % % % % % %
% % % Define prior distribution parameters% % % %
% % % % % % % % % % % % % % % % % % % % % % % % %

% Inputs for prior distributions for each parameter
alpha_prior_params = [0.15, 0.25];
delta_prior_params = [4.825, 0.25];
beta_prior_params = [0.7, 0.25];

% End points of parameter space to consider
upper_limits=[alpha_prior_params(2), 300, 6.25];
lower_limits=[alpha_prior_params(1), 35, 0.25];

% % % % % % % % % % % % % % % % % % % % % %
% % % Define experimental parameters% % % %
% % % % % % % % % % % % % % % % % % % % % %

% Design specific parameters

nchickens=40;
maxngroups = 5;

% Shape and rate parameters for the Erlang distributed latent period. 
% nclasses must be integer.
exposure_rate=2;
nclasses=2;


% Grid spacing on design space
dose_stepsize = 0.05;
time_stepsize= 0.05;

% Range of design choices
dose_limits = [dose_stepsize, 8];
time_limits = [time_stepsize, 6];

% ABC tolerance
tol=.25;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % %  Define prior distribution across grid points % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Grid (alpha, delta)

n=100;
MAX_XY=upper_limits;
MIN_XY=lower_limits;

scaling=MAX_XY-MIN_XY;
% Get same values used to evaluate posterior
[X,Y]=meshgrid(MIN_XY(1):scaling(1)/(n-1):MAX_XY(1),MIN_XY(2):scaling(2)/(n-1):MAX_XY(2));

Xctrs = X(1,:);
Yctrs = Y(:,1);

% Evaluate prior at those values
% Xprior = unifpdf(X,alpha_prior_params(1),alpha_prior_params(2));
Xprior = 1/(MAX_XY(1)-MIN_XY(1)) * ones(size(X));
Yprior = lognpdf(Y,delta_prior_params(1),delta_prior_params(2));



% Sample param indices according to weights (X/Y/Zprior)
alpha_idx = randsample(1:length(Xctrs), number_of_data, true, Xprior(1,:));
delta_idx = randsample(1:length(Yctrs), number_of_data, true, Yprior(:,1));

% Get param values corresponding to indices
alpha_params = Xctrs(alpha_idx)';
delta_params = Yctrs(delta_idx);


prior_values = [alpha_params, delta_params];

% Random beta values
zvals = MIN_XY(3):scaling(3)/(n-1):MAX_XY(3);
Zprior = lognpdf(zvals,beta_prior_params(1),beta_prior_params(2));
beta_idx = randsample(1:length(zvals), number_of_data, true, Zprior);
beta_params = zvals(beta_idx);


% Evaluate prior grid, as EKLD estimate is w.r.t. beta
prior_grid = hist(beta_params, zvals);
prior_grid = prior_grid./sum(prior_grid);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % Set INSH convergence parameters % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Number of iterations
W=60;

% Iterations at which to increase exploitation
stepdown = [30, 15, 15];

% if (sum(stepdown)~=W)
%     fprintf('Error: Make sure sum stepdown equals W. \n')
%     return;
% end

% No. of new designs to sample, and to keep, at each iteration
m=[3*ones(stepdown(1),1); 6*ones(stepdown(2),1); 15*ones(stepdown(3),1)];
ntokeep = [150*ones(stepdown(1),1); 75*ones(stepdown(2),1); 30*ones(stepdown(3),1)];

% Reduce standard deviation to hone in on designs around optimal region
time_sample_sd= [0.1*ones(stepdown(1),1); 0.075*ones(stepdown(2),1); 0.05*ones(stepdown(3),1)];
dose_sample_sd= [0.1*ones(stepdown(1),1); 0.075*ones(stepdown(2),1); 0.05*ones(stepdown(3),1)];

% If want fixed standard deviation throughout, use =sigma*ones(W,1);

%%
% Initial number of designs to be considered at each of 2,3,4,5 group
% sizes. Must be maxngroups-1 in length.
ninitialdesignspergroupsize = [50, 100, 250, 500];


% Randomly generate the initial designs.
current_wave_designs = [];

for i=2:maxngroups
    current_wave_designs = [current_wave_designs; zeros(ninitialdesignspergroupsize(i-1), maxngroups-i) repmat(floor(nchickens/i), ninitialdesignspergroupsize(i-1), i)  zeros(ninitialdesignspergroupsize(i-1), maxngroups-i)  round(unifrnd(dose_limits(1)*ones(ninitialdesignspergroupsize(i-1),i), dose_limits(2))/dose_stepsize)*dose_stepsize   zeros(ninitialdesignspergroupsize(i-1), maxngroups-i)   round(unifrnd(time_limits(1)*ones(ninitialdesignspergroupsize(i-1),i), time_limits(2))/time_stepsize)*time_stepsize   ];
end




% Order the designs according to the dose. Re-order the obs. times 
% according to the dose ordering. 

for i=1:size(current_wave_designs,1)
    [B, I] = sort(current_wave_designs(i, (maxngroups+1):(2*maxngroups)));
    current_wave_designs(i, (maxngroups+1):(2*maxngroups)) = B;
    
    ts = current_wave_designs(i, (2*maxngroups+1):(3*maxngroups));
    current_wave_designs(i, (2*maxngroups+1):(3*maxngroups)) = ts(I);
    
end



keep_all = [];


for w=1:W 
    sprintf('Wave: %u', w)
   
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % Simulate data at each of the current unique designs % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    % Manipulate the data set so that I can simulate easily. Split out each
    % group independently, and simulate the groups. Then, stick them back
    % together when evaluating the utility for each design
    
    groupsizes = current_wave_designs(:,1:maxngroups);
    doses = current_wave_designs(:,(maxngroups+1):(2*maxngroups));
    times = current_wave_designs(:,(2*maxngroups+1):(3*maxngroups));
    
    allgroupdesigns = [groupsizes(:), doses(:), times(:)];
    % Retain only those that aren't all zero
    allgroupdesigns( ~any(allgroupdesigns,2), : ) = [];
    
    % Get only unique designs
    allgroupdesigns = unique(allgroupdesigns, 'rows');
    
%     sprintf('No. unique group designs: %u', size(allgroupdesigns,1))
    
    % Get the unique number of chickens and doses to create initial number
    % of exposed chickens according to these values
    unique_n_d = unique(allgroupdesigns(:,1:2), 'rows');
    
    % Dose only comes into the equation here, the initial number of "exposed"
    prob_inf = (1-bsxfun(@power, 1+bsxfun(@ldivide, delta_params, 10.^unique_n_d(:,2)'), -alpha_params));
    initially_exposed_data = binornd(repmat(unique_n_d(:,1)',number_of_data, 1), prob_inf);
    
    clear('prob_inf')
    
    data = zeros(size(prior_values,1), size(allgroupdesigns,1));
        
    % Forward simulate each of the groups according to SEKI dynamics with
    % transmission_rate in beta_params(j). 
    for i=1:size(unique_n_d,1)
        groupsize = unique_n_d(i,1);
        
        ind = ismember(allgroupdesigns(:,1:2), unique_n_d(i,:), 'rows')';
        possibletimes = allgroupdesigns(ind, 3);
        
        tmpdat = zeros(number_of_data, sum(ind));
        parfor(j=1:size(prior_values,1), 4)
                tmpdat(j,:) = getdataSEKI(beta_params(j), exposure_rate, groupsize, initially_exposed_data(j,i), possibletimes, nclasses);
        end
        
        data(:,ind) = tmpdat;   
    end
    
    clear('initially_exposed_data', 'tmpdat')
    
    % Data now has a column corresponding to each unique group in the original
    % proposed designs. Need to now stitch these back together in some way when
    % evaluating the posterior distribution for each
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % Evaluate utility at each design using ABCdE approach  % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    kld = zeros(size(current_wave_designs,1),1);

    parfor(i=1:size(current_wave_designs,1), 4)
        tmpdesign = [current_wave_designs(i,1:maxngroups)' current_wave_designs(i,(maxngroups+1):(2*maxngroups))' current_wave_designs(i,(2*maxngroups+1):(3*maxngroups))'];
        tmpdesign(~any(tmpdesign,2),:)=[];
        
        % locb gives the rows of allgroupdesigns that together make current
        % design. These also correspond to the columns that will then be used to
        % form the observed data used for evaluating the posterior distribution
        [~, locb] = ismember(tmpdesign, allgroupdesigns, 'rows');
        
        tmp_data=data(:,locb);
        
        % Grab unique data sets, and number of times each occurs
        [uni_data, num, ~] = howmanyunique(tmp_data);
        
        tmp_kld=zeros(size(uni_data,1),1);
        
        % Get std dev's for ABC discrepancy
        std_tmp = std(tmp_data);
        
        for k=1:size(uni_data,1)
            keep = ABCDE_chunk(tmp_data, uni_data(k,:), std_tmp)<(tol*size(tmpdesign,1));
                       
            xr = beta_idx(keep)';

            lendens = length(xr);
            % Count number of each index
            % [URP, numURP] = howmanyunique([xr,yr]);
            [URP, numURP] = howmanyunique(xr);
            % Give index of prior values
            % pos = sub2ind([n,n],URP(:,1),URP(:,2));
            pos = URP;
            % Normalise
            denspos = numURP/lendens;
            % Evaluate KLD
            tmp_kld(k) = (log(denspos) - log(prior_grid(pos)'))'*denspos;
        end
        % Evaluate utility for design i, as the average utility contribution
        % from each unique data set
        kld(i) = dot(tmp_kld, num/sum(num));
        
        
    end
    
    clear('data')
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % Store all sampled design points and corresponding utility % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    keep_all = [keep_all; [current_wave_designs, kld, w*ones(length(kld),1)]];
    
    % Acceptance Criteria
    sortklds = sort(kld,'descend');
    accept_cutoff = sortklds(ntokeep(w));
    update_wave = current_wave_designs(kld >= accept_cutoff, :);
    
    % Include the optimal design thus far in the update wave so that we
    % keep exploring the space around the optimal, rather than hitting it
    % then moving away from it
    
    if(keep_all(keep_all(:,3*maxngroups+1)==max(keep_all(:,3*maxngroups+1)),end)~=w)
        update_wave = [update_wave; keep_all((keep_all(:,3*maxngroups+1)==max(keep_all(:,3*maxngroups+1))),1:3*maxngroups)];
    end

    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % Sample new designs around accepted designs % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    nonzero_indicator = update_wave>0;
    ngroupsperdesign = sum(nonzero_indicator,2)/3;
    
    current_wave_designs = [];
    for k=1:size(update_wave,1)
        current_wave_designs = [current_wave_designs; zeros(m(w), maxngroups-ngroupsperdesign(k)),  repmat(floor(nchickens/ngroupsperdesign(k)), m(w), ngroupsperdesign(k)), zeros(m(w), maxngroups-ngroupsperdesign(k)),  discrete_trunc_sample(update_wave(k,(2*maxngroups-ngroupsperdesign(k)+1):(2*maxngroups)), dose_limits, dose_sample_sd(w)*eye(ngroupsperdesign(k)) , dose_stepsize, m(w)),...
            zeros(m(w), maxngroups-ngroupsperdesign(k)),  discrete_trunc_sample(update_wave(k,(3*maxngroups-ngroupsperdesign(k)+1):(3*maxngroups)), time_limits, time_sample_sd(w)*eye(ngroupsperdesign(k)) , time_stepsize, m(w)) ];
    end
    
    % Add OED back in to be re-evaluated at next iteration
    current_wave_designs = [current_wave_designs; keep_all((keep_all(:,3*maxngroups+1)==max(keep_all(:,3*maxngroups+1))),1:3*maxngroups)];
    

    for i=1:size(current_wave_designs,1)
        [B, I] = sort(current_wave_designs(i, (maxngroups+1):(2*maxngroups)));
        current_wave_designs(i, (maxngroups+1):(2*maxngroups)) = B;
        
        ts = current_wave_designs(i, (2*maxngroups+1):(3*maxngroups));
        current_wave_designs(i, (2*maxngroups+1):(3*maxngroups)) = ts(I);
        
    end
end

runtime = toc(totaltime);

% Create a final column containing the number of groups
ng = sum(keep_all(:,1:maxngroups)>0, 2);
keep_all = [keep_all, ng];
%


save('Output_inf_b_EKLD.mat')



