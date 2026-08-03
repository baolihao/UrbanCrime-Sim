function [A, n, E, t_vec] = system_evolve(B0, n0, tspan, params)
% function [A, n, E, t_vec] = system_evolve(B0, n0, tspan, params)
% input
% output

% M. Zhong (UH)

% this function tries to implement the agent-based model described in the paper
% ``A Statistical Model of Criminal Behavior'' by Bertozzi et al, 2008.
% refer to the Figure 2 on page 1254 for the flow chart, and table 1 for the meaning of the
% parameters

% unpack the parameters
% grid spacing
ell              = params.ell;
% time step size
delta_t          = params.delta_t;
% dynamica attractiveness decay rate
omega            = params.omega;
% measures neighborhood effects ([0, 1])
eta              = params.eta;
% increase in attractiveness due to one burglary event
theta            = params.theta;
% intrinsic (or static) attractiveness of site s
Ast              = params.Ast;
% rate of burglar generation at each site
Gamma            = params.Gamma;
% output skips, how many time steps to skip then save the output
skips            = params.skips;
% compute the number of time steps needed to get to final time
num_steps        = round((tspan(2) - tspan(1))/delta_t);
% check to see if delta_t divides tspan(2) - tspan(1) exactly
if abs(tspan(2) - tspan(1) - delta_t * num_steps) > 1e-15, error(''); end
% number of sites in y-direction
num_rows         = size(B0, 1);
% number of sites in x-direction
num_cols         = size(B0, 2);
% initializing storage
A                = zeros(num_rows, num_cols, 1 + num_steps/skips);
n                = zeros(num_rows, num_cols, 1 + num_steps/skips);
E                = zeros(num_rows, num_cols, 1 + num_steps/skips); 
% the initial time
A(:, :, 1)       = Ast + B0;
n(:, :, 1)       = n0;
% prepare the time loop
Bs_t             = B0;
ns_t             = n0;
% set up the time index counter
t_idx            = 1;
t_vec            = zeros(1, 1 + num_steps/skips);
% now for the time loop
for step_count = 1 : num_steps  
% begin the criminal loop, we might want to parallelize it
% but we have to treat the sites in three groups
% corner site, boundary site, and interior site
% each site will have 4-neighbors (for the infinite domain case)
% the ordering is first from left to right, then from bottom and up
% each site s = (i, j), 1 <= i <= num_rows, 1 <= j <= num_cols
% initialize the moving agents storage
  ns_m           = zeros(num_rows, num_cols);
  ns_g           = zeros(num_rows, num_cols);
  Bs_tp1         = zeros(num_rows, num_cols);
  Es_t           = zeros(num_rows, num_cols);
  for i = 1 : num_rows
    for j = 1 : num_cols
% compute the probability to burgle for each crimnals, and they are independent of each other
% find out if Ast/Gamma is a map or not, and then rectrieve the value only once
      if params.use_map
        Astij    = Ast(i, j);
        Gammaij  = Gamma(i, j);
      else
        Astij    = Ast;
        Gammaij  = Gamma;
      end
      As_t       = Astij + Bs_t(i, j);
% each criminal has a probability of p_s(t) = 1 - e^{-As(t)\delta_t} probability to burglar
% so it's a binormial for a total trials of n_s(i, j) criminals
      ps_t       = 1 - exp(-As_t * delta_t);
% find out where the neighbors are to site s = (i, j)
        neg_idx  = get_neighbor_index(i, j, num_rows, num_cols, params.BC_type);      
% only to committe crime when there is an criminal agent
      if ns_t(i, j) > 0
% 1: to burglar; 0: not to burglar; x is the output vector
        x        = binornd(ns_t(i, j), ps_t);

% compute the number of agents to disappear having committed a burglar (number of 1s)
        n_dis    = nnz(x);
% number of burglary events (each agent committ one)      
        Es_t(i, j) = n_dis;
% n_dis will be used to re-generate the agents later
% compute the number of agents to move
        n_move   = ns_t(i, j) - n_dis;
% save them
        neg_is  = neg_idx(1, :);
        neg_js  = neg_idx(2, :);    
% attractiveness at four neighbors, the order is defined in the function: get_neighbor_indx_PBC
        As_1    = Astij + Bs_t(neg_idx(1, 1), neg_idx(2, 1));
        As_2    = Astij + Bs_t(neg_idx(1, 2), neg_idx(2, 2));
        As_3    = Astij + Bs_t(neg_idx(1, 3), neg_idx(2, 3));
        As_4    = Astij + Bs_t(neg_idx(1, 4), neg_idx(2, 4));  
% total attractivenss
        As_sum  = As_1 + As_2 + As_3 + As_4;
% probabilities to move to neighbors
        qs_1 = As_1/As_sum; qs_2 = As_2/As_sum; qs_3 = As_3/As_sum; qs_4 = As_4/As_sum;
% these four probabilities must add up to one
% the moving is a multinomial distribution for the number of criminals needing to move
% output will be 0/1 with 0: not to move, 1: to move
% the position of the vector indicate the neighbors: 
% each criminal is indepenent of each other
% each criminal having not committedd bulglar, has to move, so there is no probability to stay
        x    = mnrnd(1, [qs_1, qs_2, qs_3, qs_4], n_move);
% this is a loop
        for idx_m = 1 : n_move
% make it into a logical indexing
          idx_i    = neg_is(logical(x(idx_m, :)));
          idx_j    = neg_js(logical(x(idx_m, :)));
          ns_m(idx_i, idx_j) = ns_m(idx_i, idx_j) + 1;
        end
      else
% when on one is at the side s = (i, j), no burglary
        Es_t(i, j) = 0;
      end
% update B
      Laplace_Bst  = (Bs_t(neg_idx(1, 1), neg_idx(2, 1)) + Bs_t(neg_idx(1, 2), neg_idx(2, 2)) ...
                     + Bs_t(neg_idx(1, 3), neg_idx(2, 3)) + Bs_t(neg_idx(1, 4), neg_idx(2, 4)) ...
                     - 4 * Bs_t(i, j))/ell^2;
      Bs_tp1(i, j) = (Bs_t(i, j) + eta * ell^2/4 * Laplace_Bst) * (1 - omega * delta_t) ...
                     + theta * Es_t(i, j);
% re-generate criminal agents, with a rate at gamma
      ps_gen     = 1 - exp(-Gammaij * delta_t);
% just let one criminal to come back
      x            = binornd(1, ps_gen);
      ns_g(i, j)   = nnz(x);
% % binomial probabiliy for each hidden criminals to come back (since n_dis of them has committed
% % crime)
%       if ns_t(i,j) == 0 % no criminal, still has a chance to pop up one
%         x          = binornd(1, ps_gen);
%       else
%         x          = binornd(n_dis, ps_gen);
%       end
% % coming back
%       ns_g(i, j)   = nnz(x);
    end
  end
% update n
  ns_tp1           = ns_g + ns_m;  
% save the output  
  if mod(step_count, skips) == 0       
    t_idx          = t_idx + 1;
    A(:, :, t_idx) = Ast + Bs_tp1;
    n(:, :, t_idx) = ns_tp1;   
    E(:, :, t_idx) = Es_t;
    t_vec(t_idx)   = step_count * delta_t;
    fprintf('\nt = %.2f', step_count * delta_t);
  end  
% save the updates for next
  Bs_t             = Bs_tp1;
  ns_t             = ns_tp1;
end
end