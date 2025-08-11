% Mapping discovery and linearisation

% Clean workspace
clear all, close all 
clc
format long

% Initializations
n = 6;

% Lagrange Points
L1 = 0.83691513;
L2 = 1.15568217;
L3 = -1.00506265;
L = L1;  % Focal Lagrange point

% Normalization Constants
LU = 389703; TU = 382981;

% Name of file containing initial conditions (ICs) of UPOs
file_path = 'Lyapunov_ICs.csv';
% file_path = 'Halo_ICs.csv';

% Reading in data
init_conds = readmatrix(file_path, 'HeaderLines', 0);

no_ICs = 10;
target = 50;  % Location of desired UPO
LB = target - 5;
UB = target + 5;

init_1_10 = init_conds(LB:UB,2:7);
init = init_1_10;

% Adding augmentation perturbation
mag =  2.5e-7;  % Augmentation perturbation mag, delta v in paper.
init = [init; init_1_10 + [0,0,0,mag,0,0]];
init = [init; init_1_10 + [0,0,0,-mag,0,0]]; 
init = [init; init_1_10 + [0,0,0,0,-mag,0]];
init = [init; init_1_10 + [0,0,0,0,mag,0]];

tim = init_conds(LB:UB,9);  % Period of UPOs
times = [tim;tim;tim;tim;tim];

kfinal = size(init,1);

% Initialising vars
Psec = []; PsecNext = [];

T_des = init_conds(target,9);
map_time = T_des;
coords_des = init_conds(target,2:7);
coords_ori = coords_des;

options = odeset('RelTol',1e-12,'AbsTol',1e-14*ones(1,6));
%% S2 section - uncomment for map disc at S2
% [t,sol] = ode87(@ThreeBody,[0 T_des/2-0.001],coords_ori,options); % THIS IS WRONG
% [t1_z,x0_z] = find0(sol(end,:), t(end));
% coords_des = x0_z;

%% Data Collection
temp = [];
count = 1;

eta = 1.0;  % User-defined (Leave at 1.0)
no_sim_cross = 1;  % Number of periods to simulate for (1 for S1, 2 for S2)
for k = 1:kfinal

    % Generate trajectory starting from each different IC
    [t,sol] = ode87(@ThreeBody,[0 times(k)*no_sim_cross + 1],init(k,:),options);
    disp(['Simulating UPO no ', num2str(k)])
    
    % Adding IC as first crossing point 
    rel_init = init(k,:) - coords_des;
    if norm(rel_init) <= eta
        temp(count,:) = rel_init;
        count = count + 1;
    end

    for j = 2:length(t)-1  % Looping through all times

        if (sol(j,2) < 0 && sol(j+1,2) >= 0)  % S1
        %if (sol(j,2) > 0 && sol(j+1,2) <= 0)  % S2
            
            [t1_z,x0_z] = find0(sol(j,:), t(j));
            if norm(x0_z - coords_des) <= eta
                disp(['Adding crossing ', num2str(count)])
                temp(count,:) = x0_z - coords_des;  % nth crossing
                count = count + 1;
            end
        end
    end

    Psec = [Psec; temp(1:size(temp,1)-1,:)];
    PsecNext = [PsecNext; temp(2:size(temp,1),:)];
   	count = 1;
    temp = [];
end

xt = Psec; % x_n
xtnext = PsecNext; % x_n+1

%% Sparse Regression
polyorder = 5;  % polynomial order 
usesine = 0; 
Theta = poolData(xt,n,polyorder,usesine); 

lambda = 1e-6;  % Sparsity promoting param

Xi = sparsifyDynamics(Theta,xtnext,lambda,n);
if n == 6
[yout, newout] = poolDataLIST({'x','y','z','xd','yd','zd'},Xi,n,polyorder,usesine);
elseif n == 3
[yout, newout] = poolDataLIST({'x','y','z'},Xi,n,polyorder,usesine);
elseif n == 2
 [yout, newout] = poolDataLIST({'x','y'},Xi,n,polyorder,usesine);
elseif n == 1 
  [yout, newout] = poolDataLIST({'x'},Xi,n,polyorder,usesine);
end 

fprintf('SINDy model: \n ')
counting = 0;
for k = 2:size(newout,2) 
    SINDy_eq = newout{1,k}; 
    SINDy_eq = [SINDy_eq  ' = '];
    strings = '';
    new = 1;
    for j = 2:size(newout, 1)
       if newout{j,k} ~= 0 
           if new == 1 
             SINDy_eq = [SINDy_eq  num2str(newout{j,k}) newout{j,1} ]; 
             add = num2str(newout{j,k});
             add2 = strcat('*',string(newout{j,1}));
             strings = strcat(strings,strcat(add, add2));
             new = 0;
           else 
             SINDy_eq = [SINDy_eq  ' + ' num2str(newout{j,k}) newout{j,1} ' '];
             add = strcat(' + ', num2str(newout{j,k}));
             add2 = strcat('*',string(newout{j,1}));
             strings = strcat(strings, strcat(add,add2));
           end 
       end
    end
  fprintf(SINDy_eq)
  fprintf('\n ')
  if counting == 0
      eqns = str2sym(strings);
  else
      eqns = [eqns; str2sym(strings)];
  end
  counting = counting + 1;
end

%% Computing Jacobian at xbar
syms x y z xd yd zd xx yy zz xdxd ydyd zdzd
eqns = subs(eqns, [xx, yy, zz, xdxd, ydyd, zdzd], [x^2, y^2,z^2, xd^2, yd^2, zd^2]);

A = jacobian(eqns,[x y z xd yd zd]); jac = A;
A = double(subs(A, [x, y, z, xd, yd, zd], [0 0 0 0 0 0]))
A_old = A;

% Ensures map is 6 dim for control purposes
Anew = [A(1:2,:); zeros(1,6); A(3:4,:); zeros(1,6)]; 
A = Anew;


%% Restricted three-body ODE right-hand-side
function dx = ThreeBody(t,x)

    % Mass ration parameter
    mu = 1.215058560962404e-2;
    
    d = sqrt((x(1) + mu)^2 + x(2)^2 + x(3)^2);
    r = sqrt((x(1) -1 + mu)^2 + x(2)^2 + x(3)^2);

    dx(1) = x(4);
    dx(2) = x(5);
    dx(3) = x(6);
    dx(4) = 2*x(5) + x(1) - (1-mu)*(x(1)+mu) / d^3 - mu * (x(1) - 1 + mu) / r^3;
    dx(5) = -2*x(4) + x(2) - (1-mu)*x(2) / d^3 - mu * x(2) / r^3;
    dx(6) = -(1-mu) * x(3)/ d^3 - mu * x(3) / r^3;
end

function [t1_z, x1_z] = find0(x0, tGuess)
% Find where y crosses 0
% Based on code by Shane Ross, 2024, "cr3bp matlab"

% tols & options
tolzero = 1e-12;
options = optimset('DiffMaxChange', tolzero, 'TolCon', tolzero, 'TolFun', tolzero);
optionsInt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12 * ones(1,6));

t0_z = tGuess;
x0_zz = x0;

function y = orbY_local(t1)
    if t0_z == t1
        fin = x0_zz;
    else
        [~, sol] = ode87(@ThreeBody, [t0_z t1], x0_zz, optionsInt);
        fin = sol(end, :);
    end
    y = fin(2);  
end

t1_z = fzero(@orbY_local, t0_z, options);

% compute state at found time
if t0_z == t1_z
    x1_z = x0_zz;
else
    [~, sol] = ode87(@ThreeBody, [t0_z t1_z], x0_zz, optionsInt);
    x1_z = sol(end, :);
end
end
