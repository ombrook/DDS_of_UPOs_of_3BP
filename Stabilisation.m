% Satellite control script for 1 cycle orbits

% Control constraint matrix
R = zeros([6,3]); R(4,1) = 1; R(5,2) = 1; R(6,3) = 1;

% Lagrange Points - As from JPL
L1 = 0.83691513;
L2 = 1.15568217;
L3 = -1.00506265;
L = L1; % Focal Lagrange point
LU = 389703;
TU = 382981;


dt = 0.00001;
tspan = [0 map_time+2];
K = K_sat;
t0 = 0;
xc = [];
tc = [];


x0 = coords_des;
options = odeset('RelTol',1e-12,'AbsTol',1e-14*ones(1,6));
[t,x] = ode87(@ThreeBody,[0 map_time + 0.1],x0,options);
xc = [];
tc = [];
c_hist =[];
impulse =[];
t_imp = [];
nu =[]; ps =[]; qs=[];

% Loop to detect first crossing of axis
for j = 2:size(x,1)-1
    if (x(j,2) < 0 && x(j+1,2) >= 0)
        disp('Crossing Detected')
        xc = [xc; x(1:j-1,:)];
        tc = [tc; t(1:j-1)];

        [t1_z,x0new] = find0(x(j,:), t(j));
        [ts,sol] = ode87(@ThreeBody,[t(j) t1_z],x0new,options);
        
        xc = [xc; sol(1:end-1,:)];
        tc = [tc; ts(1:end-1)];
        
        % Check inside eta
        if norm(x0new - coords_des) <= eta
            disp('Within control region ')
            disp(['Norm = ' num2str(norm(x0new - coords_des))])
            c_hist = [c_hist, R*K_sat*(x0new - coords_des)'];
            impulse = [impulse; x0new];
            t_imp = [t_imp; t(j)];
            x0new = (x0new' + R*K_sat*(x0new - coords_des)')';
        else
            disp('First passing not within control region')
            disp(['Norm = ' num2str(norm(x0new - coords_des))])
        end

        break
    end

end


%%
idx = 5; % Specify no. of periods to simulate for
for p = 1:idx
    % New simulation for around one period
    [t,x] = ode87(@ThreeBody,[t1_z t1_z + map_time*2],x0new, options);
    
    disp('Starting new sim')

    % Check for axis crossing
    for j = 2:length(x(:,1))-1
        if (x(j,2) < 0 && x(j+1,2) >= 0)
            disp('Crossing Detected')
            xc = [xc; x(1:j-1,:)];
            tc = [tc; t(1:j-1)];
    
            [t1_z,x0new] = find0(x(j,:), t(j));
            [ts,sol] = ode87(@ThreeBody,[t(j) t1_z],x0new,options);
            % Select sim value closer to y=0 section
            
            xc = [xc; sol];
            tc = [tc; ts];
            tc_old = tc(end);

            % Check inside eta
            if norm(x0new - coords_des) <= eta
                disp('Within control region, difference =')
                disp(x0new - coords_des)
                impulse = [impulse; x0new];
                t_imp = [t_imp; tc_old];
                c_hist = [c_hist, R*K_sat*(x0new - coords_des)'];
                x0new = (x0new' + R*K_sat*(x0new - coords_des)')';
                
            else
                disp('Outside of control region')
                disp( ['Norm = ' num2str(norm(x0new - coords_des))])
            end

            break % Breaks out of j loop
        end
    end
end 

tuspan = tc(end);
% coords_ori specified in Mapping_discovery.m
[tu, xu] = ode87(@ThreeBody,[0 tuspan],coords_ori,options); % Sim uncontrolled

blue = "#0000a2";
red = "#bc272d";
yellow = "#e9c716";
teal = "#50ad9f";


fsize = 15;

idxs = tc < map_time;
mult = round(tc(end)/map_time);
des = xc(idxs,:); tcd=tc(idxs);
for i = 1:mult-1
    des = [des; xc(idxs,:)];
    tcd = [tcd; tcd(end)+ tc(idxs)];
end

%% Proper x vs t
fsize = 25;
f1 = figure; hold on;
set(gcf,'renderer','Painters')
plot(tu, xu(:,1),'Color',red,'LineStyle','-','LineWidth',2,'DisplayName','Uncontrolled');
plot(tc,xc(:,1),'Color',blue,'LineStyle','-','LineWidth',2,'DisplayName','Controlled');

%plot(tcd,des(:,1),'--','LineWidth',1,'DisplayName','Desired','Color',	"#EDB120"); 

plot(t_imp, impulse(:,1), 'Color', teal, ...
    'Marker', '.', 'MarkerSize', 25, 'LineStyle', 'none', 'DisplayName', 'Impulse');
hold off;
legend('show', 'FontSize', fsize, 'Interpreter', 'latex', 'location', 'southwest');
xlabel('t [TU]','Interpreter','latex', 'FontSize', fsize); ylabel('x [LU]','Interpreter','latex', 'FontSize', fsize);
set(gca, 'TickLabelInterpreter', 'Latex')
ax = gca;
ax.XLim = [0 100];

ylim(ax, [-0.5, 1.5])

ax.XAxis.FontSize = fsize;
ax.YAxis.FontSize = fsize;
box on
grid minor

% Set the figure size in centimeters
width = 30; % Width in cm
height = 12; % Height in cm
set(f1, 'Units', 'centimeters');
set(f1, 'Position', [10, 10, width, height]);


%% 

[t,x] = ode87(@ThreeBody,[0 map_time],coords_des,options);
offsets = impulse(1,:);
del_v = c_hist(4:6,:) * LU/TU * 1000; % Converting to m/s

del_v_cost = sum(abs(del_v(1:2,1:14)), 'All') % Calculates cost

%% Eigenvectors plot
del_v = c_hist(4:6,:) * LU/TU * 1000;
v_hist = c_hist(4:5,:);
v_hist = normc(v_hist);

% Calculating Monodromy matrix - PHIgetHT code is modified from Shane Ross
[~,~,M,Phi] = PHIgetHT(coords_des, map_time);
[V,D] = eigs(M);

stab_eig1 = V(:,5);
stab_eig2 = V(:,6);
pos_stab1 = normc(stab_eig1(1:2)); 
pos_stab2 = normc(stab_eig2(1:2));
ustab_eig2 = V(:,1);
pos_ustab2 = normc(ustab_eig2(1:2));

figure; % Plotting the delta_v vectors
set(gcf,'renderer','Painters')
hold on;

% Local stable manifold
plot(offsets(1)+[-pos_stab2(1),0, pos_stab2(1)], offsets(2)+[-pos_stab2(2),0, pos_stab2(2)], ...
    'Color', yellow,'LineStyle','-','LineWidth',3, 'DisplayName', 'Stable Manifold')

% Local unstable manifold
plot(offsets(1)+[-pos_ustab2(1),0, pos_ustab2(1)],offsets(2)+ [-pos_ustab2(2),0, pos_ustab2(2)], ...
    'Color', red,'LineStyle','-','LineWidth',3,'DisplayName', 'Unstable Manifold')

% Velocity impulse
plot(offsets(1)+[0, v_hist(1,1)],offsets(2)+[0, v_hist(2,1)],'Color',teal, ...
    'LineStyle','-.','LineWidth',2, 'DisplayName', '$\Delta v, R = 10^{-?}$');

plot(x(:,1), x(:,2),'Color',blue,'linestyle',':','LineWidth', 2, 'HandleVisibility','off')
hold off;
legend('show', 'FontSize', fsize, 'Interpreter', 'latex','location','southeast');
ylabel('y [LU]','Interpreter','latex', 'FontSize', fsize)
xlabel('x [LU]','Interpreter','latex', 'FontSize', fsize)

set(gca, 'TickLabelInterpreter', 'Latex')
ax = gca;
ax.XAxis.FontSize = fsize;
ax.YAxis.FontSize = fsize;
box on
grid

xlim(ax,[0.35,0.5])
ylim(ax,[-0.3, 0.3])

%% 

v1 = [v_hist(1,1), v_hist(2,1)];
v2 = [pos_stab2(1),pos_stab2(2)];

% Display the angle in degrees
fprintf('The angle between the vectors is %.2f degrees\n', subspace(v1',v2')*180/pi);


del_v_cost = sum(abs(del_v(1:2,1:14)), 'All')  % Comment out as neccesarry



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


function C = Jacobi(x)
    mu = 1.215058560962404e-2;
    d = sqrt((x(:,1) + mu).^2 + x(:,2).^2 + x(:,3).^2);
    r = sqrt((x(:,1) -1 + mu).^2 + x(:,2).^2 + x(:,3).^2);

    U = (x(:,1).^2 + x(:,2).^2)./2 + mu./r + (1-mu)./d;
    C = 2.*U - (x(:,4).^2 + x(:,5).^2 + x(:,6).^2);

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


%% Compute the state transition matrix

function [x,t,phi_T,PHI]=PHIgetHT(x0,tf)
% [x,t,phi_T,PHI]=PHIget(x0,tf)
% 
% Gets state transition matrix, phi_T, and the trajectory (x,t) for a length 
% of time, tf (2*pi is 1 period). In particular, for periodic solutions of 
% period tf=T, one can obtain the monodromy matrix, PHI(0,T).
%
% Function code from Shane Ross (revised 9.7.97) 
% Modified by Owen Brook 

OPTIONS = odeset('RelTol',1e-12,'AbsTol',1e-14*ones(1,6));

N=6;

PHI_0(1:N^2) 	   = reshape(eye(N),N^2,1);
PHI_0(1+N^2:N+N^2) = x0;

[t,PHI] = ode87(@(t,PHI) var3D(t,PHI), [0 tf], PHI_0, OPTIONS); 

x = PHI(:,1+N^2:N+N^2); % trajectory

phi_T = reshape(PHI(length(t),1:N^2),N,N); % monodromy matrix, PHI(O,T)

end



function [PHIdot] = var3D(t,PHI)
% function [PHIdot]=var3D(t,PHI);
%
% This here is a preliminary state transition, PHI(t,t0),
% matrix equation attempt for the CR3BP, based on...
%
%        d PHI(t, t0)
%        ------------ =  F(t) * PHI(t, t0)
%             dt
%
%  CONVENTION
%
%                 L4
%
%    L3-----M1-------L1---M2---L2         M1=1-mu, M2=mu
%
%                 L5
% Function code from Shane Ross (revised 7.1.97)

FORWARD = 1;
mu = 1.215058560962404e-2;
 
mu2 = 1-mu;
x(1:6) = PHI(37:42);
phi  = reshape(PHI(1:36), 6, 6);

r2= (x(1)+mu )^2 + x(2)^2 + x(3)^2;	% r: distance to m1, LARGER MASS
R2= (x(1)-mu2)^2 + x(2)^2 + x(3)^2;	% R: distance to m2, smaller mass
r3= r2^1.5; r5= r2^2.5;
R3= R2^1.5; R5= R2^2.5;

omgxx= 1+(mu2/r5)*(3*(x(1)+mu)^2)+(mu/R5)*(3*(x(1)-mu2)^2)-(mu2/r3+mu/R3);
omgyy= 1+(mu2/r5)*(3* x(2)^2    )+(mu/R5)*(3* x(2)^2     )-(mu2/r3+mu/R3);
omgzz=   (mu2/r5)*(3* x(3)^2    )+(mu/R5)*(3* x(3)^2     )-(mu2/r3+mu/R3);

omgxy= 3*x(2)*     (mu2*(x(1)+mu)/r5+mu*(x(1)-mu2)/R5); 
omgxz= 3*x(3)*     (mu2*(x(1)+mu)/r5+mu*(x(1)-mu2)/R5); 
omgyz= 3*x(2)*x(3)*(mu2          /r5+mu           /R5);

% Df is the jacobian matrix of the vector field f
   Df     =[   0     0     0     1     0	 0 ; 
	           0     0     0     0 	   1 	 0 ; 
        	   0	 0     0     0     0     1 ;
		    omgxx omgxy omgxz    0     2 	 0 ; 
         	omgxy omgyy omgyz   -2     0 	 0 ;
		    omgxz omgyz omgzz    0	   0	 0 ];

phidot = Df * phi;

PHIdot        = zeros(42,1);
PHIdot(1:36)  = reshape(phidot, 36, 1);
PHIdot(37)    = x(4);
PHIdot(38)    = x(5);
PHIdot(39)    = x(6);
PHIdot(40)    = x(1)-(mu2*(x(1)+mu)/r3) -(mu*(x(1)-mu2)/R3) + 2*x(5);
PHIdot(41)    = x(2)-(mu2* x(2)    /r3) -(mu* x(2)     /R3) - 2*x(4);
PHIdot(42)    =     -(mu2* x(3)    /r3) -(mu* x(3)     /R3);
PHIdot(37:42) = PHIdot(37:42)*FORWARD;

end