% Finds a controller by solving LMIs and semidefinite program
% Requires robust control toolbox

% Only velocity perts
R = zeros([6,3]);
R(4,1) = 1;
R(5,2) = 1;
R(6,3) = 1;
B = A * R;

% Specify matrix variables
setlmis([])
Q = lmivar(1,[6,1]); % Type 1 - Symmetric 6x6 matrix
Y = lmivar(2,[3,6]); % Type 2 - 3x6 matrix


% Specifying LMI 1 - Big matrix
LMI_1 = newlmi;

lmiterm([-LMI_1, 1, 1, Q], 1, 1); % Q - pos(1,1), term 1, factor=1, outer factor=1
lmiterm([-LMI_1, 1, 2, Q], 1, A');
lmiterm([-LMI_1, 1, 2, -Y], 1, B');
lmiterm([-LMI_1, 2, 2, Q], 1, 1);

% Specifying LMI 2 - Q > 0 OR 0 < Q
LMI_2 = newlmi;
lmiterm([-LMI_2, 1, 1, Q], 1, 1);
LMISYS = getlmis;

fR = 1e-10;  % Feasibility radius 
options = [0,0,fR,0,0];
[tmin, xfeas] = feasp(LMISYS,options,0);

Q_sol = dec2mat(LMISYS,xfeas,Q);
Y_sol = dec2mat(LMISYS,xfeas,Y);

K_sat = Y_sol * inv(Q_sol)

C_sat = A + B*K_sat;
[~,sat_eig] = eig(C_sat);
abs(sat_eig)
