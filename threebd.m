clear all
close all
clc

%% General Data
m_Earth = 6.0477e24;        % kg
m_Sun = 1.9891e30;          % kg

mass = m_Earth + m_Sun;     % kg
time = 5022415;             % s
dist = 149597870700;        % m
mu = m_Earth/mass;


%% Lagrangian Points
syms x
eqn1 = x - (1-mu)/abs((mu+x))^3*(mu+x) + mu/abs((1-mu-x))^3*(1-mu-x) == 0;
L1x = vpasolve(eqn1, x,[eps;1-mu]);
L2x = vpasolve(eqn1, x,[1-mu;Inf]);
L3x = vpasolve(eqn1, x,[-Inf;-eps]);

L1 = double([L1x;0;0])
L2 = double([L2x;0;0])
L3 = double([L3x;0;0])

%% Shift Lagrangian Points
L1bis = L1-0.011*[1;0;0];
L3bis = L3+0.011*[1;0;0];

[nL1bis, betaL1Bis] = shiftLagrangian(L1bis)
[nL3bis, betaL3Bis] = shiftLagrangian(L3bis)
[nVerify, betaVerify] = shiftLagrangian([0.95;0;0.1])

sideL1bis = sqrt(solarSailArea(betaL1Bis, L1bis, mu))
sideL3bis = sqrt(solarSailArea(betaL3Bis, L3bis, mu))

%% Lagrangian Points Stability

AL1 = jacobianLagrangian(L1, mu);
AL3 = jacobianLagrangian(L3, mu);
format short
[vectL1, valL1] = eig(AL1)
[vectL3, valL3] = eig(AL3)

%% Manifolds
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
tspanF = [0 5*365.25*24*3600/time];
tspanB = [0 -5*365.25*24*3600/time];
pert = 1e-5;

% L1
x0L1S  = [L1;0;0;0] + real(pert*[vectL1(1:2,2);0;vectL1(4:5,2);0]);
x0L1S2 = [L1;0;0;0] - real(pert*[vectL1(1:2,2);0;vectL1(4:5,2);0]);
x0L1U  = [L1;0;0;0] + real(pert*[vectL1(1:2,1);0;vectL1(4:5,1);0]);
x0L1U2 = [L1;0;0;0] - real(pert*[vectL1(1:2,1);0;vectL1(4:5,1);0]);
[tL1S,xL1S] = ode45(@threebd_ode,tspanB,x0L1S,options);
[tL1S2,xL1S2] = ode45(@threebd_ode,tspanB,x0L1S2,options);
[tL1U,xL1U] = ode45(@threebd_ode,tspanF,x0L1U,options);
[tL1U2,xL1U2] = ode45(@threebd_ode,tspanF,x0L1U2,options);

figure
plot(1-mu, 0, 'bo', 'LineWidth', 2)
hold on
plot(L1(1), L1(2), 'bx', 'LineWidth', 2)
plot(xL1S(:,1), xL1S(:,2), 'k', 'LineWidth', 2)
plot(xL1S2(:,1), xL1S2(:,2), 'k--', 'LineWidth', 2)
plot(xL1U(:,1), xL1U(:,2), 'r', 'LineWidth', 2)
plot(xL1U2(:,1), xL1U2(:,2), 'r--', 'LineWidth', 2)
plot(L2(1), L2(2), 'bx', 'LineWidth', 2)
plot(L3(1), L3(2), 'bx', 'LineWidth', 2)
hold off
axis equal
set(gca, 'FontSize', 14)
title('$L_1$ Manifolds ($\epsilon = 10^{-5}$, 5 years)','interpreter','latex')
title('$L_1$ Manifolds Detailed','interpreter','latex')
xlabel('X [-]');
ylabel('Y [-]');
% legend({'Earth','Lagrangian Point','Stable Pos.', 'Stable Neg.', 'Unstable Pos.', 'Unstable Neg.'},'Location','best')


% L3
% x0L3S  = [L3;0;0;0] + real(pert*[vectL3(1:2,3);0;vectL3(3:4,3);0]);
% x0L3S2 = [L3;0;0;0] - real(pert*[vectL3(1:2,3);0;vectL3(3:4,3);0]);
% x0L3U  = [L3;0;0;0] + real(pert*[vectL3(1:2,4);0;vectL3(3:4,4);0]);
% x0L3U2 = [L3;0;0;0] - real(pert*[vectL3(1:2,4);0;vectL3(3:4,4);0]);
% [tL3S,xL3S] = ode45(@threebd_ode,tspanB,x0L3S,options);
% [tL3S2,xL3S2] = ode45(@threebd_ode,tspanB,x0L3S2,options);
% [tL3U,xL3U] = ode45(@threebd_ode,tspanF,x0L3U,options);
% [tL3U2,xL3U2] = ode45(@threebd_ode,tspanF,x0L3U2,options);
% 
% figure
% plot(L3(1), L3(2), 'bx', 'LineWidth', 2)
% hold on
% plot(xL3S(:,1), xL3S(:,2), 'k', 'LineWidth', 2)
% plot(xL3S2(:,1), xL3S2(:,2), 'k', 'LineWidth', 2)
% plot(xL3U(:,1), xL3U(:,2), 'r', 'LineWidth', 2)
% plot(xL3U2(:,1), xL3U2(:,2), 'r', 'LineWidth', 2)
% 
% axis equal
% hold off

%% Solar Sail Manifolds
time = 5022415;             % s
tspanF = [0 5*365.25*24*3600/time];
tspanB = [0 -5*365.25*24*3600/time];
pert = 1e-5;
x0U = [L1;0;0;0] + real(pert*vectL1(:,1));
x0S = [L3;0;0;0] + real(pert*vectL3(:,3));

i = 1;
eucl  = zeros(1,21);
pos   = zeros(1,57);
vel   = zeros(1,57);
timeU = zeros(1,57);
timeS = zeros(1,57);

% For the ANN
nData = 20000;
data  = zeros(57*nData,4);
genData = 0;
for alpha = -70:2.5:70
    [tU,xU,tS,xS] = solarSailTransfer(alpha*pi/180, x0U, x0S, tspanF, tspanB);
    
    [small, idxU] = pdist2(xU,xS, 'euclidean', 'Smallest', 1);
    
    if genData
        euclData = pdist2(xU,xS, 'euclidean');
        [array,order] = sort(euclData(:));
        for k = 1:nData
            [idxUData,idxSData] = ind2sub(size(euclData),order(k*1000));
            data_k = [alpha*pi/180,tU(idxUData),tS(idxSData),array(k*1000)];
            data(nData*(i-1)+k,:) = data_k;
        end
    end
    
    % Select minimum for each Alpha
    [small, idxS] = min(small);
    idxU = idxU(idxS);
    
    eucl(i)  = small;
    pos(i)   = norm(xU(idxU,1:3)-xS(idxS,1:3))*dist/1000;
    vel(i)   = norm(xU(idxU,4:6)-xS(idxS,4:6))*dist/time/1000;
    timeU(i) = tU(idxU)*time/(365.25*24*3600);
    timeS(i) = tS(idxS)*time/(365.25*24*3600);
    i = i+1;
end

if genData
    save('ANNdata.txt', 'data', '-ascii', '-double', '-tabs');
end

% Plot Cone Angle vs. Euclidean Norm
alpha = -70:2.5:70;
figure
plot(alpha, eucl, 'LineWidth', 2)
title('Minimum Euclidean Norm vs. Cone Angle (\alpha)')
xlabel('\alpha (^\circ)')
ylabel('Euclidean Norm (-)')
set(gca, 'FontSize', 14)

% Select Overall Best
[euclBest, idx] = min(eucl);
alphaBest = alpha(idx);
posBest   = pos(idx);
velBest   = vel(idx);
timeUBest = timeU(idx);
timeSBest = timeS(idx);

% Compute and Plot Overall Best Manifolds
tspanF = [0 timeUBest*365.25*24*3600/time];
tspanB = [0 timeSBest*365.25*24*3600/time];
[tU,xU,tS,xS] = solarSailTransfer(alphaBest*pi/180, x0U, x0S, tspanF, tspanB);

% In 3D
figure
plot3(L1(1), L1(2), L1(3), 'bo', 'LineWidth', 2)
hold on
plot3(L3(1), L3(2), L3(3), 'bo', 'LineWidth', 2)
plot3(xU(:,1), xU(:,2), xU(:,3), 'r', 'LineWidth', 2)
plot3(xS(:,1), xS(:,2), xS(:,3), 'k', 'LineWidth', 2)
axis equal
hold off

% In 2D
figure
plot(L1(1), L1(2), 'bo', 'LineWidth', 2)
hold on
plot(L3(1), L3(2), 'bx', 'LineWidth', 2)
plot(xU(:,1), xU(:,2), 'r', 'LineWidth', 2)
plot(xS(:,1), xS(:,2), 'k', 'LineWidth', 2)
hold off
axis equal
set(gca, 'FontSize', 14)
title('$L_1-L_3$ Heteroclinic Connection ($\epsilon = 10^{-5}$, $\beta = 0.025$, $\alpha  = -22.5^\circ$)','interpreter','latex')
xlabel('X [-]');
ylabel('Y [-]');
legend({'L_1','L_3','Unstable Leg', 'Stable Leg'},'Location','best')



function [tU,xU,tS,xS] = solarSailTransfer(alpha, x0U, x0S, tspanF, tspanB)

options = odeset('RelTol',1e-12,'AbsTol',1e-12);

[tU,xU] = ode45(@solarsail_ode,tspanF,x0U,options);
[tS,xS] = ode45(@solarsail_ode,tspanB,x0S,options);

function Y = solarsail_ode(~,y)

Y = zeros(6,1);
m_Earth = 6.0477e24;        % kg
m_Sun = 1.9891e30;          % kg

mass = m_Earth + m_Sun;     % kg
time = 5022415;             % s
dist = 149597870700;        % m
mu = m_Earth/mass;

r  = y(1:3);
dr = y(4:6);

r1 = [y(1)+mu; y(2); y(3)];
r2 = [y(1)-(1-mu); y(2); y(3)];
w  = [0;0;1];
th  = cross(w,r1/norm(r1));
eta = cross(r1/norm(r1),th);

beta = 0.025;
delta = 0;
n = [r1/norm(r1) th eta]*[cos(alpha); sin(alpha)*sin(delta); sin(alpha)*cos(delta)];
a_sail = beta*(1-mu)/norm(r1)^2*(dot(r1/norm(r1),n)^2)*n;

Y(1:3) = dr;
Y(4:6) = a_sail - ((1-mu)/norm(r1)^3*r1 + mu/norm(r2)^3*r2) - 2*cross(w,dr) - cross(w,cross(w,r));
end

end



%% Auxiliary Functions
function [n, beta, a] = shiftLagrangian(p)
m_Earth = 6.04723e24;   % kg
m_Sun = 1.989e30;       % kg
mass = m_Earth + m_Sun; % kg
mu = m_Earth/mass;

x   = p(1); y = p(2); z = p(3);
r1  = p-[-mu;0;0];
r2  = p-[1-mu;0;0];
dUx = x - (1-mu)/norm(r1)^3*(mu+x) + mu/norm(r2)^3*(1-mu-x);
dUy = y - (1-mu)/norm(r1)^3*y - mu/norm(r2)^3*y;
dUz = -(1-mu)/norm(r1)^3*z - mu/norm(r2)^3*z;
dU  = [dUx;dUy;dUz];

n       = dU/norm(dU);
beta    = norm(r1)^2*dot(dU,n)/((1-mu)*dot(r1/norm(r1),n)^2);
a       = beta*((1-mu)/norm(r1)^2)*dot(r1/norm(r1),n)^2*n;
end

function area = solarSailArea(beta, p, mu)
P = 4.54e-6;            % N/m2
m = 10;                 % kg
muS = 1.327178e20;      % m3/s2
dist = 149597870700;    % m
R  = norm(dist*(p-[-mu;0;0]));

area = beta*muS*m/(2*P*R^2);
end

function A = jacobianLagrangian(p, mu)
    r1  = norm(p-[-mu;0;0]);
    r2  = norm(p-[1-mu;0;0]);
    K = (1 - mu)/r1^3 + mu/r2^3;

    Uxx = 1 + 2*K;
    Uyy = 1 - K ;
    Uzz = -K;

    A = [zeros(3), eye(3);
        Uxx, 0, 0, 0, 2, 0;
        0, Uyy, 0, -2, 0, 0;
        0, 0, Uzz, 0, 0, 0];

%     A = [zeros(2), eye(2);
%         Uxx, -Uxy, 0, 2;
%         -Uxy, Uyy, -2, 0];
end

function Y = threebd_ode(~,y)
Y = zeros(6,1);
m_Earth = 6.0477e24;        % kg
m_Sun = 1.9891e30;          % kg

mass = m_Earth + m_Sun;     % kg
time = 5022415;             % s
dist = 149597870700;        % m
mu = m_Earth/mass;

r  = y(1:3);
dr = y(4:6);
r1 = [y(1)+mu; y(2); y(3)];
r2 = [y(1)-(1-mu); y(2); y(3)];
w  = [0;0;1];

Y(1:3) = dr;
Y(4:6) = -((1-mu)/norm(r1)^3*r1 + mu/norm(r2)^3*r2) - 2*cross(w,dr) - cross(w,cross(w,r));

end