clc
close all
clear all

%%%%%-------------------------------------------------%%%%%
%%%%%-------------------------------------------------%%%%%
% This program deals with the quasi-one-dimensional
% subsonic-supersonic isentropic
% flow through a convergent-divergent nozzle.

% Based on the report titled "quasi 1D flow - without shock"
% The finite difference expression is set up using
% MacCormack's explicit technique for the numerical solution
% of governing equations given by Eqs (46),(48) & (50)


%%%%%--------------------------------------------------%%%%%
%%%%%--------------------------------------------------%%%%%

% Note : All parameters are nondimensionalized.


dx = 0.1;                               % grid size
c  = 0.5;                               % courant number
gamma = 1.4;                            % Ratio of specific heats

% Discretization of nondimensional distance along the nozzle
x = 0:dx:3.0;
k = length(x);
n = 25;

% Parabolic area distrubition along the nozzle
A = 1+2.2.*(x-1.5).^2;                         % Area,A/A*


% Grid distribution in the nozzle

for i = 1:k
    y = linspace(-A(i)/2,A(i)/2,n);
    for j = 1:n
        xx(i,j) = x(i);
        yy(i,j) = y(j);
    end
end

figure(1)
plot(xx,yy,'-r',xx',yy','-r')

% % Initialization

drhodt_av = zeros(1,k);
dVdt_av = zeros(1,k);
dTdt_av = zeros(1,k);
dta = zeros(1,k);

% Initial Condition
rho = 1-0.3146.*x;                      % density,rho/rho0, Eq (66a)
T = 1-0.2314.*x;                        % Temperature,T/T0, Eq (66b)
V = (0.1+1.09.*x).*T.^0.5;              % Velocity,V/a0,    Eq (66c)
P = rho.*T;                             % Pressure, P/P0
M = V.*T.^-0.5;                         % Mach number,M = V/a0/a/a0
mf = rho.*A.*V;                         % mass flow rate

mf_0 = mf;                              % mass flow rate @ t = 0
gamma1 = 1/gamma;
gamma2 = gamma-1;

resmax = 10^-8;
res = 1;
t = 0;
nstep = 0;

while res > resmax
    prev1 = drhodt_av;
    prev2 = dVdt_av;
    prev3 = dTdt_av;

    for i = 1:k
        dta(i) = (c*dx)/(T(i)^0.5+V(i));
    end
    dt = min(dta);

    % Predictor Step
    for i = 2:k-1;
        V1 = (V(i+1)-V(i))/dx;
        lnA1 = (log(A(i+1))-log(A(i)))/dx;
        rho1 = (rho(i+1)-rho(i))/dx;
        T1 = (T(i+1)-T(i))/dx;
        drhodt_1(i) = -rho(i)*V1-rho(i)*V(i)*lnA1-V(i)*rho1;                   % eq. (43)
        dVdt_1(i) = -V(i)*V1-gamma1*T1-gamma1*rho1*T(i)/rho(i);                % eq. (44)
        dTdt_1(i) = -V(i)*T1-gamma2*T(i)*V1-gamma2*T(i)*V(i)*lnA1;             % eq. (45)
    end

    % predicted values
    for i = 2:k-1
        rho_p(i) = rho(i)+drhodt_1(i)*dt;                                      % eq. (46)
        V_p(i) = V(i)+dVdt_1(i)*dt;                                            % eq. (47)
        T_p(i) = T(i)+dTdt_1(i)*dt;                                            % eq. (48)
    end

    t = t+dt;
    nstep = nstep +1;

    rho_p(1) = rho(1);
    V_p(1) = V(1);
    T_p(1) = T(1);

    % corrector step
    for i = 2:k-1;
        V1 = (V_p(i)-V_p(i-1))/dx;
        lnA1 = (log(A(i))-log(A(i-1)))/dx;
        rho1 = (rho_p(i)-rho_p(i-1))/dx;
        T1 = (T_p(i)-T_p(i-1))/dx;
        drhodt_2(i) = -rho_p(i)*V1-rho_p(i)*V_p(i)*lnA1-V_p(i)*rho1;        % eq. (49)
        dVdt_2(i) = -V_p(i)*V1-gamma1*T1-gamma1*rho1*T_p(i)/rho_p(i);       % eq. (50)
        dTdt_2(i) = -V_p(i)*T1-gamma2*T_p(i)*V1-gamma2*T_p(i)*V_p(i)*lnA1;  % eq. (51)
    end

    % Average time derivatives

    for i = 2:k-1;
        drhodt_av(i) = 0.5*(drhodt_1(i)+drhodt_2(i));                       % eq. (52)
        dVdt_av(i) = 0.5*(dVdt_1(i)+dVdt_2(i));                             % eq. (53)
        dTdt_av(i) = 0.5*(dTdt_1(i)+dTdt_2(i));                             % eq. (54)
    end

    % corrected values of flow field variables at time t+dt

    for i = 2:k-1
        rho(i) = rho(i)+drhodt_av(i)*dt;                                     % eq. (55)
        V(i) = V(i)+dVdt_av(i)*dt;                                           % eq. (56)
        T(i) = T(i)+dTdt_av(i)*dt;                                           % eq. (57)
        P(i) = rho(i)*T(i);
        M(i) = V(i)*T(i)^-0.5;
    end

    % Boundary condition @ first node
    V(1) = 2*V(2)-V(3);
    rho(1) = 1;
    T(1) = 1;
    P(1) = 1;
    M(1) = V(1)*T(1)^-0.5;

    % Boundary condition @ last node
    V(end) = 2*V(end-1)-V(end-2);
    rho(end) = 2*rho(end-1)-rho(end-2);
    T(end) = 2*T(end-1)-T(end-2);
    P(end) = rho(end)*T(end);
    M(end) = V(end)*T(end)^-0.5;

    % mass flow rate after t+dt
    mf = rho.*A.*V;

    if nstep == 50
        mf_50 = mf;
    elseif nstep == 100
        mf_100 = mf;

    elseif nstep == 150
        mf_150 = mf;
    elseif nstep == 200
        mf_200 = mf;
    elseif nstep == 700
        mf_700 = mf;
    end

    res1 = abs(drhodt_av(16));
    res2 = abs(dVdt_av(16));

    if res1 > res2
        res = res1;
    else
        res = res2;
    end

    rho_16(nstep) = rho(16);
    T_16(nstep) = T(16);
    P_16(nstep) = P(16);
    M_16(nstep) = M(16);

end

for i = 1:k
    for j = 1:n
        rho1(i,j) = rho(i);
        T1(i,j) = T(i);
        P1(i,j) = P(i);
        M1(i,j) = M(i);
        V1(i,j) = V(i);
    end
end

figure (2)
plot(1:nstep,rho_16(1:nstep),'-r','LineWidth',2);
xlabel('Number of time steps','Fontsize',12)
h = ylabel('$\frac{\rho}{\rho_{0}}$','Interpreter','latex');
h.FontSize=20;
set(h,'position',get(h,'position')+[-50 0 0]);
set(h,'rotation',0);

figure (3)
plot(1:nstep,T_16(1:nstep),'-k','LineWidth',2);
xlabel('Number of time steps','Fontsize',12)
h = ylabel('$\frac{T}{T_{0}}$','Interpreter','latex');
h.FontSize=20;
set(h,'position',get(h,'position')+[-50 0 0]);
set(h,'rotation',0)


figure (4)
plot(1:nstep,P_16(1:nstep),'-b','LineWidth',2);
xlabel('Number of time steps','Fontsize',12)
h = ylabel('$\frac{p}{p_{0}}$','Interpreter','latex');
h.FontSize=20;
set(h,'position',get(h,'position')+[-50 0 0]);
set(h,'rotation',0);

figure (5)
plot(1:nstep,M_16(1:nstep),'-g','LineWidth',2);
xlabel('Number of time steps','Fontsize',12)
h = ylabel('$M$','Interpreter','latex');
h.FontSize=20;
set(h,'position',get(h,'position')+[-50 0 0]);
set(h,'rotation',0);


% creating file and writing steady state values on it

fid = fopen( 'flow field variables', 'wt' );
fprintf(fid,'x/L     A/A*   rho/rho0   V/a0   T/T0   p/p0     M   mass flow\n')

%fprintf(fid,'%s'/t,'');
for j = 1:k

    B = [x(j);A(j);rho(j);V(j);T(j);P(j);M(j);mf(j)];
    fprintf(fid,'%.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n',B);

end

figure(6)
plot(x,rho,'-k','LineWidth',2)
hold on
plot(x,P,'-g','LineWidth',2)
hold on
plot(x,T,'-b','LineWidth',2)
legend('\rho/\rho_{0}','p/p_{0}','T/T_{0}')
xlabel('Nondimensional distance through nozzle (x)','FontSize',12)
ylabel('Density ratio,Pressure ratio,Temperature ratio')

figure(7)
plot(x,mf_0,'--r',x,mf_50,'r',x,mf_100,'r',x,mf_150,'r',x,mf_200,'r',x,mf_700,'r')
xt = [2.4 2.4 2.4 2.4 2.4 2.4 ];
i = find(x==2.4);
yt = [mf_0(i)+0.08 mf_50(i)+0.08 mf_100(i)+0.08 mf_150(i)+0.08 mf_200(i)+0.08 mf_700(i)+0.08];
str = {'0\Deltat','50\Deltat','100\Deltat','150\Deltat','200\Deltat','700\Deltat'};
text(xt,yt,str)
grid on
xlabel('Nondimensional distance through nozzle (x)','Fontsize',12)
h = ylabel('Nondimensional mass flow (\rhoVA)/(\rho_{0}a_{0}A^{*})');
h.FontSize=12;

figure(8)
surf(xx,yy,P1);view(2)
shading interp
colormap(jet(256))
colorbar
title('Presure distrubtion across the nozzle')

