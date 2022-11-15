clc
clear all
close all
%%%%%-------------------------------------------------%%%%%
%%%%%-------------------------------------------------%%%%%

% This program deals with the quasi-one-dimensional
% flow through a convergent-divergent nozzle using shock 
% capturing method.
% The finite difference expression is set up using
% MacCormack's explicit technique for the numerical solution
% of governing equations in conservation form.

%%%%%--------------------------------------------------%%%%%
%%%%%--------------------------------------------------%%%%%

% Note : All parameters are nondimensionalized.

dx = 0.05;             % grid size                   
c  = 0.5;              % courant number                 
gamma = 1.4;           % Ratio of specific heats

% Discretization of nondimensional distance along the nozzle
x = 0:dx:3.0;
k = length(x);
n = 25;

Cx = 0.2;               % artificial viscosity arbitrary parameter

% Parabolic area distrubition along the nozzle
A = 1+2.2.*(x-1.5).^2;   % Area,A/A*                    

% Grid distribution in the nozzle
for i = 1:k
    y = linspace(-A(i)/2,A(i)/2,n);
    for j = 1:n
        xx(i,j) = x(i);
        yy(i,j) = y(j);
    end
end

% Initial Condition 

for i = 1:k
    if x(i) >= 0 && x(i) < 0.5
        rho(i) = 1;                          % Eq. (58a)      
        T(i) = 1;                            % Eq. (58b)
    elseif x(i) >= 0.5 && x(i) < 1.5
        rho(i) = 1-0.366*(x(i)-0.5);         % Eq. (58c)
        T(i) = 1-0.167*(x(i)-0.5);           % Eq. (58d)
    elseif x(i) >= 1.5 && x(i) < 2.1
        rho(i) = 0.634-0.702*(x(i)-1.5);     % Eq. (58e)
        T(i) = 0.833-0.4908*(x(i)-1.5);      % Eq. (58f)
    else
        rho(i) = 0.5892+0.10228*(x(i)-2.1);  % Eq. (58g)
        T(i) = 0.93968+0.0622*(x(i)-2.1);    % Eq. (58h)
    end
end

V = 0.59./(rho.*A);                          % Velocity,V/a0,[Eq.(59)]
P = rho.*T;                                  % pressure,p/p0
mf_0 = rho.*A.*V;                            % mass flow rate @ t = 0

gamma1 = 1/gamma;
gamma2 = gamma-1;

% Initial condition of solution vectors 
U1 = rho.*A;                                           
U2 = rho.*A.*V;                                        
U3 = rho.*A.*((T./gamma2)+(0.5*gamma.*V.^2));          

% Initial condition of flux vectors 
F1 = U2;                                                    % Eq.(23)
F2 = (U2.^2./U1)+(1-gamma1).*(U3-0.5*gamma.*U2.^2./U1);     % Eq.(25)
F3 = (gamma.*U2.*U3./U1)-0.5*gamma*gamma2.*U2.^3./U1.^2;    % Eq.(27)

resmax = 10^-6;      % maximum error
res = 1;
t = 0;       
nstep = 0;          

while res > resmax
    dta = (c*dx)./(V+T.^0.5);
    dt = min(dta);
    
    nstep = nstep+1;
     
    % Predictor Step
     for i = 2:k-1
         J2_p = gamma1*rho(i)*T(i)*((A(i+1)-A(i))/dx);
         dU1dt_p(i) = -(F1(i+1)-F1(i))/dx;
         dU2dt_p(i) = -(F2(i+1)-F2(i))/dx + J2_p;
         dU3dt_p(i) = -(F3(i+1)-F3(i))/dx;
     end
     
     for i = 2:k-1
         num = abs(P(i+1)-2*P(i)+P(i-1));
         den = P(i+1)+2*P(i)+P(i-1);
         S1_p(i) = (Cx*num*(U1(i+1)-2*U1(i)+U1(i-1)))/den;
         S2_p(i) = (Cx*num*(U2(i+1)-2*U2(i)+U2(i-1)))/den;
         S3_p(i) = (Cx*num*(U3(i+1)-2*U3(i)+U3(i-1)))/den;
     end

     for i = 2:k-1
         U1_p(i) = U1(i)+(dU1dt_p(i)*dt)+S1_p(i);      % Eq.(61)
         U2_p(i) = U2(i)+(dU2dt_p(i)*dt)+S2_p(i);      % Eq.(62)
         U3_p(i) = U3(i)+(dU3dt_p(i)*dt)+S3_p(i);      % Eq.(63)
     end
     
     U1_p(1) = U1(1);
     U2_p(1) = U2(1);
     U3_p(1) = U3(1);
     U1_p(k) = U1(k);
     U2_p(k) = U2(k);
     U3_p(k) = U3(k);

     
for i = 1:k
     rho_p(i) = U1_p(i)/A(i);                               % Eq.(17)
     V_p(i) = U2_p(i)/U1_p(i);                              % Eq.(18)
     T_p(i) = gamma2*((U3_p(i)/U1_p(i))-0.5*gamma*V_p(i)^2);% Eq.(19)
     P_p(i) = rho_p(i)*T_p(i);                              % Eq.(20)
     F1_p(i) = U2_p(i);                                     % Eq.(23)
     F2_p(i) = (U2_p(i)^2/U1_p(i))+(1-gamma1)*...
               (U3_p(i)-0.5*gamma*U2_p(i)^2/U1_p(i));       % Eq.(25)
     F3_p(i) = (gamma*U2_p(i)*U3_p(i)/U1_p(i))-...
               (0.5*gamma*gamma2*U2_p(i)^3/U1_p(i)^2);      % Eq.(27)
end
     

     % corrector Step
     for i = 2:k-1
         J2_c = gamma1*rho_p(i)*T_p(i)*((A(i)-A(i-1))/dx);
         dU1dt_c(i) = -(F1_p(i)-F1_p(i-1))/dx;
         dU2dt_c(i) = -(F2_p(i)-F2_p(i-1))/dx + J2_c;
         dU3dt_c(i) = -(F3_p(i)-F3_p(i-1))/dx;
     end
     
     for i = 2:k-1
         num = abs(P_p(i+1)-2*P_p(i)+P_p(i-1));
         den = P_p(i+1)+2*P_p(i)+P_p(i-1);
         S1_c(i) = (Cx*num*(U1_p(i+1)-2*U1_p(i)+U1_p(i-1)))/den;
         S2_c(i) = (Cx*num*(U2_p(i+1)-2*U2_p(i)+U2_p(i-1)))/den;
         S3_c(i) = (Cx*num*(U3_p(i+1)-2*U3_p(i)+U3_p(i-1)))/den;
     end
     
      % Average time derivatives
      
      for i = 2:k-1
          dU1dt_av(i) = 0.5*(dU1dt_p(i)+ dU1dt_c(i));
          dU2dt_av(i) = 0.5*(dU2dt_p(i)+ dU2dt_c(i));
          dU3dt_av(i) = 0.5*(dU3dt_p(i)+ dU3dt_c(i));
      end
      
      for i = 2:k-1
          U1(i) = U1(i)+ (dU1dt_av(i)*dt)+S1_c(i);  % Eq.(64)
          U2(i) = U2(i)+ (dU2dt_av(i)*dt)+S2_c(i);  % Eq.(65)
          U3(i) = U3(i)+ (dU3dt_av(i)*dt)+S3_c(i);  % Eq.(66) 
      end
      
  % Boundary condition @ first node
      U1(1) = A(1);                    % Eq. (47)
      U2(1) = 2*U2(2)-U2(3);           % Eq. (48)
      ve = U2(1)/U1(1);
      U3(1) = U1(1)*((1/gamma2)+0.5*gamma*ve^2); % Eq.(50)
     
  % Boundary condition @ last node
    U1(k) = 2*U1(k-1)-U1(k-2);           % Eq. (51a)
    U2(k) = 2*U2(k-1)-U2(k-2);           % Eq. (51b)
    U3(k) = (0.6784*A(k)/gamma2)+...     % Eq. (55)
            (0.5*gamma*U2(k)^2/U1(k));
    
    for i = 1:k
        F1(i) = U2(i);
        F2(i) = (U2(i)^2/U1(i))+(1-gamma1)*(U3(i)-0.5*gamma...
                *U2(i)^2/U1(i));
        F3(i) = (gamma*U2(i)*U3(i)/U1(i))-0.5*gamma*gamma2...
                *U2(i)^3/U1(i)^2;
    end
    
    
    % corrected values of primitive variables
    for i = 1:k
    rho(i) = U1(i)/A(i);                            % Eq.(17) 
    V(i) = U2(i)/U1(i);                             % Eq.(18)
    T(i) = gamma2*((U3(i)/U1(i))-0.5*gamma*V(i)^2); % Eq.(19)
    M(i) = V(i)*T(i)^-0.5;
    P(i) = rho(i)*T(i);                             % Eq.(20)
    end

    % mass flow rate after t+dt
    mf = rho.*A.*V;
    
    res1 = abs(dU1dt_av(16));
    res2 = abs(dU2dt_av(16));
    
    if res1 > res2
        res = res1;
    else
        res = res2;
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

figure(1)
if mod(nstep,5) == 0
surf(xx,yy,M1);view(2)
shading interp
colormap(jet(256))
colorbar
axis off
drawnow
title('Mach number contour')
end
      
if nstep == 2600,break,end
end
figure(2)
plot(x,P,'-r','Linewidth',2)
xlabel('x/L')
ylabel('Nondimensional pressure p/p_{0}')

 figure(3)
 plot(x,M,'-b','Linewidth',2)
 xlabel('x/L')
ylabel('Mach number')

 figure(4)
 plot(x,mf,'-k','Linewidth',2)
 xlabel('x/L')
ylabel('\rhoAV/\rho_{0}A*a_{0}')
axis([0 3.0 0.56 0.68])
