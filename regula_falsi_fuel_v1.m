function [x,t,k] = regula_falsi_fuel_v1(P,a,b,c,l,u,altlam)
%Bisection method:
% a -> gamma_j \in U(1,5)
% l -> left side box constraint
% u -> right side box constraint
% c -> right side of knapsack constraint
%   Detailed explanation goes here
tic;
eps1 = 10^(-4);
t = inf;
[A,B] = busca_intervalo_fuel(P,a,b,c,l,u);
for k=1:10^3
    xA = xis_fuel_v1(A,P,a,b,l,u);
    xB = xis_fuel_v1(B,P,a,b,l,u);
    gA = b'*(xA) - c;
    gB = b'*(xB) - c;
    L = (A*gB - B*gA)/(gB - gA);
    if (nargin==7)&&(k==1)
     xL = xis_fuel_v1(altlam,P,a,b,l,u);   
    else
     xL = xis_fuel_v1(L,P,a,b,l,u);
    end
    gL = b'*(xL) - c;
    if abs(gL) < eps1
        t = toc();
        break;
    end
    if sign(gL) == sign(gA)
        A = L;
    else
        B = L;
    end
end
x = xis_fuel_v1(L,P,a,b,l,u);
end