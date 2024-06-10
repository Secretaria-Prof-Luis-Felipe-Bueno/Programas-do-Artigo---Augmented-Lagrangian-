function [x,t,k] = secant_fuel_v1(P,a,b,c,l,u)
% secant method:
% g -> gamma_j \in U(1,5)
% l -> left side box constraint
% u -> right side box constraint
% c -> right side of knapsack constraint
%   Detailed explanation goes here
tic;
eps1 = 10^(-4);
[A,B] = busca_intervalo_fuel(P,a,b,c,l,u);
%[A,B] = refina_intervalo_fixo(A,B,a,c,l,u);
t = inf;
L0 = A/2;
L1 = B/2;
for k=1:10^3
    xL0 = xis_fuel_v1(L0,P,a,b,l,u);
    xL1 = xis_fuel_v1(L1,P,a,b,l,u);
    gL0 = b'*(xL0) - c;
    gL1 = b'*(xL1) - c;
    L2 = L1 - ((L1 - L0)/(gL1 - gL0))*gL1;
    xL2 = xis_fuel_v1(L2,P,a,b,l,u);
    gL2 = b'*(xL2) - c;
    if abs(gL2) < eps1
        t = toc();
        break;
    end
    L0 = L1;
    L1 = L2;
end
x = xis_fuel_v1(L2,P,a,b,l,u);
end