function soma = solut_fuel(x,P,a,c,l,u)
n = length(x);
df = P\a; % derivada de f
pj = x - df;
options = optimoptions(@quadprog,'Algorithm','interior-point-convex','StepTolerance',10^(-4));
projx = quadprog(speye(n),-pj,[],[],ones(1,n),c,l,u,[],options);
soma = norm(x - projx);
end