%os dados precisam ja estar salvos... naturalmente...
clear all
clc

n=100; m=100;

times=333*ones(30,1); %se der 333 eh pq ele n salvou por algum motivo estranho
iters=6002*ones(30,1); %iter max deveria ser 6000, se der 6002...
norms = 333*ones(30,1);
conds = zeros(30,1); lambs = -ones(30,1); 
erres=conds; epslons = conds;

timesprog=zeros(30,1);

i= 5;
%for i=1:30
mystr = strcat('.\Salvar\matrizes',num2str(n),num2str(m),num2str(i));
                
dadfile = matfile(mystr);
P = dadfile.P;
a = dadfile.a;
b = ones(n,1);
c = n;
l = zeros(n,1); %-(n*ones(n,1))/2; %
u = n*ones(n,1);  %(n*ones(n,1))/2;   
%lambda = -37.0092;
lambda=10;

[k,timen,xk0,lamb,r,epslon,xquadprog,timeprog] = func_teste_nda_aleatorio(6,P,a,b,c,l,u,lambda,1e-4,1);
times(i) = timen;
timesprog(i) = timeprog;
norms(i) = norm(xk0-xquadprog);
iters(i) = k;
conds(i) = cond(P);
lambs(i) = lamb;
erres(i) = r;
epslons(i) = epslon;

%disp('Problema '); disp(num2str(i));
%pause

%end  %end for do 30
