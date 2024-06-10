clear all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dados do problema  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%resolve o problema [x y z][1 -1 0; -1 2 0; 0 0 1][x y z] sa 
%z=0, [x y z] in [1 0 -1, 5 1 1]

prec = 1e-4;
vaiplot=1;
%p = ones(n,1);

%% gracao de dados
% p= matriz do problema]

problema=8;
switch problema
    case 1 
        P = [1, - 1, 0; -1, 2, 0; 0, 0, 1]; b=[0;0;1]; a=[0;0;0];
         l=[1;0;-1]; u=[5;1;1]; c=0; n=length(u);
    case 2 
        P = [1, - 2, 0; -2, 5, 0; 0, 0, 1]; b = [1;1;1]; a=[0;0;0];
         l=[1;0;-10^20]; u=[5;1;10^20]; c=0; n=length(u);
    case 3 
        P = [2,-1;-1,6]; b=[0;0]; a=[0;0];
        l=[1;0]; u=[5;1]; c=0; n=length(u);
    case 4
        P = [4,1, 1;1, 4, 1; 1, 1, 4]; b = [1;1;1]; a=[0;0;0];
         l=[1;0;-5]; u=[5;1;5]; c=0; n=length(u);
    case 5
        P = [4,1; 1, 4]; b = [1;1]; a=[0;0];
         l=[0;0]; u=[2;2]; c=0; n=length(u);
    case 6
        P = [3/4, -1;-1, 3/2]; b=[0;1]; a=[-1;2];
        l=[0;0]; u=[1;1]; c=1; n=length(u);
    case 7
        P = [3, 1;1, 3]; b=[1;1]; a=[4;-4];
        l=[0;-1]; u=[5;1]; c=0; n=length(u);
    case 8
        P = 10*[9/2, -2; -2,1]; b=[0;1]; a = 10*[-3;2];
        l=[0;0]; u=[1;1]; c=1; n=length(u);
        
end

%b = p;

% p=10 + (25-10).*rand(10,1) U[10,25];
%a = (2:n+1)';

% p=10 + (25-10).*rand(10,1) U[10,25];

btl=l'*b;
btu=u'*b;

disp('gerou dados')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dados inciais  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xk0=[0;1]; % Ponto inicial
lamb=0;        % Multiplicador inicial a restricao de igualdade b'*x-c
v=ones(n,1);   % multiplicador associado a restricao l-x
w=v;         % multiplicador associado a restricao x-u
r=1;          % Parametro de penalidade inicial
%r=100;
if vaiplot==1
   xplot(1) = xk0(1); yplot(1) = xk0(2);
   xpplot(1) = xk0(1); ypplot(1) = xk0(2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop Principal  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tt=cputime; k=1;  aux=0;  % Time
%calculando P inversa b:
pinvb = P\b; pinva = P\a; rest0 = b'*xk0 - c; rest=rest0;
while (k<6000)&&(aux==0)
  
% %   bt=sqrt(r)*b;
% %   alfa=1/(1+bt'*bt);
% %   IP= eye(n)-(alfa*(bt*bt')); % inversa de P formula de Sherman-Morrisson
% %   xk=IP*(a-((lamb-r*c)*b)); % Solucao de grad(lag)=0

alfa=r/(1+r*b'*pinvb);
lar = r*c - lamb;
xk = pinva + (lar - alfa*(b'*pinva + lar*(b'*pinvb)))*pinvb
pause
sexkoutra = quadprog(P,-2*xk ,[],[],[],[],l,u)

if vaiplot==1
   xplot(k+1) = xk(1); yplot(k+1) = xk(2);
end

xk=max(l,min(xk,u))
pause

if vaiplot==1
   xpplot(k+1) = xk(1); ypplot(k+1) = xk(2);
end

 %% Proje????o da solucao na caixa
% %     for i=1:n
% %         if xk(i) < l(i) 
% %             xk(i)=l(i);   % Componente na solu????o fora da caixa do lado esquerdo
% %         elseif xk(i) > u(i) 
% %             xk(i)=u(i);   % Componente na solu????o fora da caixa do lado direito
% % %         else
% % %             xk(i)=xk(i);  % Nao ?? necessario??
% %         end
% %     end
    
  % Atualizamos os multiplicadores das restri????es da caixa mul, muu
  %gradlagx2=xk-a+(lamb+r*(b'*xk-c))*b; %avaliamos o gradiente da lagrangeana em mul
   
%    for i=1:n
%         if xk(i) <= l(i) 
%           v(i)=xk(i)-a(i)+(lamb+r*(b'*xk-c))*b(i); % multiplicador v inviavel
%         else 
%             v(i)=0; % N??o precissa associado a restricao x-u
%         end
%     end
%    
%    for i=1:n
%         if xk(i) <= u(i) 
%           w(i)=0 ; % multiplicador associado a restricao x-u
%         %elseif l(i)<=xk(i) && xk(i)<=u(i) %????????????
%          %       w(i)=0;
%          else
%             w(i)=-(xk(i)-a(i)+(lamb+r*(b'*xk-c))*b(i)); % N??o precissa
%         end
%     end
   %% Atualizamos o multiplicador da restricao de igualdade
    %rest0=b'*xk0-c; %%
    rest=b'*xk-c;
    
    
    %if abs(rest)>0.1*abs(rest0)
    %    r=2*r;
    %end   
    rest0 = rest;
    
    difproj = xk - (P*xk + (b'*xk)*r*b - a - lar*b)
    difproj=max(l,min(difproj,u))
    erro3=norm(xk-difproj)
        
    %r=n/(4*k+1);x0
    
    %lamb0=lamb;
    lamb=lamb+r*rest
    %r=abs(lamb-lamb0)/abs(rest0)
     % k=k+1;
      
      %% Criterio de parada
    erro1=norm(xk-xk0);
    erro2=abs(rest);
    %pause
    
    if erro1 <= prec && erro2 <= prec
        aux=1;
    else
      xk0=xk;
      %pause
      %r=10;
    end
    k=k+1;
    %pause
end
timen=cputime-tt;
lamb
difproj = xk - (P*xk + (b'*xk)*r*b - a - lar*b);
difproj=max(l,min(difproj,u));
erro3=norm(xk-difproj)
%xk
k;
%xk=disp([xk(1); xk(2)])
%disp([erro1 erro2 timen])
%errot = norm(p.*xs - a + lambxs.ineqlin*b + lambxs.upper - lambxs.lower)
%disp([k]);

hold on
%plot(xplot,yplot,'.', 'markersize', 8);
plot(xplot,yplot,'.',xpplot,ypplot,'.');
hold off



%comparando quadprog:
%[xquadprog,fquad,~,out]=fmincon(@(x) x'*(P*x)-a'*x,ones(n,1),[],[],b',[c],l,u);
xquadprog = quadprog(P,-a,[],[],b',c,l,u);


