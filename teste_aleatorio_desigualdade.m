clear all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dados do problema  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=1000000;
prec = 1e-4;
%p = ones(n,1);

%% gracao de dados
% p= matiz diadonal problema
p=10 + 15.*rand(n,1); % p=10 + (25-10).*rand(10,1) U[10,25];
%b = p;
b=10 + 15.*rand(n,1); % p=10 + (25-10).*rand(10,1) U[10,25];
%a = (2:n+1)';
a=10 + 15.*rand(n,1); % p=10 + (25-10).*rand(10,1) U[10,25];

B=1 + 14.*rand(2,n);
l=min(B)';
u=max(B)';
btl=l'*b;
btu=u'*b;
c=btl+(btu-btl)*rand(1);

disp('gerou dados')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dados inciais  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xk0=ones(n,1); % Ponto inicial
lamb=10;        % Multiplicador inicial a restricao de igualdade b'*x-c
v=ones(n,1);   % multiplicador associado a restricao l-x
w=v;         % multiplicador associado a restricao x-u
r=1;          % Parametro de penalidade inicial
%r=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop Principal  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tt=cputime;     % Time
for k = 1:600000
  
% %   bt=sqrt(r)*b;
% %   alfa=1/(1+bt'*bt);
% %   IP= eye(n)-(alfa*(bt*bt')); % inversa de P formula de Sherman-Morrisson
% %   xk=IP*(a-((lamb-r*c)*b)); % Solucao de grad(lag)=0
s=max(-lamb/r-b'*xk0+c,0);
if s == 0
    bt=sqrt(r)*b;
    vec=bt./p;
    alfa=1/(1+bt'*vec);
    lar=lamb-r*c;
    xk=(a-lar*b)./p-alfa*((vec'*a)*vec-lar*(vec'*b)*vec);
else
    xk=a./p;
end


xk=max(l,min(xk,u));

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
    rest=b'*xk-c+s;
    
    
    %if abs(rest)>0.1*abs(rest0)
     %   r=2*r
    %end   
        
    %r=n/(4*k+1);
    
    %lamb0=lamb;
    lamb=lamb+r*rest;
    %r=abs(lamb-lamb0)/abs(rest0)
      k=k+1;
      
      %% Criterio de parada
    erro1=norm(xk-xk0);
    erro2=abs(rest)
    if erro1 <= prec && erro2 <= prec
        break
    else
      xk0=xk;  
      %r=10;
    end
    
end
timen=cputime-tt
s
lamb
%xk
k
%xk=disp([xk(1); xk(2)])
%disp([erro1 erro2 timen])
%errot = norm(p.*xs - a + lambxs.ineqlin*b + lambxs.upper - lambxs.lower)
%disp([k]);