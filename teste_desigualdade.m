clear all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dados do problema  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=100;
prec = 1e-5;
prec2= 1e-4;
p = ones(n,1);
b = p;
a = (2:n+1)';
l = zeros(n,1);
u = 10*p;
c = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dados inciais  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xk0=ones(n,1); % Ponto inicial
lamb=-100;        % Multiplicador inicial a restricao de igualdade b'*x-c
v=ones(n,1);   % multiplicador associado a restricao l-x
w=v;         % multiplicador associado a restricao x-u
%r=1;          % Parametro de penalidade inicial
r=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop Principal  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tt=cputime;     % Time
for k = 1:6000000 
% %   bt=sqrt(r)*b;
% %   alfa=1/(1+bt'*bt);
% %   IP= eye(n)-(alfa*(bt*bt')); % inversa de P formula de Sherman-Morrisson
% %   xk=IP*(a-((lamb-r*c)*b)); % Solucao de grad(lag)=0
s=max(-lamb/r-b'*xk0+c,0);
if s == 0
  bt=sqrt(r)*b;
  alfa=1/(1+bt'*bt);
  lar=lamb-r*c;
  xk=a-lar*b-alfa*((bt'*a)*bt-lar*(bt'*b)*bt);
else
  xk=a;
end

 %% Proje????o da solucao na caixa
 xk=max(l,min(xk,u));
    %for i=1:n
     %   if xk(i) < l(i) 
      %      xk(i)=l(i);   % Componente na solu????o fora da caixa do lado esquerdo
       % elseif xk(i) > u(i) 
        %    xk(i)=u(i);   % Componente na solu????o fora da caixa do lado direito
        %end
    %end
    
  % Atualizamos os multiplicadores das restri????es da caixa mul, muu
  %gradlagx2=xk-a+(lamb+r*(b'*xk-c))*b; %avaliamos o gradiente da lagrangeana em mul
   
   %for i=1:n
    %    if xk(i) <= l(i) 
     %     v(i)=xk(i)-a(i)+(lamb+r*(b'*xk-c))*b(i); % multiplicador v inviavel
      %  else 
       %     v(i)=0; % N??o precissa associado a restricao x-u
        %end
   % end
   
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
 
    rest=n/10*(b'*xk-c);
    lamb=lamb+r*rest;
      k=k+1;
      
      %% Criterio de parada
    erro1=norm(xk-xk0);
    erro2=abs(rest);
    if erro1 <= prec && erro2 <= prec2
        break
    else
      xk0=xk;  
     % r=min(r*2,100);
    end
    
end
%xk
timen=cputime-tt
k
%xk=disp([xk(1); xk(2)])
%disp([erro1 erro2 timen])
%errot = norm(p.*xs - a + lambxs.ineqlin*b + lambxs.upper - lambxs.lower)
%disp([k]);