function [X,fval,lambda,iter] = wright(P,a,l,u,X0,b,c) 
%   tentando wright p resolver
%   problem:
%
%            min x'*H*x + f'*x   subject to:  b'x = c 
%             l<=x<=u    
%
prec = 1e-4; n=length(X0); 
x=X0; s=u-X0; lambda=1; mul=ones(n,1); muu=mul;
iter=0;
gap = (x'*mul-s'*muu)/(2*n+iter); tau=0.05*gap;
if nargin>5
    r1 = 2*(P*x)-a + lambda*b - mul + muu;
else
    r1 = 2*(P*x)-a - mul + muu;
end
r2 = diag(x-l)*mul;
r3 = diag(s)*muu;
if nargin>5
 r4 = b'*x-c;
else
    %dai o r4 eh so p passar no teste
    r4=0;
end

while ( (norm(r1)>prec)||(norm(r2)>prec)||(norm(r3)>prec)||(norm(r4)>prec)  )&&(iter<6000)
    

    % computa vetor
    if nargin>5
        rf = zeros(n+1,1);
    else
        rf = zeros(n,1);
    end
    rf(1:n)=r1 - r2./(x-l) + r3./s;
    
    %computa matriz (so do x e lambda)
    if nargin>5
        matrizk = zeros(n+1,n+1);
        matrizk(1:n,1:n) = 2*P + diag(mul./(x-l) - muu./s);
        matrizk(1:n,n+1) = b;
        matrizk(n+1,1:n) = b';
    else
        matrizk(1:n,1:n) = 2*P + diag(mul./(x-l) - muu./s);
    end
    
    %matrizk(n+1:2*n,1:n) = -eye(n,n);
    %matrizk(n+1:2*n,n+2:2*n+1) = eye(n,n);
    %matrizk(n+1,1:n) = b';
    %matrizk(n+2:2*n+1,1:n) = -diag(muu);
    %matrizk(n+2:2*n+1,1:n) = -diag(muu);
    %matrizk(2*n+2:3*n+1,n+2:2*n+1) = diag(s);
    %matrizk(2*n+2:3*n+1,n+2:2*n+1) = diag(x-l);
    %pause
    
    vetordeltas = (matrizk)\(rf);
    
    deltax = vetordeltas(1:n);
    if nargin>5
        deltarho = vetordeltas(n+1);
    end
    
    deltas = -deltax;
    
    deltamul = (-r3-muu.*deltax)./s;
    
    deltamuu = (-r2-mul.*deltax)./(x-l);

    
    alfa=1; aux=0;
    while (aux==0)&&(alfa>1e-300)
        
        xplus = x + alfa*deltax;
        splus = s + alfa*deltas;
        mulplus = mul + alfa*deltamul;
        muuplus = muu + alfa*deltamuu;
        if nargin>5
            lambdaplus = lambda + alfa*deltarho;
        end
        aux=1;
        for j=1:n
         if (aux==1)&&( (xplus(j)<l(j))||(splus(j)<0)||(mulplus(j)<0)||(muuplus(j)<0)   )
            aux=0; alfa=alfa/2;
         end        
        end
    end %while
    x=xplus; 
    s=splus; 
    muu=muuplus; 
    mul=mulplus; 
    if nargin>5
        lambda=lambdaplus;
    end
    
    gap = (x'*mul-s'*muu)/(2*n+iter); tau=0.05*gap;
    if nargin>5
        r1 = 2*(P*x)-a + lambda*b - mul + muu;
    else
        r1 = 2*(P*x)-a - mul + muu;
    end
    r2 = diag(x-l)*mul-tau*ones(n,1);
    r3 = diag(s)*muu-tau*ones(n,1);
    if nargin>5
        r4 = b'*x-c;
    else
        r4 = 0;
    end
    iter=iter+1;
end %end while

X = x;
fval = x'*((P*x)-a);

end %endfunction