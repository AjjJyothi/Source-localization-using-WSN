function F = mlobjfun(x,locSen,Noise_Var,d)
N=length(d);
sig=Noise_Var*(((1/2)*ones(N))+((1/2)*eye(N)));
k=inv(sig);
mod=sqrt(det(sig));
 for i=1:N
     h(i)=(norm(x-locSen(i+1,:)))-(norm(x-locSen(1,:)))-d(i);
 end
x1=h*k*h';
%F=(exp((-1/2)*x1)/(2*(pi^(N/2))*mod));
F=x1;
end