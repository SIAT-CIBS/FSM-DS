function [ total ] = estimate_accuracy( x_de,x_esde,r,q,epsion)
%x_de and x_esde is d*N matrix where d is dimenson and N is size

[s_1,s_2]=size(x_de);
[s_3,s_4]=size(x_esde);
if s_1~=s_3||s_2~=s_4
    error('input wrong');
end

T=s_4;
total=0;
for i=1:T
    tmp=x_de(:,i)'*x_esde(:,i);
    div=(norm(x_de(:,i))*norm(x_esde(:,i))+epsion);
    tmp=(1-tmp/div)^2;
    total_tmp=r*tmp;
    tmp=(x_de(:,i)-x_esde(:,i))'*(x_de(:,i)-x_esde(:,i));
    tmp=tmp/div;
    total_tmp=total_tmp+q*tmp;
    total=total+sqrt(total_tmp);
end
total=total/T;
end

