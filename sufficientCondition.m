function [B, accuracy]=sufficientCondition(H,T,W,B0,P0,options)
% This function provide a solution to the suffcient constraints to ensure 
% the gloablly asymptotical stability of the dynamical system.
%Inputs:
% options: A structure to set the optional parameters of the solver.
%              The following parameters can be set in the options:
%
%       - .tol_stopping:     A small positive scalar defining the stoppping
%                            tolerance for the optimization solver [default: 10^-10]
%
%       - .max_iter:            maximum number of iteration for the solver [default: i_max=1000]
%       - .alpha:               Sylverster's criterion [default: alpha=0.05]
%
% H: The temporary hidden matrix with dimension N*nhidden
% T: Targets from demonstrations with dimension N*d
% W: Input matrix with dimension d*nhidden
% B0: Initial output matrix with dimension d*nhidden
% P0: The adjust matrix which changes the energy function into f=x^T*P0*x
%
%Outputs:
% B: The resulted output matrix
% accuracy: the computed norm error.
%



if ~isfield(options,'max_iter') % maximum number of iterations
    options.max_iter=1000;
end
if ~isfield(options,'alpha') % tolerance for constraint
    options.alpha=0.05;
end
if ~isfield(options,'tol') % tolerance for objective value
    options.tol=0.05;
end
if ~isfield(options,'dt') % step value
    options.dt=0.05;
end
if options.display
        str = 'iter';
    else
        str = 'off';
end


[Dimension,hidden]=size(W);
[s_3,s_4]=size(B0);

%ensure dimension equals to 2
if hidden~=s_4||Dimension~=s_3
    error('input wrong in sufficientCondition');
end
[s_1,s_2]=size(H);
[s_3,s_4]=size(T);
if s_1~=s_3||s_2~=hidden||s_4~=Dimension
    error('input wrong in sufficientCondition');
end
if isempty(P0)
    P0=eye(Dimension);
end


%min||HB-O|| s.t (b_i)*(w_i)'+(w_i)*(b_i)' negative semi-definite for i=1,...,hidden

%prepare initial value
b=reshape(B0,hidden*Dimension,1);

optNLP = optimset( 'Algorithm', 'interior-point', 'LargeScale', 'off',...
        'GradObj', 'on',  'DerivativeCheck', 'off', ...
        'Display', 'iter', 'TolX', 1e-12, 'TolFun', 1e-12, 'TolCon', 1e-12, ...
        'MaxFunEval', 200000, 'MaxIter', options.max_iter, 'DiffMinChange', ...
         options.tol_stopping, 'Hessian','off','display',str);
     function_handle= @(x) obj(x,H,T);
     ctr_handle = @(x) ctr_eigenvalue(x,W,P0,Dimension,hidden,options);
     B = fmincon(function_handle, b,[],[],[],[],[],[],ctr_handle,optNLP);
     
     optNLP = optimset( 'Algorithm', 'sqp', 'LargeScale', 'off',...
        'GradObj', 'on',  'DerivativeCheck', 'off', ...
        'Display', 'iter', 'TolX', 1e-12, 'TolFun', 1e-12, 'TolCon', 1e-12, ...
        'MaxFunEval', 200000, 'MaxIter', options.max_iter, 'DiffMinChange', ...
         options.tol_stopping, 'Hessian','off','display',str);
     B = fmincon(function_handle, B,[],[],[],[],[],[],ctr_handle,optNLP);
     
     accuracy=obj(B,H,T);%B: The resulted output matrix %H: The temporary hidden matrix with dimension N*nhidden % T: Targets from demonstrations with dimension N*d
     B=reshape(B,Dimension,hidden);
     
     
    function [c ceq dc dceq]=ctr_eigenvalue(B,W,P0,Dimenson,hidden,options)
        b=reshape(B,Dimenson,hidden);
        W=100*W';
        c=zeros(Dimenson*hidden,1);
        for i=1:hidden
            lambda=eig(P0*b(:,i)*W(i,:)+(P0*b(:,i)*W(i,:))')/2;
            c((i-1)*Dimenson+1:i*Dimenson) = lambda -options.alpha;
        end
        ceq=[];
        dc=[];
        dceq=[];
     
     function [val d_beta]=obj(beta,H,O)
         Dimension=size(O,2);
         hidden=size(H,2);
         beta=reshape(beta,Dimension,hidden);
         beta=beta';
         val=norm(H*beta-O);
         dJ=H'*(H*beta-O);
         dJ=dJ';
         d_beta=reshape(dJ,Dimension*hidden,1);