function [B, accuracy]=weakCondition(H,T,W,B0,P0,options)
% This function provide solution to the weak condition to make a stable
% sense to the system.
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
%
%Outputs:
% B: The resulted output matrix
% accuracy: the computed norm error.
%
%
B0=B0';
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
if ~isfield(options,'tol_mat_bias')
    options.tol_mat_bias=10^(-15);
end
if options.display
        str = 'iter';
    else
        str = 'off';
end


[s_1,s_2]=size(W);
[s_3,s_4]=size(B0);

%ensure dimension equals to 2
if s_2~=s_3||s_1~=s_4
    error('input wrong');
end
Dimenson=s_1;

hidden=s_2;

[s_1,s_2]=size(H);
[s_3,s_4]=size(T);
if s_1~=s_3||s_2~=hidden||s_4~=Dimenson
    error('input wrong');
end

%min||HB-O|| s.t (w_i)'*(b_i)<0 for i=1,...,hidden
%prepare initial value
b=B0';
b=reshape(b,1,hidden*Dimenson);
b=b';
obj_handle = @(p) obj(p,H,T);
if isempty(P0)
    ctr_handle = @(p) ctr_eigenvalue(W,p,options.tol_mat_bias);
else
    %ctr_handle = @(p) ctr_eigenvalue_P(W,p,P0,options.tol_mat_bias);
    ctr_handle = @(p) ctr_eigenvalue(W,p,options.tol_mat_bias);
end
 
 % Options for NLP Solvers
    optNLP = optimset( 'Algorithm', 'interior-point', 'LargeScale', 'off',...
        'GradObj', 'on','GradConstr', 'on',  'DerivativeCheck', 'off', ...
        'Display', 'iter', 'TolX', 1e-12, 'TolFun', 1e-12, 'TolCon', 1e-12, ...
        'MaxFunEval', 200000, 'MaxIter', options.max_iter, 'DiffMinChange', ...
         options.tol_stopping, 'Hessian','off','display',str);
     B = fmincon(obj_handle, b,[],[],[],[],[],[],ctr_handle,optNLP);
     accuracy=obj(B,H,T);
     B=B';
     B=reshape(B,Dimenson,hidden);
     
     
     function [val d_beta]=obj(beta,H,O)
         Dimension=size(O,2);
         hidden=size(H,2);
         beta=beta';
         beta=reshape(beta,Dimension,hidden);
         beta=beta';
         val=norm(H*beta-O);
         dJ=H'*(H*beta-O);
         [c_1,c_2]=size(dJ);
         dJ=dJ';
         d_beta=reshape(dJ,1,c_1*c_2);
         d_beta=d_beta';

    function [c ceq dc dceq]=ctr_eigenvalue(W,beta,tol_mat_bias)
        ceq=[];
        dceq=[];
        
        [Dimension,hidden]=size(W);
        c=zeros(Dimension,1);
        dc=zeros(Dimension*hidden,Dimension);
        beta=reshape(beta,Dimension,hidden);
        beta=beta';
        beta_copy=zeros(hidden,Dimension);
        B=W*beta+(W*beta)';
        rArs=zeros(Dimension,Dimension);
        for j=1:Dimension
             dJ=zeros(j,j);
            if (-1)^(j+1)*det(B(1:j,1:j))+(tol_mat_bias)^(j/Dimension) > 0
                c(j)=(-1)^(j+1)*det(B(1:j,1:j))+(tol_mat_bias)^(j/Dimension);
            end
            
            for i1=1:j
                for i2=1:j
                    rArs = rArs * 0;
                    rArs (i2,i1) = 1;
                    rBrs = rArs + rArs';
                    if j==1
                        dJ(1,1) = rBrs(1,1);
                    else
                        tmp = det(B(1:j,1:j));
                        if abs(tmp) > 1e-10
                            term = trace(B(1:j,1:j)\rBrs(1:j,1:j))*tmp;
                        else
                            term = trace(adjugate(B(1:j,1:j))*rBrs(1:j,1:j));
                        end
                        dJ(i1,i2)=(-1)^(j+1)*term;
                    end
                end
            end
            beta_copy=beta_copy*0;
            beta_copy(:,1:j)=W(1:j,:)'*dJ;
            dc(:,j)=reshape(beta_copy',Dimension*hidden,1);
        end
        

    function [c ceq dc dceq]=ctr_eigenvalue_P(W,beta,P,tol_mat_bias)
        ceq=[];
        dceq=[];
        dc=[];
        
        [Dimension,hidden]=size(W);
        %c=zeros(Dimension,1);
        beta=reshape(beta,Dimension,hidden);
        beta=beta';
        B=P*W*beta+(P*W*beta)';
        c=eig(B)+(tol_mat_bias);
        
  
        
