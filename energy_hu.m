function P0=energy_hu(Data)
% Entrance function to learn the P0 which allows the energy function to be
% x^T*P0*x. We use the CLF-DM codes directly.

%% Setting parameters of the Control Lyapunov Function
clear Vxf0
Vxf0.L = 0; %the number of asymmetric quadratic components L=0.
Vxf0.d = size(Data,1)/2;
Vxf0.w = 1e-4; %A positive scalar weight regulating the priority between the 
               %two objectives of the opitmization. Please refer to the
               %page 7 of the paper for further information.

% A set of options that will be passed to the solver. Please type 
% 'doc preprocess_demos' in the MATLAB command window to get detailed
% information about other possible options.
options.tol_mat_bias = 10^-1; % a very small positive scalar to avoid
                              % having a zero eigen value in matrices P^l [default: 10^-15]
                              
options.display = 1;          % An option to control whether the algorithm
                              % displays the output of each iterations [default: true]
                              
options.tol_stopping=10^-10;  % A small positive scalar defining the stoppping
                              % tolerance for the optimization solver [default: 10^-10]

options.max_iter = 500;       % Maximum number of iteration for the solver [default: i_max=1000]

options.optimizePriors = true;% This is an added feature that is not reported in the paper. In fact
                              % the new CLFDM model now allows to add a prior weight to each quadratic
                              % energy term. IF optimizePriors sets to false, unifrom weight is considered;
                              % otherwise, it will be optimized by the sovler.
                              
options.upperBoundEigenValue = true; %This is also another added feature that is impelemnted recently.
                                     %When set to true, it forces the sum of eigenvalues of each P^l 
                                     %matrix to be equal one. 



%% Estimating an initial guess for the Lyapunov function
b_initRandom = false;
if b_initRandom
    lengthScale = sqrt(var(Data(1:Vxf0.d,:)'));
    lengthScaleMatrix = sqrtm(cov(Data(1:Vxf0.d,:)'));
    lengthScale = lengthScale(:);
    Vxf0.Priors = rand(Vxf0.L+1,1);
    for l=1:Vxf0.L+1
        tmpMat = randn(Vxf0.d,Vxf0.d);
        Vxf0.Mu(:,l) = randn(Vxf0.d,1).*lengthScale;
        Vxf0.P(:,:,l) = lengthScaleMatrix*(tmpMat*tmpMat')*lengthScaleMatrix;
    end
else
    Vxf0.Priors = ones(Vxf0.L+1,1);
    Vxf0.Priors = Vxf0.Priors/sum(Vxf0.Priors);
    Vxf0.Mu = zeros(Vxf0.d,Vxf0.L+1);
    for l=1:Vxf0.L+1
        Vxf0.P(:,:,l) = eye(Vxf0.d);
    end
end

% Solving the optimization
Vxf = learnEnergy(Vxf0,Data,options);
P0=Vxf.P(:,:,1);
if Vxf.d~=2
    return
end
%Vxf.P(:,:,1)=eye(2); % comparison with quadratic energy function.
%% Plotting the result
fig = figure;
sp = gca;
hold on
h(1) = plot(sp,Data(1,:),Data(2,:),'r.');
axis tight
ax=get(gca);
axis([ax.XLim(1)-(ax.XLim(2)-ax.XLim(1))/10 ax.XLim(2)+(ax.XLim(2)-ax.XLim(1))/10 ...
      ax.YLim(1)-(ax.YLim(2)-ax.YLim(1))/10 ax.YLim(2)+(ax.YLim(2)-ax.YLim(1))/10]);

h(3) = EnergyContour(Vxf,axis,[],[],sp, [], false);
h(2) = plot(0,0,'g*','markersize',15,'linewidth',3);

%xlabel('x (mm)','fontsize',15);
%ylabel('y (mm)','fontsize',15);
%title('Energy Levels of the learned Lyapunov Functions')
legend(h,'demonstrations','target','energy levels','location','southwest')