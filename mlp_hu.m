function [net accuracy] = mlp_hu(TrainingData, nhidden,ActivationFunction,options)
% Entrance function for sufficient or weak conditions to ensure the
% globally asymptotical stability of the dynamical system. 
% Inputs:
%   TraniningData:  A 2d x N_Total matrix containing all demonstration data points.
%   nhidden:    The number of the hidden neural layer
%   ActivationFunction: The active function in the neuron
%   options: A structure specify whther quadratic energy functions are used
%   and sufficient or weak conditions are applied.
%       - energy: on - using energy function x^T*P*x
%                 off - using 0.5x^Tx
%       - suf: on - using sufficient condtions


%% compute the energy matrix
if ~isfield(options,'energy') % recompute energy function
    options.energy='off';
end
if ~isfield(options,'suf') % impose sufficient conditions
    options.suf='on';
end

%%
N=size(TrainingData,2);
nin=size(TrainingData,1)/2;

Inputs=TrainingData(1:nin,:)';
TargetData=TrainingData(nin+1:end,:)';
if strcmp(options.energy,'on')
    P0=energy_hu(TrainingData);
else
    P0=[];
end

%set pre-data by randomly
InputWeight=rand(nin,nhidden)-1;
adjust=max(abs(Inputs),[],1);
InputWeight=InputWeight./repmat(adjust',1,nhidden);
bias=0;% 
tempH=Inputs*InputWeight+bias;

net.type='elm';
net.nin=nin;
net.nout=nin;
net.w1=InputWeight;
net.b1=bias;
net.outfn='linear';

%%%%%%%%%%% Calculate hidden neuron output matrix H
switch lower(ActivationFunction)
    case {'sig','sigmoid'}
        %%%%%%%% Sigmoid 
        H = 2 ./ (1 + exp(-tempH))-1;
        net.activefn='sig';
    case {'sin','sine'}
        %%%%%%%% Sine
        H = sin(tempH);  
        net.activefn='sin';
    case {'tanh'}
        %%%%%%%% tanh
        H=tanh(tempH);
        net.activefn='tanh';
    case {'sinh'}
        H=sinh(tempH);
        net.activefn='sinh';
    case {'hardlim'}
        %%%%%%%% Hard Limit
        H = hardlim(tempH); 
        net.activefn='hardlim';
        %%%%%%%% More activation functions can be added here 
    otherwise
        error('wrong active function!');
end
clear tempH;                                        %   Release the temparary array for calculation of hidden neuron output matrix H


%% set parameters for optimization
options.max_iter=1000;
options.alpha=0;
options.tol_stopping=10^-10;  % A small positive scalar defining the stoppping
                              % tolerance for the optimization solver [default: 10^-10]
options.display = 1;          % An option to control whether the algorithm
                              % displays the output of each iterations [default: true]
options.alpha=0.05;           % tolerance
                              
%%compute the initial output weights
OutputWeight=pinv(H) * TargetData;
OutputWeight=OutputWeight';
lamda=-sum(sum(InputWeight.*OutputWeight))/sum(sum(InputWeight.*InputWeight));
if lamda>0
    OutputWeight=-lamda*InputWeight;
else
    OutputWeight=zeros(size(InputWeight));
end
%OutputWeight=-InputWeight;

[pre_B,pre_accuracy]= sufficientCondition(H,TargetData,InputWeight,OutputWeight,P0,options);
if ~strcmp(options.suf,'on')
    [pre_B,pre_accuracy]= weakCondition(H,TargetData,InputWeight,pre_B,P0,options);
end

net.w2=pre_B;
net.b2=0;
accuracy=pre_accuracy;
