close all
clear;clc;
load('Cshape','demos');
% Pre-processing
dt = 0.1; %The time step of the demonstrations
tol_cutting = 1; % A threshold on velocity that will be used for trimming demos
options.energy='on';
options.suf='off';
[x0 , xT, Data, index] = preprocess_demos(demos,dt,tol_cutting); %preprocessing datas

%% compute the energy matrix
options.energy='on';
if ~isfield(options,'energy') % recompute energy function
    options.energy='off';
end
if ~isfield(options,'suf') % impose sufficient conditions
    options.suf='on';
end

if strcmp(options.energy,'on')
    P0=energy_hu(Data);
else
    P0=[];
end
%% set parameters for optimization
options.max_iter=1000;
options.alpha=0;
options.tol_stopping=10^-10;  % A small positive scalar defining the stoppping
                              % tolerance for the optimization solver [default: 10^-10]
options.display = 1;          % An option to control whether the algorithm
                              % displays the output of each iterations [default: true]
options.alpha=0.05;           % tolerance

start_time_train=cputime;

nhidden=150;
[net,accuracy,H]=mlp_hu_free(Data,nhidden,'sin');
inputweight=net.w1;
InputWeight=inputweight;
outputweight=net.w2;
OutputWeight=outputweight;
TargetData=Data;
lamda=OutputWeight*InputWeight'+InputWeight*OutputWeight';
if lamda(1,1)<0,det(lamda)>0
   InputWeight=InputWeight; 
else
    [InputWeight,pre_accuracy]= LocalCondition(H,TargetData,InputWeight,OutputWeight,P0,options);
end
net.w1=InputWeight;
end_time_train=cputime;
end_time_train=end_time_train-start_time_train;  
d = size(Data,1)/2; %dimension of data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%simulation
%% Simulation
opt_sim.dt = 0.1;
opt_sim.i_max = 3000;
opt_sim.tol = 0.1;
x0_all = Data(1:d,index(1:end-1)); %finding initial points of all demonstrations
fn_handle = @(x) mlpfwd(net,x);
[x xd]=Simulation(x0_all,[],fn_handle,opt_sim); %running the simulator
x_all=[];
for i=1:size(x,3)
    x_all=[x_all x(:,:,i)];
end

%% Interpolation of demostrations
N_rep=size(x,2);
N_dem=size(demos{1,1},2);
x_dem=1:1:N_dem;
xq=N_dem/N_rep:N_dem/N_rep:N_dem;
Data_vq_x=[];Data_vq_y=[];
for i=1:size(x,3)

    traj_x(i,:)=demos{1,i}(1,:);
    traj_y(i,:)=demos{1,i}(2,:);
    vq_x(i,:)=interp1(x_dem,traj_x(i,:),xq,'spline');
    vq_y(i,:)=interp1(x_dem,traj_y(i,:),xq,'spline');
    demos_vq{i}=[vq_x(i,:),vq_y(i,:)];
Data_vq_x=[Data_vq_x,vq_x(i,:)];
Data_vq_y=[Data_vq_y,vq_y(i,:)];
end
Data_vq=[Data_vq_x;Data_vq_y];
n=size(x_all,2)/size(x,3);
for i=1:size(x,3)
    x_all_SEA(:,:,i)=x_all(:,n*(i-1)+1:n*i) 
end
x=x_all_SEA;
%%
% plotting the result
figure('name','Results from Simulation','position',[265   200   520   720])
sp(1)=subplot(3,1,1);
hold on; box on
plot(Data(1,:),Data(2,:),'r.')
xlabel('$\xi_1 (mm)$','interpreter','latex','fontsize',25);
ylabel('$\xi_2 (mm)$','interpreter','latex','fontsize',25);
title('Simulation Results')

sp(2)=subplot(3,1,2);
hold on; box on
plot(Data(1,:),Data(3,:),'r.')
xlabel('$\xi_1 (mm)$','interpreter','latex','fontsize',25);
ylabel('$\dot{\xi}_1 (mm/s)$','interpreter','latex','fontsize',25);

sp(3)=subplot(3,1,3);
hold on; box on
plot(Data(2,:),Data(4,:),'r.')
xlabel('$\xi_2 (mm)$','interpreter','latex','fontsize',25);
ylabel('$\dot{\xi}_2 (mm/s)$','interpreter','latex','fontsize',15);


for i=1:size(x,3)
    plot(sp(1),x(1,:,i),x(2,:,i),'linewidth',2)
    plot(sp(2),x(1,:,i),xd(1,:,i),'linewidth',2)
    plot(sp(3),x(2,:,i),xd(2,:,i),'linewidth',2)
    plot(sp(1),x(1,1,i),x(2,1,i),'ok','markersize',5,'linewidth',5)
    plot(sp(2),x(1,1,i),xd(1,1,i),'ok','markersize',5,'linewidth',5)
    plot(sp(3),x(2,1,i),xd(2,1,i),'ok','markersize',5,'linewidth',5)
end
for i=1:3
    axis(sp(i),'tight')
    ax=get(sp(i));
    axis(sp(i),...
        [ax.XLim(1)-(ax.XLim(2)-ax.XLim(1))/10 ax.XLim(2)+(ax.XLim(2)-ax.XLim(1))/10 ...
        ax.YLim(1)-(ax.YLim(2)-ax.YLim(1))/10 ax.YLim(2)+(ax.YLim(2)-ax.YLim(1))/10]);
    plot(sp(i),0,0,'k*','markersize',15,'linewidth',3)
    if i==1
        D = axis(sp(i));
    end
end

%% plotting streamlines
figure('name','Streamlines','position',[800   90   560   320])
% plotStreamLines(net,D);
hold on
plot(Data_vq(1,:),Data_vq(2,:),'r.','linewidth',1.5);
legend('Training Date')
hold on
plot(0,0,'k*','markersize',15,'linewidth',3);
n=size(x_all,2)/size(x,3);
for i=1:size(x,3)
    plot(x_all(1,n*(i-1)+1:n*i),x_all(2,n*(i-1)+1:n*i),'k','linewidth',2)
    plot(x(1,1,i),x(2,1,i),'ok','markersize',5,'linewidth',5);
end
set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[])

