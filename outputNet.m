function []=outputNet(net)
% Syntax:
%       outputNet(net)
% Output the net structure into a file.

fid=fopen('netData.txt','wt');
%%%%%%%%%%% Calculate hidden neuron output matrix H
switch lower(net.activefn)
    case {'sig','sigmoid'}
        %%%%%%%% Sigmoid 
        fprintf(fid,'%d\n',0);
    case {'sin','sine'}
        %%%%%%%% Sine
        fprintf(fid,'%d\n',1);
    case {'tanh'}
        %%%%%%%% tanh
        fprintf(fid,'%d\n',2);
    case {'sinh'}
        fprintf(fid,'%d\n',3);
    case {'hardlim'}
        %%%%%%%% Hard Limit
        fprintf(fid,'%d\n',4);
    otherwise
        error('wrong active function!');
end
[m,n]=size(net.w1);
fprintf(fid,'%d %d %d\n',m,n,m);
for i=1:m
    for j=1:n
        fprintf(fid,'%f ',net.w1(i,j));
    end
    fprintf(fid,'\n');
end

for i=1:m
    for j=1:n
        fprintf(fid,'%f ',net.w2(i,j));
    end
    fprintf(fid,'\n');
end

fprintf(fid,'%f',net.b1);
fclose(fid);