clear all
nCompIn1=3;
timeN = 800;

%%%%%%%%%%%%%%%%%%%%%
num = 50000;    %number of samples to be generated
OutputMatrix=zeros(num,800);
rng("shuffle") %Set the seed for random number generator %rng(3) %Set the seed for random number generator
%parameterMatrix = rand(num,10);  %generate a matrix of random numbers in [0,1]
parameterMatrix = rand(num,15);  %generate a matrix of random numbers in [0,1]

% parameterMatrix(1:num,1:8) = parameterMatrix(1:num,1:8)*100;
% parameterMatrix(1:num,9:10) = parameterMatrix(1:num,9:10)*30;
parameterMatrix(1:num,1:12) = parameterMatrix(1:num,1:12)*100;
parameterMatrix(1:num,13:15) = parameterMatrix(1:num,13:15)*30;
%%%%%%%%%%%%%%%%%%%%%

core_num = 8; %number of cpu cores to be used for parallel computing
parfor (i=1:num, core_num)  
% parameterIn1= parameterMatrix(i,1:8);
% Injection= parameterMatrix(i,9:10)+0.0001;
% parameterIn2=[parameterIn1(1:4);parameterIn1(5:8)];
parameterIn1= parameterMatrix(i,1:12);
Injection= parameterMatrix(i,13:15)+0.0001;
parameterIn2=[parameterIn1(1:4);parameterIn1(5:8);parameterIn1(9:12)];

[t1, R_obs1, Tsol] = parameter_to_noisyData(nCompIn1,parameterIn2,Injection);
t_new=linspace(2,Tsol/60,timeN);
R_new = interp1(t1,R_obs1,t_new);
OutputMatrix(i,:)=R_new;   %density curve of the ith sample 
end

xlswrite('N3-param15.xlsx',parameterMatrix)
xlswrite('N3-Out15.xlsx',OutputMatrix)





