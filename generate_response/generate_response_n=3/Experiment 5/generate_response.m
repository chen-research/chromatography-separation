
clear all
nCompIn1=3;
timeN = 800;
num = 1;
%The ith set of random number generated are always the same set
run_times = 1; 

OutputMatrix=zeros(num,800);
%parameterMatrix = rand(num,10);
%parameterMatrix = rand(num,6);
parameterMatrix = rand(num,9);

%%%%%%%%%%%%%%%%%%%%%%%%
%Adsorption Isotherm Parameters
%parameterMatrix(num,1:8) = [89.509 11.44  89.66  15.213 10.009 47.647 23.39  58.153];
%parameterMatrix(num,1:4) = [10.91270638 17.30455837 89.74561594 14.81492469];
%parameterMatrix(num,1:6) = [89.98, 53.30, 10.0, 54.27, 49.87, 53.96];
parameterMatrix(num,1:6) = [89.95  11.72  52.60  89.91  10.79   29.87];

%Injection Profile
%parameterMatrix(num,9:10) = [5.000     5.001]; 
%parameterMatrix(num,5:6) = [5.00158084  5.00244998]; 
%parameterMatrix(num,7:9) = [13.72, 15.83, 5.44];
parameterMatrix(num,7:9) = [5.51  5.00  9.12];
%%%%%%%%%%%%%%%%%%%%%%%%


i = 1;

% parameterIn1= parameterMatrix(i,1:4);
% Injection= parameterMatrix(i,5:6);
parameterIn1= parameterMatrix(i,1:6);
Injection= parameterMatrix(i,7:9);


%parameterIn2=[parameterIn1(1:2);parameterIn1(3:4)];
parameterIn2=[parameterIn1(1:2);parameterIn1(3:4);parameterIn1(5:6)];

[t1, R_obs1, Tsol] = parameter_to_noisyData(nCompIn1,parameterIn2,Injection);

t_new=linspace(2,Tsol/60,timeN);
R_new = interp1(t1,R_obs1,t_new);

OutputMatrix(i,:)=R_new;

%parameterMatrix
plot(OutputMatrix)
xlabel('time')




