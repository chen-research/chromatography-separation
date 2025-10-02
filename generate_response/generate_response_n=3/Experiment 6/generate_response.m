
clear all
nCompIn1=3;
timeN = 800;
num = 1;
%The ith set of random number generated are always the same set
run_times = 1; 

OutputMatrix=zeros(num,800);
% parameterMatrix = rand(num,10);
parameterMatrix = rand(num,15);

%%%%%%%%%%%%%%%%%%%%%%%%
%Adsorption Isotherm Parameters
%parameterMatrix(num,1:8) = [89.509 11.44  89.66  15.213 10.009 47.647 23.39  58.153];
%parameterMatrix(num,1:8) = [12.08945091 41.22400857 10.22242404 59.19095629 81.01037426 46.68415624 85.66629336 44.77811191];
%parameterMatrix(num,1:12) = [48.26, 48.02, 60.39, 45.03, 89.97, 53.40, 89.95, 41.06, 16.11,  34.46, 13.24, 51.48];
parameterMatrix(num,1:12) = [18.42, 41.42, 14.24,  51.91, 83.34, 49.62, 84.26,  62.25, 45.78, 42.86, 47.72, 52.67];
%Injection Profile
%parameterMatrix(num,9:10) = [5.000     5.001]; 
%parameterMatrix(num,9:10) = [5.01434745 18.73013587];
%parameterMatrix(num,13:15) = [5.0, 16.94, 5.0];
parameterMatrix(num,13:15) = [7.07, 16.37,  5.93];
%%%%%%%%%%%%%%%%%%%%%%%%


i = 1;

% parameterIn1= parameterMatrix(i,1:8);
% Injection= parameterMatrix(i,9:10);
parameterIn1= parameterMatrix(i,1:12);
Injection= parameterMatrix(i,13:15);


%parameterIn2=[parameterIn1(1:4);parameterIn1(5:8)];
parameterIn2=[parameterIn1(1:4);parameterIn1(5:8);parameterIn1(9:12)];
[t1, R_obs1, Tsol] = parameter_to_noisyData(nCompIn1,parameterIn2,Injection);

t_new=linspace(2,Tsol/60,timeN);
R_new = interp1(t1,R_obs1,t_new);

OutputMatrix(i,:)=R_new;

%parameterMatrix
plot(OutputMatrix)
xlabel('time')




