
clear all
nCompIn1=2;
timeN = 800;
num = 1;
%The ith set of random number generated are always the same set
run_times = 1; 

OutputMatrix=zeros(num,800);
parameterMatrix = rand(num,10);

%%%%%%%%%%%%%%%%%%%%%%%%
%Adsorption Isotherm Parameters
%parameterMatrix(num,1:8) = [89.509 11.44  89.66  15.213 10.009 47.647 23.39  58.153];
%parameterMatrix(num,1:8) = [12.08945091 41.22400857 10.22242404 59.19095629 81.01037426 46.68415624 85.66629336 44.77811191];
parameterMatrix(num,1:8) = [88.23  36.56  89.32  31.53  15.63  35.68  18.38  34.63 ];
%Injection Profile
%parameterMatrix(num,9:10) = [5.000     5.001]; 
%parameterMatrix(num,9:10) = [5.01434745 18.73013587];
parameterMatrix(num,9:10) = [5.25  5.02];
%%%%%%%%%%%%%%%%%%%%%%%%


i = 1;

parameterIn1= parameterMatrix(i,1:8);
Injection= parameterMatrix(i,9:10);


parameterIn2=[parameterIn1(1:4);parameterIn1(5:8)];
[t1, R_obs1, Tsol] = parameter_to_noisyData(nCompIn1,parameterIn2,Injection);

t_new=linspace(2,Tsol/60,timeN);
R_new = interp1(t1,R_obs1,t_new);

OutputMatrix(i,:)=R_new;

%parameterMatrix
plot(OutputMatrix)
xlabel('time')




