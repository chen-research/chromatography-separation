
clear all
nCompIn1=2;
timeN = 800;
num = 1;
%The ith set of random number generated are always the same set
run_times = 1; 

OutputMatrix=zeros(num,800);
%parameterMatrix = rand(num,10);
parameterMatrix = rand(num,6);

%%%%%%%%%%%%%%%%%%%%%%%%
%Adsorption Isotherm Parameters
%parameterMatrix(num,1:8) = [89.509 11.44  89.66  15.213 10.009 47.647 23.39  58.153];
%parameterMatrix(num,1:4) = [89.854 45.256 10.002 53.23];
parameterMatrix(num,1:4) = [89.78   45.06   10.84   27.49 ];
%Injection Profile
%parameterMatrix(num,9:10) = [5.000     5.001]; 
%parameterMatrix(num,5:6) = [5.0     5.899]; 
parameterMatrix(num,5:6) = [5.0  5.06]; 
%%%%%%%%%%%%%%%%%%%%%%%%

i = 1;

parameterIn1= parameterMatrix(i,1:4);
Injection= parameterMatrix(i,5:6);


parameterIn2=[parameterIn1(1:2);parameterIn1(3:4)];
[t1, R_obs1, Tsol] = parameter_to_noisyData(nCompIn1,parameterIn2,Injection);

t_new=linspace(2,Tsol/60,timeN);
R_new = interp1(t1,R_obs1,t_new);

OutputMatrix(i,:)=R_new;

%parameterMatrix
plot(OutputMatrix)
xlabel('time')




