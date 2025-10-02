
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
%parameterMatrix(num,1:8) = [89.87686235 50.77844858 88.37793873 43.65747114 13.60154007 61.64685252 16.36673544 43.99300176];
parameterMatrix(num,1:8) = [88.35 33.03 89.21 53.24  10.69  36.73 15.23 65.12  ];

%Injection Profile
%parameterMatrix(num,9:10) = [5.000     5.001];  
%parameterMatrix(num,9:10) = [5.70819304  5.00771512]; 
parameterMatrix(num,9:10) = [6.06   5.01]; 
%%%%%%%%%%%%%%%%%%%%%%%%
%[89.87686235 50.77844858 88.37793873 43.65747114 13.60154007 61.64685252 16.36673544 43.99300176]  5.70819304  5.00771512

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




