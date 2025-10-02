clear all
nCompIn1=3;
timeN = 800;

%% 8 parameters [aI_i, bI_i, aII_i, bII_i],i=2 have been reduced to 4 [aI_i, bI_i],i=2. The changes made include ...
%  ...dataGenerateInjection_NonRepeat line 17,21,22,34,35,36


%%%%%%%%%%%%%%%%%%%%%
num = 50000;    %number of samples to be generated
OutputMatrix=zeros(num,800);
rng("shuffle") %Set the seed for random number generator %rng(3) %Set the seed for random number generator

%%%%%%%%%%%%%%%%%%%%%
% down
%parameterMatrix = rand(num,10);  %generate a matrix of random numbers in [0,1]
%parameterMatrix = rand(num,6);  %generate a matrix of random numbers in [0,1]
parameterMatrix = rand(num,9);  %generate a matrix of random numbers in [0,1]

%parameterMatrix(1:num,1:8) = parameterMatrix(1:num,1:8)*100;
%parameterMatrix(1:num,9:10) = parameterMatrix(1:num,9:10)*30;
% parameterMatrix(1:num,1:4) = parameterMatrix(1:num,1:4)*100;
% parameterMatrix(1:num,5:6) = parameterMatrix(1:num,5:6)*30;
parameterMatrix(1:num,1:6) = parameterMatrix(1:num,1:6)*100;
parameterMatrix(1:num,7:9) = parameterMatrix(1:num,7:9)*30;
% up
%%%%%%%%%%%%%%%%%%%%%

core_num = 8; %number of cpu cores to be used for parallel computing


    parfor (i=1:num, core_num)

        %%%%%%%%%%%%%%%%%%%%%
        %  down
        %parameterIn1= parameterMatrix(i,1:8);
        %Injection= parameterMatrix(i,9:10)+0.0001;
        %parameterIn2=[parameterIn1(1:4);parameterIn1(5:8)];
        % parameterIn1= parameterMatrix(i,1:4);
        % Injection= parameterMatrix(i,5:6)+0.0001;
        % parameterIn2=[parameterIn1(1:2);parameterIn1(3:4)];
        parameterIn1= parameterMatrix(i,1:6);
        Injection= parameterMatrix(i,7:9)+0.0001;
        parameterIn2=[parameterIn1(1:2);parameterIn1(3:4);parameterIn1(5:6)];
        % up
        %%%%%%%%%%%%%%%%%%%%%

        [t1, R_obs1, Tsol] = parameter_to_noisyData(nCompIn1,parameterIn2,Injection);
        t_new=linspace(2,Tsol/60,timeN);
        R_new = interp1(t1,R_obs1,t_new);
        OutputMatrix(i,:)=R_new;   %density curve of the ith sample
    end


xlswrite('N3-param9.xlsx',parameterMatrix)
xlswrite('N3-Out9.xlsx',OutputMatrix)





