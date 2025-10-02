function [t, R_obs, Tsol] = parameter_to_noisyData(nCompIn,parameterIn,Injection)

[t, C, Tsol] = parameter_to_solutions(nCompIn,parameterIn,Injection);

noise_level=ones(nCompIn,1);
for i=1:nCompIn
    noise_level(i)=0.0*nCompIn;
end

C_noise=0*C;
for i=1:nCompIn
    C_noise(:,i)=C(:,i).*(1+noise_level(i)*randn(size(C,1),1));
end

C_noise=(C_noise+abs(C_noise))/2;

noise_level_model=0.00;

R_obs=sum(C_noise,2).*(1+noise_level_model*randn(size(C,1),1));

