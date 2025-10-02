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



% function [t, R_obs, Tsol] = parameter_to_noisyData(nCompIn, parameterIn, Injection)
%     % 调用求解器
%     [t, C, Tsol] = parameter_to_solutions(nCompIn, parameterIn, Injection);
% 
%     % ---- 防越界护栏 ----
%     nComp = size(C, 2);
%     if nComp ~= nCompIn
%         % 若维度不匹配，优先用实际返回的列数，避免索引越界
%         warning('parameter_to_noisyData:CompMismatch', ...
%             'nCompIn=%d，但解返回的列数=%d。请检查 parameterIn/Injection 的维度。此处按 nComp=%d 继续。', ...
%             nCompIn, nComp, nComp);
%     end
% 
%     % ---- 噪声参数（相对噪声幅度，可按需调）----
%     comp_noise   = 0.00;  % 每个组分的相对噪声，例如 0.1
%     model_noise  = 0.00;  % 总信号的相对噪声，例如 0.002
% 
%     % ---- 组件噪声：矢量化写法 ----
%     if comp_noise > 0
%         % 为每个组分生成独立高斯噪声，形状与 C 相同
%         C_noise = C .* (1 + comp_noise * randn(size(C)));
%     else
%         C_noise = C;
%     end
% 
%     % 负值截断到 0（物理上不允许负浓度）
%     C_noise = max(C_noise, 0);
% 
%     % ---- 合成总信号并加“检测器噪声” ----
%     R_obs = sum(C_noise, 2);
%     if model_noise > 0
%         R_obs = R_obs .* (1 + model_noise * randn(size(R_obs)));
%     end
% end
