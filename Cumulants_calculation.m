function [c] = Cumulants_calculation(k,DoF,Noncentral_parameters,weight,size_SNR)
c_1 = zeros(1,size_SNR);
c_2 = zeros(1,size_SNR);

for i = 1:size_SNR
    for j = 1:length(DoF)
        c_1(i) = c_1(i) + weight(i,j)^k .* DoF(j);
    end
end

for i = 1:size_SNR
    for j = 1:length(DoF)
        c_1(i) = c_1(i) + weight(i,j)^k .* Noncentral_parameters(i,j);
    end
end

c = c_1 + k.*c_2;