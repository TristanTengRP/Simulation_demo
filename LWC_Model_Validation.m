clear all;

SNR_dB=[-20:2:10]; 
simulations=3000;
thresh_gau=[1:1:20]; % thresh of authentication 5:0.2:15 for ROC curve   
thresh_lap=[0:0.1:2];
thresh_MC=[1:1:20];
l_snr=length(SNR_dB);
l_thresh_gau=length(thresh_gau);
l_thresh_MC=length(thresh_MC);
l_thresh_lap=length(thresh_lap);
PF_DOA_sim=zeros(l_snr,l_thresh_gau);
PD_DOA_sim=zeros(l_snr,l_thresh_gau);
PF_lap_DOA_sim=zeros(l_snr,l_thresh_lap);
PD_lap_DOA_sim=zeros(l_snr,l_thresh_lap);
PF_MC_sim=zeros(l_snr,l_thresh_MC);
PD_MC_sim=zeros(l_snr,l_thresh_MC);

%% Simulation and theoretical analysis of scheme using Gaussian DOA
DOA_est_E=load("Eve_Gaussian_5000_52.mat","theta_est","theta_true"); 
DOA_est_A=load("Alice_Gaussian_3000.mat","theta_est","theta_true");

var_DOA_A = reshape(var(DOA_est_A.theta_est),l_snr,1);
var_DOA_E = reshape(var(DOA_est_E.theta_est(:,:,1:16)),l_snr,1);
mean_DOA_A=[55];
mean_DOA_E=[52];

var_DOA_c=(var_DOA_A.*var_DOA_E)./(var_DOA_E-var_DOA_A);
m_theta=var_DOA_c.*((mean_DOA_A./var_DOA_A)-(mean_DOA_E./var_DOA_E));

%% PD_simulation using Gaussian DOA
LRT_DOA_1=zeros(l_snr,simulations);
LRT_DOA_PD=zeros(l_snr,simulations);
Theta_delta_theata_pd=zeros(l_snr,simulations);

for i=1:l_snr
    for j=1:l_thresh_gau
        pd_count=0;
        for sim=1:simulations
            Theta_delta_theata_pd(i,sim)=(DOA_est_E.theta_est(:,sim,i)-m_theta(i).')./sqrt(var_DOA_c(i));
            LRT_DOA_1(i,sim)=0.5.*(Theta_delta_theata_pd(i,sim).'*Theta_delta_theata_pd(i,sim)); 
            LRT_DOA_PD(i,sim)=(2.*var_DOA_c(i)./var_DOA_E(i).*LRT_DOA_1(i,sim));
            if LRT_DOA_PD(i,sim)>2.*var_DOA_c(i)./var_DOA_E(i).*thresh_gau(j) 
                pd_count=pd_count+1;
            end
        end
        PD_DOA_sim(i,j)=pd_count/simulations;
    end
end

%% PD Theory using Gaussian DOA
PD_DOA_Theo = zeros(l_snr,l_thresh_gau);
lamda_1_theta = zeros(l_snr,1);

for i=1:l_snr
    for j=1:l_thresh_gau
        lamda_1_theta(i) = ((mean_DOA_E-m_theta(i))./sqrt(var_DOA_E(i))).'*((mean_DOA_E-m_theta(i))./sqrt(var_DOA_E(i))); % H1条件下的中心尺度参数
        PD_DOA_Theo(i,j)=1-ncx2cdf(2.*var_DOA_c(i)./var_DOA_E(i).*thresh_gau(j),length(mean_DOA_A),lamda_1_theta(i));
    end
end

%% PF simulation using Gaussian DOA
LRT_DOA_2=zeros(l_snr,simulations);
LRT_DOA_PF=zeros(l_snr,simulations);
Theta_delta_theata_pf=zeros(l_snr,simulations);

for i=1:l_snr
    for j=1:l_thresh_gau
        pf_count=0;
        for sim=1:simulations
            Theta_delta_theata_pf(i,sim)=(DOA_est_A.theta_est(:,sim,i)-m_theta(i).')./sqrt(var_DOA_c(i));
            LRT_DOA_2(i,sim)=0.5.*(Theta_delta_theata_pf(i,sim).'*Theta_delta_theata_pf(i,sim));
            LRT_DOA_PF(i,sim)=(2.*var_DOA_c(i)./var_DOA_A(i).*LRT_DOA_2(i,sim));
            if LRT_DOA_PF(i,sim)>2.*var_DOA_c(i)./var_DOA_A(i).*thresh_gau(j)
                pf_count=pf_count+1;
            end
        end
        PF_DOA_sim(i,j)=pf_count/simulations;
    end
end

%% PF Theory using Gaussian DOA
PF_DOA_Theo = zeros(l_snr,l_thresh_gau);
lamda_0_theta = zeros(l_snr,1);
for i=1:l_snr
    for j=1:l_thresh_gau
        lamda_0_theta(i) = ((mean_DOA_A-m_theta(i))./sqrt(var_DOA_A(i))).'*((mean_DOA_A-m_theta(i))./sqrt(var_DOA_A(i))); % H1条件下的中心尺度参数
        PF_DOA_Theo(i,j)=1-ncx2cdf(2.*var_DOA_c(i)./var_DOA_A(i).*thresh_gau(j),length(mean_DOA_E),lamda_0_theta(i));
    end
end

%% Simulation and theoretical analysis of scheme using MC
MC_est_E=load("Eve_Gaussian_3000_MC_1.mat","MC_est","MC_true"); %load estimation of DOA under SNR=[0,30]dB
MC_est_A=load("Alice_Gaussian_3000.mat","MC_est","MC_true");
P = length(MC_est_A.MC_est(:,1,1)); % calculate the number of MC coefficients
var_MC_A = 2.*(reshape(var(real(MC_est_A.MC_est(3,:,:))),l_snr,1)+reshape(var(real(MC_est_A.MC_est(3,:,:))),l_snr,1))./2;% The average MC estimation variance of Alice;
var_MC_E = 2.*(reshape(var(real(MC_est_E.MC_est(3,:,:))),l_snr,1)+reshape(var(real(MC_est_E.MC_est(3,:,:))),l_snr,1))./2;% The average MC estimation variance of Eve;
mean_MC_A = [0.167459095433972	0.0837295477169860];% mean of MC, Alice
mean_MC_E = [0.666666666666667	0.333333333333333];% mean of MC, Eve 

var_MC_c = (var_MC_A.*var_MC_E)./(var_MC_E-var_MC_A);
m_MC = var_MC_c.*((mean_MC_A./var_MC_A)-(mean_MC_E./var_MC_E));

%% PD simulation using MC
LRT_MC_1=zeros(l_snr,simulations);
LRT_MC_PD=zeros(l_snr,simulations);

for i=1:l_snr
    for j=1:l_thresh_MC
        pd_MC_count=0;
        for sim=1:simulations
            C_delta_MC_pd=(MC_est_E.MC_est(2:P,sim,i)-m_MC(i,:).')./sqrt(var_MC_c(i));
            LRT_MC_1(i,sim)=(C_delta_MC_pd'*C_delta_MC_pd); 
            LRT_MC_PD(i,sim)=(2.*var_MC_c(i)./var_MC_E(i).*LRT_MC_1(i,sim));
            if LRT_MC_PD(i,sim)>2.*var_MC_c(i)./var_MC_E(i).*thresh_MC(j) 
                pd_MC_count = pd_MC_count+1;
            end
        end
        PD_MC_sim(i,j) = pd_MC_count/simulations;
    end
end

%% PD Theory using MC
PD_MC_Theo = zeros(l_snr,l_thresh_MC);
lamda_1_mc = zeros(l_snr,1);

for i=1:l_snr
    for j=1:l_thresh_MC
        lamda_1_mc(i) = (sqrt(2).*(mean_MC_E-m_MC(i,:)).'./sqrt(var_MC_E(i)))'*(sqrt(2).*(mean_MC_E-m_MC(i,:)).'./sqrt(var_MC_E(i))); 
        PD_MC_Theo(i,j)=1-ncx2cdf(2.*var_MC_c(i)./var_MC_E(i).*thresh_MC(j),2*(P-1),lamda_1_mc(i));
    end
end

%% PF simulation using MC
LRT_MC_2=zeros(l_snr,simulations);
LRT_MC_PF=zeros(l_snr,simulations);

for i=1:l_snr
    for j=1:l_thresh_MC
        pf_MC_count=0;
        for sim=1:simulations
            C_delta_MC_pf=(MC_est_A.MC_est(2:P,sim,i)-m_MC(i,:).')./sqrt(var_MC_c(i));
            LRT_MC_2(i,sim)=(C_delta_MC_pf'*C_delta_MC_pf); 
            LRT_MC_PF(i,sim)=(2.*var_MC_c(i)./var_MC_A(i).*LRT_MC_2(i,sim));
            if LRT_MC_PF(i,sim)>2.*var_MC_c(i)./var_MC_A(i).*thresh_MC(j) 
                pf_MC_count=pf_MC_count+1;
            end
        end
        PF_MC_sim(i,j)=pf_MC_count/simulations;
    end
end

%% PF Theory using MC
PF_MC_Theo = zeros(l_snr,l_thresh_MC);
lamda_0_mc = zeros(l_snr,1);

for i=1:l_snr
    for j=1:l_thresh_MC
        lamda_0_mc(i) = (sqrt(2).*(mean_MC_A-m_MC(i,:)).'./sqrt(var_MC_A(i)))'*(sqrt(2).*(mean_MC_A-m_MC(i,:)).'./sqrt(var_MC_A(i))); 
        PF_MC_Theo(i,j)=1-ncx2cdf(2.*var_MC_c(i)./var_MC_A(i).*thresh_MC(j),2*(P-1),lamda_0_mc(i));
    end
end


%% Simulation and theoretical analysis of linear weighted combination (lwc) method 
weight = 0.5; % initial weight
thresh_lwc = [0.1:0.5:20]; % threshold of lwc method
l_thresh_lwc = length(thresh_lwc);

%% PD simulation using lwc
LRT_lwc_PD = zeros(l_snr,simulations);
PD_LWC_sim = zeros(l_snr,l_thresh_lwc);
t_1 = zeros(l_snr,l_thresh_lwc);

for i=1:l_snr
    for j=1:l_thresh_lwc
        pd_lmc_count=0;
        for sim=1:simulations
            LRT_lwc_PD(i,sim) = (weight./var_MC_E(i)).* LRT_DOA_PD(i,sim) + ((1-weight)./var_DOA_E(i)).* LRT_MC_PD(i,sim);
            t_1(i,j) = thresh_lwc(j).*(1./(var_MC_E(i).*var_DOA_E(i)));
            if LRT_lwc_PD(i,sim) > t_1(i,j) 
                pd_lmc_count=pd_lmc_count+1;
            end
        end
        PD_LWC_sim(i,j)=pd_lmc_count/simulations;
    end
end

%% PD theory using lwc
[c_pd_1] = Cumulants_calculation(1,[length(mean_DOA_E),2*(P-1)],[lamda_1_theta,lamda_1_mc],[weight./var_MC_E,(1-weight)./var_DOA_E],l_snr);
[c_pd_2] = Cumulants_calculation(2,[length(mean_DOA_E),2*(P-1)],[lamda_1_theta,lamda_1_mc],[weight./var_MC_E,(1-weight)./var_DOA_E],l_snr);
[c_pd_3] = Cumulants_calculation(3,[length(mean_DOA_E),2*(P-1)],[lamda_1_theta,lamda_1_mc],[weight./var_MC_E,(1-weight)./var_DOA_E],l_snr);
[c_pd_4] = Cumulants_calculation(4,[length(mean_DOA_E),2*(P-1)],[lamda_1_theta,lamda_1_mc],[weight./var_MC_E,(1-weight)./var_DOA_E],l_snr);

s_pd_1 = c_pd_3./(c_pd_2.^(3/2));
s_pd_2 = c_pd_4./(c_pd_2.^2);
flag_pd = (s_pd_1.^2-s_pd_2);
if flag_pd(11) > 0  
    a_pd = 1./(s_pd_1-(s_pd_1.^2-s_pd_2).^(1/2));
    Non_central_pd = s_pd_1.*(a_pd.^3) - a_pd.^2;
    DoF_pd = a_pd.^2-2.*Non_central_pd;
else
    a_pd = 1./s_pd_1;
    Non_central_pd = s_pd_1.*(a_pd).^3-a_pd.^2;
    DoF_pd = a_pd.^2 -2.*Non_central_pd;
end

mu_chi_pd = Non_central_pd + DoF_pd;
sigma_chi_pd = sqrt(2).* a_pd;
mu_Q_pd = c_pd_1;
sigma_Q_pd = (2.*c_pd_2).^(1/2);

% PD theory formula
PD_LMC_Theo = zeros(l_snr,l_thresh_lwc);
t_1_final = zeros(l_snr,l_thresh_lwc);

t_1_tmp = (t_1 - mu_Q_pd');
for i=1:l_snr
    t_1_final(i,:) = (t_1_tmp(i,:)./sigma_Q_pd(i)).*sigma_chi_pd(i) + mu_chi_pd(i);
end

for i=1:l_snr
    for j=1:l_thresh_lwc
        PD_LMC_Theo(i,j)=1-ncx2cdf(t_1_final(i,j),DoF_pd(i),Non_central_pd(i));
    end
end

%% PF simulation using lwc
LRT_lwc_PF = zeros(l_snr,simulations);
PF_LWC_sim = zeros(l_snr,l_thresh_lwc);
t_0 = zeros(l_snr,l_thresh_lwc);

for i=1:l_snr
    for j=1:l_thresh_lwc
        pf_lmc_count=0;
        for sim=1:simulations
            LRT_lwc_PF(i,sim) = (weight./var_MC_A(i)).* LRT_DOA_PF(i,sim) + ((1-weight)./var_DOA_A(i)).* LRT_MC_PF(i,sim);
            t_0(i,j) = thresh_lwc(j).*(1./(var_MC_A(i).*var_DOA_A(i)));
            if LRT_lwc_PF(i,sim) > t_0(i,j)
                pf_lmc_count=pf_lmc_count+1;
            end
        end
        PF_LWC_sim(i,j)=pf_lmc_count/simulations;
    end
end

%% PF theory
% Cumulants_calculation
[c_pf_1] = Cumulants_calculation(1,[length(mean_DOA_A),2*(P-1)],[lamda_0_theta,lamda_0_mc],[weight./var_MC_A,(1-weight)./var_DOA_A],l_snr);
[c_pf_2] = Cumulants_calculation(2,[length(mean_DOA_A),2*(P-1)],[lamda_0_theta,lamda_0_mc],[weight./var_MC_A,(1-weight)./var_DOA_A],l_snr);
[c_pf_3] = Cumulants_calculation(3,[length(mean_DOA_A),2*(P-1)],[lamda_0_theta,lamda_0_mc],[weight./var_MC_A,(1-weight)./var_DOA_A],l_snr);
[c_pf_4] = Cumulants_calculation(4,[length(mean_DOA_A),2*(P-1)],[lamda_0_theta,lamda_0_mc],[weight./var_MC_A,(1-weight)./var_DOA_A],l_snr);
s_pf_1 = c_pf_3./(c_pf_2.^(3/2));
s_pf_2 = c_pf_4./(c_pf_2.^2);
flag_pf = (s_pf_1.^2-s_pf_2);
if flag_pf(11) > 0  
    a_pf = 1./(s_pf_1-(s_pf_1.^2-s_pf_2).^(1/2));
    Non_central_pf = s_pf_1.*(a_pf.^3) - a_pf.^2;
    DoF_pf = a_pf.^2-2.*Non_central_pf;
else
    a_pf = 1./s_pf_1;
    Non_central_pf = s_pf_1.*(a_pf).^3-a_pf.^2;
    DoF_pf = a_pf.^2 -2.*Non_central_pf;
end

mu_chi_pf = Non_central_pf + DoF_pf;
sigma_chi_pf = sqrt(2).* a_pf;
mu_Q_pf = c_pf_1;
sigma_Q_pf = (2.*c_pf_2).^(1/2);

% PF theory formula
PF_LMC_Theo = zeros(l_snr,l_thresh_lwc);
t_0_final = zeros(l_snr,l_thresh_lwc);

t_0_tmp = (t_0 - mu_Q_pf');
for i=1:l_snr
    t_0_final(i,:) = (t_0_tmp(i,:)./sigma_Q_pf(i)).*sigma_chi_pf(i) + mu_chi_pf(i);
end


for i=1:l_snr
    for j=1:l_thresh_lwc
        PF_LMC_Theo(i,j)=1-ncx2cdf(t_0_final(i,j),DoF_pf(i),Non_central_pf(i));
    end
end

%% PF and PM for scheme using LMC.
figure(3);
hold on;
grid on;
plot(thresh_lwc,PF_LMC_Theo(11,:),'-.',LineWidth=1.5,Color=[0.93 0.69 0.13]);
plot(thresh_lwc,PF_LWC_sim(11,:),'s',LineWidth=1.5,Color=[0.93 0.69 0.13],MarkerSize=6);
ylabel("PF");
yyaxis right;
plot(thresh_lwc,1-PD_LMC_Theo(11,:),'-',LineWidth=1.5,Color=[0 0.4470 0.7410]);
plot(thresh_lwc,1-PD_LWC_sim(11,:),'o',LineWidth=1.5,Color=[0 0.4470 0.7410],MarkerSize=6);
ylabel('PM',Color='k');
hold off;
xlabel("\tau");
legend("SNR = 0 dB, Approximate PF (theory)","SNR = 0 dB, Simulations",...
    "SNR = 0 dB, Approximate PM (theory)","SNR = 0 dB, Simulations");

