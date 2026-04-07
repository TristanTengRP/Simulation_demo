clear all;

SNR_dB=[-20:2:10]; 
simulations=3000;
thresh_gau=[1:1:20]; 
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

%% Simulation and theoretical analysis of scheme using Gaussian AoA
DOA_est_E=load("Eve_Gaussian_5000_52.mat","theta_est","theta_true");
DOA_est_A=load("Alice_Gaussian_3000.mat","theta_est","theta_true");


var_DOA_A = reshape(var(DOA_est_A.theta_est),l_snr,1);
var_DOA_E = reshape(var(DOA_est_E.theta_est(:,:,1:16)),l_snr,1);
mean_DOA_A=[55];
mean_DOA_E=[52];

var_DOA_c=(var_DOA_A.*var_DOA_E)./(var_DOA_E-var_DOA_A);
m_theta=var_DOA_c.*((mean_DOA_A./var_DOA_A)-(mean_DOA_E./var_DOA_E));

%% PD_simulation using Gaussian AoA
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

%% PD Theory using Gaussian AoA
PD_DOA_Theo = zeros(l_snr,l_thresh_gau);
lamda_1_theta = zeros(l_snr,1);

for i=1:l_snr
    for j=1:l_thresh_gau
        lamda_1_theta(i) = ((mean_DOA_E-m_theta(i))./sqrt(var_DOA_E(i))).'*((mean_DOA_E-m_theta(i))./sqrt(var_DOA_E(i))); 
        PD_DOA_Theo(i,j)=1-ncx2cdf(2.*var_DOA_c(i)./var_DOA_E(i).*thresh_gau(j),length(mean_DOA_A),lamda_1_theta(i));
    end
end

%% PF simulation using Gaussian AoA
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

%% PF Theory using Gaussian AoA
PF_DOA_Theo = zeros(l_snr,l_thresh_gau);
lamda_0_theta = zeros(l_snr,1);
for i=1:l_snr
    for j=1:l_thresh_gau
        lamda_0_theta(i) = ((mean_DOA_A-m_theta(i))./sqrt(var_DOA_A(i))).'*((mean_DOA_A-m_theta(i))./sqrt(var_DOA_A(i))); 
        PF_DOA_Theo(i,j)=1-ncx2cdf(2.*var_DOA_c(i)./var_DOA_A(i).*thresh_gau(j),length(mean_DOA_E),lamda_0_theta(i));
    end
end

%% Simulation and theoretical analysis of scheme using Laplace AoA
DOA_est_E_lap=load("Eve_Laplace_3000_1.mat","theta_est","theta_true"); 
DOA_est_A_lap=load("Alice_Laplace_3000.mat","theta_est","theta_true");

for i=1:l_snr
    [sigma_e,mu_e] = histfitlaplace(DOA_est_E_lap.theta_est(1,:,i).');
    [sigma_a,mu_a] = histfitlaplace(DOA_est_A_lap.theta_est(1,:,i).');
    sigma_E(i) = sigma_e;
    sigma_A(i) = sigma_a;
end

%% PF simulation
LRT_lap_DOA_PF=zeros(l_snr,simulations);

for i=1:l_snr
    for j=1:l_thresh_lap
        pf_count=0;
        for sim=1:simulations
            LRT_lap_DOA_PF(i,sim)=2./(sigma_A(i)).*sum(abs(DOA_est_A_lap.theta_est(:,sim,i)),1);
            if LRT_lap_DOA_PF(i,sim)>2./(sigma_A(i)).*thresh_lap(j) 
                pf_count=pf_count+1;
            end
        end
        PF_lap_DOA_sim(i,j)=pf_count/simulations;
    end
end

%% PF theoretical
PF_lap_DOA_Theo=zeros(l_snr,l_thresh_lap);
for i=1:l_snr
    for j=1:l_thresh_lap
        PF_lap_DOA_Theo(i,j)=Qchipr2(2.*length(DOA_est_A_lap.theta_est(:,1,1)),0,2./(sigma_A(i)).*thresh_lap(j),1e-5);
    end
end

%% PD simulation
LRT_lap_DOA_PD=zeros(l_snr,simulations);

for i=1:l_snr
    for j=1:length(thresh_lap)
        pd_count=0;
        for sim=1:simulations
            LRT_lap_DOA_PD(i,sim)=2./(sigma_E(i)).*sum(abs(DOA_est_E_lap.theta_est(:,sim,i)),1);
            if LRT_lap_DOA_PD(i,sim)>2./(sigma_E(i)).*thresh_lap(j) 
                pd_count=pd_count+1;
            end
        end
        PD_lap_DOA_sim(i,j)=pd_count/simulations;
    end
end

PD_lap_DOA_Theo=zeros(l_snr,l_thresh_lap);
for i=1:l_snr
    for j=1:l_thresh_lap
        PD_lap_DOA_Theo(i,j)=Qchipr2(2.*length(DOA_est_E_lap.theta_est(:,1,1)),0,2./(sigma_E(i)).*thresh_lap(j),1e-5);
    end
end

%% Simulation and theoretical analysis of scheme using MC
MC_est_E=load("Eve_Gaussian_3000_MC_1.mat","MC_est","MC_true"); 
MC_est_A=load("Alice_Gaussian_3000.mat","MC_est","MC_true");
P = length(MC_est_A.MC_est(:,1,1)); 
var_MC_A = 2.*(reshape(var(real(MC_est_A.MC_est(3,:,:))),l_snr,1)+reshape(var(real(MC_est_A.MC_est(3,:,:))),l_snr,1))./2;
var_MC_E = 2.*(reshape(var(real(MC_est_E.MC_est(3,:,:))),l_snr,1)+reshape(var(real(MC_est_E.MC_est(3,:,:))),l_snr,1))./2;
mean_MC_A = [0.167459095433972	0.0837295477169860];
mean_MC_E = [0.666666666666667	0.333333333333333];

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
                pd_MC_count=pd_MC_count+1;
            end
        end
        PD_MC_sim(i,j)=pd_MC_count/simulations;
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
lamda_2_mc = zeros(l_snr,1);

for i=1:l_snr
    for j=1:l_thresh_MC
        lamda_2_mc(i) = (sqrt(2).*(mean_MC_A-m_MC(i,:)).'./sqrt(var_MC_A(i)))'*(sqrt(2).*(mean_MC_A-m_MC(i,:)).'./sqrt(var_MC_A(i))); 
        PF_MC_Theo(i,j)=1-ncx2cdf(2.*var_MC_c(i)./var_MC_A(i).*thresh_MC(j),2*(P-1),lamda_2_mc(i));
    end
end


%% ROC curve for the Authentication scheme jointly using MC and Gaussian DOA
PF_sim = PF_DOA_sim.*PF_MC_sim + (1-PF_DOA_sim).*PF_MC_sim + (1-PF_MC_sim).*PF_DOA_sim;
PD_sim = 1-(1-PD_MC_sim).*(1-PD_DOA_sim);
PF_Theo = PF_DOA_Theo.*PF_MC_Theo + (1-PF_DOA_Theo).*PF_MC_Theo + (1-PF_MC_Theo).*PF_DOA_Theo;
PD_Theo = 1-(1-PD_MC_Theo).*(1-PD_DOA_Theo);

figure(1);
hold on;
grid on;
plot(PF_Theo(11,:),PD_Theo(11,:),'r-',LineWidth=1.5);% 0dB
plot(PF_Theo(13,:),PD_Theo(13,:),'b:',LineWidth=1.5);% 4dB
plot(PF_Theo(15,:),PD_Theo(15,:),'g--',LineWidth=1.5);% 8dB
plot(PF_sim(11,:),PD_sim(11,:),'ko',LineWidth=1.5);% 0dB
plot(PF_sim(13,:),PD_sim(13,:),'ko',LineWidth=1.5);% 4dB
plot(PF_sim(15,:),PD_sim(15,:),'ko',LineWidth=1.5);% 8dB
hold off;
xlim([0 0.05]);
ylim([0.99 1]);
xlabel("P_f");
ylabel("P_d");
legend("Theo: SNR = 0 dB", "Theo: SNR = 4 dB", "Theo: SNR = 8 dB","Simulation",'Location', 'southeast','EdgeColor', 'none');
 
