clc
clear all
close all

f_v=100*1e-6;
radius=1*10^-9;
thickness=500*10^-6;

lamda=(900:300:900)'*10^-9;

photon_number=10^5;
n_cdf_random=1000;
nang_gid=500;

lamda_nm=lamda*10^9;
% polar_angle=linspace(0,89.99999,3); %30 makul
polar_angle=0; %30 makul
polar_angle_rad=polar_angle*pi/180;

n_pigment=sio2_n(lamda); 
k_pigment=sio2_k(lamda); 
n_air=ones(length(lamda),1); 
k_air=zeros(length(lamda),1); 
n_medium=n_air; 
k_medium=k_air; 

n_substrat=ones(length(lamda),1); 
k_substrat=zeros(length(lamda),1); 


thickness_um=thickness*10^6;


teta_prime=zeros(length(lamda),length(polar_angle));
sur_reflection=zeros(length(lamda),length(polar_angle));
for i=1:length(lamda)
    for j=1:length(polar_angle_rad)
        teta_prime(i,j)=F_fresnel_2(n_medium(i),k_medium(i),polar_angle_rad(j))*180/pi;
        cos_teta=cosd(polar_angle(j));
        sin_teta=sqrt(1-cos_teta*cos_teta);
        carpan2=1/(n_medium(i)-1i*k_medium(i));
        sin_x2=sin_teta*carpan2;
        cos_x2=sqrt(1-sin_x2*sin_x2);
        carpan1=cos_teta/cos_x2;
        carpan3=cos_x2/cos_teta;
        E_parallel=(carpan1-carpan2)/(carpan1+carpan2);
        R_parallel=E_parallel*conj(E_parallel);
        E_orth=(carpan3-carpan2)/(carpan3+carpan2);
        R_orth=E_orth*conj(E_orth);
        reflectance=real(R_parallel+R_orth)*0.5;
        sur_reflection(i,j)=reflectance;
    end
end

ref_lamda=zeros(length(lamda),length(polar_angle));
tra_lamda=zeros(length(lamda),length(polar_angle));
abs_lamda=zeros(length(lamda),length(polar_angle));
inv_cdf_cos=zeros(n_cdf_random,length(lamda));
HG_invcdf=zeros(n_cdf_random,1);
Qsca=zeros(length(lamda),1);
Qabs=zeros(length(lamda),1);
mu_tot=zeros(length(lamda),1);
scat_prob=zeros(length(lamda),1);
g=zeros(length(lamda),1);

V=4*pi*radius^3/3;
x2=cos(transpose(linspace(0,pi,nang_gid)));
cdf=zeros(nang_gid,1);
zero_to_one=linspace(0,1,n_cdf_random)';
zero_to_pi=linspace(0,pi,length(cdf))';
tic
for i=1:length(lamda)
    x=2*pi*n_medium(i)*radius/lamda(i);
    m=(n_pigment(i)+1i*k_pigment(i))/n_medium(i);
    fonksiyon=Mie(m,x);
    Qsca(i)=fonksiyon(2);
    Qabs(i)=fonksiyon(3);
    Csca=pi*radius^2*Qsca(i);
    Cabs=pi*radius^2*Qabs(i);
    alfa=f_v*Csca/V;
    beta=f_v*Cabs/V+(1-f_v)*4*pi*k_medium(i)/lamda(i);
    mu_tot(i)=alfa+beta;
    scat_prob(i)=alfa/mu_tot(i);
    g(i)=fonksiyon(5);
    scat_ph_fn=Mie_phasefn(m, x, nang_gid)*pi/(x*x*Qsca(i));
    
    %figure
    %polarplot(teta,scat_ph_fn)

    y2=scat_ph_fn./trapz(x2,scat_ph_fn);
    for i2=2:nang_gid
        cdf(i2)=trapz(x2(1:i2),y2(1:i2));
    end  
    inv_cdf_cos(:,i)=cos(interp1(cdf,zero_to_pi,zero_to_one,'linear','extrap'));
    
    for i3=1:n_cdf_random
        rnd=zero_to_one(i3);
        %carpan=(1 - g(i)*g(i))/(1 - g(i) + 2*g(i)*rnd);
        %HG_invcdf(n_cdf_random-i3+1)=(1 + g(i)*g(i) - carpan*carpan)/(2*g(i));
        r=4*rnd-2;
        HG_invcdf(n_cdf_random-i3+1)=(r+sqrt(1 + r*r))^(1/3) - (r + sqrt(1 + r*r))^(-1/3) ;
    end
    
    figure
    plot(zero_to_one,inv_cdf_cos,zero_to_one,HG_invcdf)
    legend('Mie','Rayleigh')
    
end    

% for i=1:length(lamda)
%     for j=1:length(polar_angle)
%         [ref_lamda(i,j),tra_lamda(i,j),abs_lamda(i,j)]=monte_carlo(photon_number,sur_reflection(i,j),cosd(teta_prime(i,j)),thickness,scat_prob(i),mu_tot(i),n_medium(i),k_medium(i),n_substrat(i),k_substrat(i),inv_cdf_cos(:,i),n_cdf_random);
%     end
% end
% toc

% figure
% plot(lamda_nm,ref_lamda(:,1),lamda_nm,tra_lamda(:,1),lamda_nm,abs_lamda(:,1),'LineWidth',2)
% ylim([0 1])
