clc
clear all
close all

%solution of radiative transfer equation in one layer plane parallel medium
%ray is incident from air to coating. coating is coated on a substrate
%substrate could be air or other material such as silver, glass etc.
%the code estimates spectral hemispherical reflectance, transmittance and absorptance
%can handle; independent scattering, boundary reflections, absorption in
%medium. can't handle; coherent backscattering, dependent scattering and polarized ray tracing.
%while calculating the scattering direction the code uses cumulative inverse
%relation of exact single scattering pahse function or henyey greenstein 
%function approximation depending on choice 

lambda=(300:20:800)*10^-9; %freespace wavelength of incident ray in meter 
thickness=100*10^-6;  %thickness of coating in meter 
radius=120*10^-9; %radius of particle in meter 
f_v=0.01; %volume fraction. 0.01 corresponds to 1% 
polar_angle=linspace(0,89.99999,9); %incident angles. 0 = perpendicular to slab face. 90 parallel and should be avoided.
use_HG=1; %if 0 use exact scattering phase function, if 1 uses henyey greenstein phase function approximation


n_pigment=sio2_n(lambda);  %real refractive index of pigment
k_pigment=sio2_k(lambda);  %imaginary refractive index of pigment
n_medium=ones(length(lambda),1); %real refractive index of medium
k_medium=zeros(length(lambda),1); %imaginary refractive index of medium
n_substrat=ones(length(lambda),1); %real refractive index of substrate
k_substrat=zeros(length(lambda),1); %imaginary refractive index of substrate

photon_number=10^5; %number of rays that will be traced, higher the number more accurate the result
n_cdf_random=1000; %how many pieces will be between (0,1) in random number relation with inverse cumulative function, higher the number accurate the phase function. useless if HG is used
nang_gid=1000; %how many pieces will be between (0,pi) in cumulative distribution function, higher the number more accurate the result. useless if HG is used


lambda_nm=lambda*10^9; %for plot
polar_angle_rad=polar_angle*pi/180;

% calculate surface reflection from air to medium. 
% medium to air and medium to substrate is calculated within snell.m.
% air to medium is calculated seperately since we need refraction angle
teta_prime=zeros(length(lambda),length(polar_angle));
sur_reflection=zeros(length(lambda),length(polar_angle));
for i=1:length(lambda)
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

%initialize variables
ref_lambda=zeros(length(lambda),length(polar_angle));
tra_lambda=zeros(length(lambda),length(polar_angle));
abs_lambda=zeros(length(lambda),length(polar_angle));
Qsca=zeros(length(lambda),1);
Qabs=zeros(length(lambda),1);
mu_tot=zeros(length(lambda),1);
scat_prob=zeros(length(lambda),1);
g=zeros(length(lambda),1);
if use_HG==0
    inv_cdf_cos=zeros(n_cdf_random,length(lambda));
    x2=flip(cos(transpose(linspace(0,pi,nang_gid))));
    zero_to_one=linspace(0,1,n_cdf_random)';
    zero_to_pi=linspace(0,pi,nang_gid)';
else
    inv_cdf_cos=zeros(1,length(lambda));
end
V=4*pi*radius^3/3; %volume of a single sphere
parfor i=1:length(lambda)
    x=2*pi*n_medium(i)*radius/lambda(i); %size parameter
    m=(n_pigment(i)+1i*k_pigment(i))/n_medium(i);
    fonksiyon=Mie(m,x); %calculate Lorenz Mie theory
    Qsca(i)=fonksiyon(2);%scattering efficiency
    Qabs(i)=fonksiyon(3);%absorption efficiency
    g(i)=fonksiyon(5);
    Csca=pi*radius^2*Qsca(i);%scattering crossection[m^2]
    Cabs=pi*radius^2*Qabs(i);%absorption crossection[m^2]
    alfa=f_v*Csca/V;%scattering coefficient[1/m]
    beta=f_v*Cabs/V+(1-f_v)*4*pi*k_medium(i)/lambda(i);%absorption coefficient[1/m] absorption of medium is implemented as a bulk property
    mu_tot(i)=alfa+beta;%extinction coefficient[1/m]
    scat_prob(i)=alfa/mu_tot(i);%scattering albedo also scattering probability
    if use_HG==0
        cdf=zeros(nang_gid,1);
        scat_ph_fn=Mie_phasefn(m, x, nang_gid)/(x*x*Qsca(i));%scattering phase function
        for i2=2:nang_gid
            cdf(i2)=trapz(x2(1:i2),scat_ph_fn(1:i2));%calculate cumulative distribution function
        end  
        inv_cdf_cos(:,i)=cos(interp1(cdf,zero_to_pi,zero_to_one,'linear','extrap'));%calculate inverse cumulative distribution function
    end
    
end   
tic
%loop the monte carlo code for lambda and polar_angle
for j=1:length(polar_angle)
    parfor i=1:length(lambda)
        [ref_lambda(i,j),tra_lambda(i,j),abs_lambda(i,j)]=monte_carlo(photon_number,sur_reflection(i,j),cosd(teta_prime(i,j)),thickness,scat_prob(i),mu_tot(i),n_medium(i),k_medium(i),n_substrat(i),k_substrat(i),inv_cdf_cos(:,i),n_cdf_random,g(i),use_HG);
    end
end
toc

figure %draw normal to diffuse R, T and A for normal incidence (first index in my case)
plot(lambda_nm,ref_lambda(:,1),lambda_nm,tra_lambda(:,1),lambda_nm,abs_lambda(:,1),'LineWidth',2)
ylim([0 1])
xlim([min(lambda_nm) max(lambda_nm)])
legend('Reflectance','Transmittance','Absorptance','Location', 'Best')
xlabel('Wavelength [nm]')
ylabel({'Normal to Hemispherical';'Reflectance, Transmittance, Absorptance'})
