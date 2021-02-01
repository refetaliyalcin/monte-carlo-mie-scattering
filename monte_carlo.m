function [r_tot,t_tot,a_tot]=monte_carlo(photon_number,s_ref,cos_gecen,h,scat_prob,mu_tot,n_medium,k_medium,n_subs,k_subs,inv_cdf_cos,n_cdf_random)
r_tot=0;
t_tot=0;
a_tot=0;
for i=1:photon_number
    r_no=0;
    t_no=0;
    a_no=0;
    Rastgele=rand();
    if Rastgele>s_ref
        %gecer
        alive=1;
    else
        %yansir
        alive=0;
        r_no=1;
    end
    x=0;
    y=0;
    z=0;
    s_x=0;s_y=sqrt(1-cos_gecen*cos_gecen);s_z=cos_gecen; %direction vector
    l_beta=-log(rand())/mu_tot; %ext length
    while alive   
        if (s_z>0)
            l_w = (h - z)/s_z; %distance to lower boundary
        else
            l_w = -z/s_z; %distance to upper boundary
        end
        if l_w<l_beta
            min_index=1;
            min_l=l_w;
        else
            min_index=2;
            min_l=l_beta;
        end
    %         [min_l,min_index]=min([l1 l2]); %select minimum
        x=x+min_l*s_x;
        y=y+min_l*s_y;
        z=z+min_l*s_z;
        if (min_index==1)
    %            disp('hit boundary');
            alive=snell(s_z,n_medium,k_medium,n_subs,k_subs);
            if (alive==0)
                if s_z>0
                    t_no=1;
                else
                    r_no=1;
                end
            else
                l_beta=l_beta-l_w;
                s_z=-s_z;
            end
        else
            random_no=rand();
            if random_no<scat_prob
    %               disp('scattering');
                cos_theta= inv_cdf_cos(ceil(rand()*n_cdf_random));
                [s_x,s_y,s_z]=scatter_mc(cos_theta,s_x,s_y,s_z);
                l_beta=-log(rand())/mu_tot;
            else
    %               disp('absorption');
                alive=0;
                a_no=1;
            end
        end
    end
    r_tot=r_tot+r_no;
    t_tot=t_tot+t_no;
    a_tot=a_tot+a_no;
end
r_tot=r_tot/photon_number;
t_tot=t_tot/photon_number;
a_tot=a_tot/photon_number;