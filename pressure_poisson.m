function [P] = pressure_poisson(Nx, uhatk_old, vhatk_old, whatk_old, kx, ky, kz, k_mag)

[u_pssp_if, v_pssp_if, w_pssp_if, Nx_pssp] = pre_product_psuedo_spect(Nx, uhatk_old, vhatk_old, whatk_old);

int_summ_terms = 0;

ext_summ_terms = 0;

%%
k2 = k_mag.^2;
k2(k2 == 0) = 1;

%%
    for i = 1 : 3
        
        if i == 1
            vel_hat_A_pssp = u_pssp_if;
            k_A = kx;
        elseif i == 2
            vel_hat_A_pssp = v_pssp_if;
            k_A = ky;
        elseif i == 3
            vel_hat_A_pssp = w_pssp_if;
            k_A = kz;
        end
    
        for j = 1 : 3
            
            if j == 1
                vel_hat_B_pssp = u_pssp_if;
                k_B = kx;
            elseif j == 2 
                vel_hat_B_pssp = v_pssp_if;
                k_B = ky;
            elseif j == 3
                vel_hat_B_pssp = w_pssp_if;
                k_B = kz;
            end
    
            post_product_ff = deAliased_product_psuedo_spect(Nx, Nx_pssp, vel_hat_A_pssp, vel_hat_B_pssp);
    
            int_P = k_B .*post_product_ff;
    
            int_summ_terms = int_summ_terms + int_P;
            
        end
        
        ext_P = (k_A ./ k2) .* int_summ_terms ;
        
        ext_summ_terms = ext_summ_terms + ext_P;

        int_summ_terms = 0; 
    
    end

    P = - ext_summ_terms;

    P(k_mag == 0) = 0;
end