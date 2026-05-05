function [u_rhs , v_rhs , w_rhs] =main_function_code(Nx, uhatk_old, vhatk_old, whatk_old, kx, ky, kz, k_mag, diff_gamma)

    int_summ_conv_terms = 0;
    
    [P] = pressure_poisson(Nx, uhatk_old, vhatk_old, whatk_old, kx, ky, kz, k_mag);

    [u_pssp_if, v_pssp_if, w_pssp_if, Nx_pssp] = pre_product_psuedo_spect(Nx, uhatk_old, vhatk_old, whatk_old);
    
    for i = 1 : 3
    
        if i == 1
            vel_hat_A_pssp = u_pssp_if;
            vel_hat_A = uhatk_old;
            k_A = kx;
        elseif i == 2
            vel_hat_A_pssp = v_pssp_if;
            vel_hat_A = vhatk_old;
            k_A = ky;
        elseif i == 3
            vel_hat_A_pssp = w_pssp_if;
            vel_hat_A = whatk_old;
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

            int_conv = -1i .* k_B .*post_product_ff;

            int_summ_conv_terms = int_summ_conv_terms + int_conv;

        end

        rhs_1 = int_summ_conv_terms;

        rhs_2 = - 1i .* k_A .* P;

        rhs_3 = -diff_gamma .* (k_mag.^2) .* vel_hat_A;

        rhs = rhs_1 + rhs_2 + rhs_3;

        if i == 1
            u_rhs = rhs;
        elseif i == 2
            v_rhs = rhs;
        elseif i==3
            w_rhs = rhs;
        end

        int_summ_conv_terms = 0; 
    end

end