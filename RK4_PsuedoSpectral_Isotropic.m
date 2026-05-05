function [uhatk_np1_RK4 , vhatk_np1_RK4 , whatk_np1_RK4] = RK4_PsuedoSpectral_Isotropic(Nx, uhatk_old, vhatk_old, whatk_old, kx, ky, kz, k_mag,delta_t_max,diff_gamma)

    alpha_const = 1/2 ;
    
    [u_rhs_f1 , v_rhs_f1 , w_rhs_f1] = main_function_code(Nx, uhatk_old, vhatk_old, whatk_old, kx, ky, kz, k_mag, diff_gamma);
    u_kbar1 = delta_t_max * u_rhs_f1;
    v_kbar1 = delta_t_max * v_rhs_f1;
    w_kbar1 = delta_t_max * w_rhs_f1;

    [u_rhs_f2 , v_rhs_f2 , w_rhs_f2] = main_function_code(Nx, (uhatk_old + (alpha_const * u_kbar1)), (vhatk_old + (alpha_const * v_kbar1)), (whatk_old + (alpha_const * w_kbar1)), kx, ky, kz, k_mag, diff_gamma);
    u_kbar2 = delta_t_max * u_rhs_f2;
    v_kbar2 = delta_t_max * v_rhs_f2;
    w_kbar2 = delta_t_max * w_rhs_f2;
    
    [u_rhs_f3 , v_rhs_f3 , w_rhs_f3] = main_function_code(Nx, (uhatk_old + (alpha_const * u_kbar2)), (vhatk_old + (alpha_const * v_kbar2)), (whatk_old + (alpha_const * w_kbar2)), kx, ky, kz, k_mag,diff_gamma);
    u_kbar3 = delta_t_max * u_rhs_f3;
    v_kbar3 = delta_t_max * v_rhs_f3;
    w_kbar3 = delta_t_max * w_rhs_f3;
    
    [u_rhs_f4 , v_rhs_f4 , w_rhs_f4] = main_function_code(Nx, (uhatk_old + (alpha_const * u_kbar3)), (vhatk_old + (alpha_const * v_kbar3)), (whatk_old + (alpha_const * w_kbar3)), kx, ky, kz, k_mag,diff_gamma);
    u_kbar4 = delta_t_max * u_rhs_f4;
    v_kbar4 = delta_t_max * v_rhs_f4;
    w_kbar4 = delta_t_max * w_rhs_f4;

    uhatk_np1_RK4 = uhatk_old + (1/6)*u_kbar1 + ((1/3)*(u_kbar2 + u_kbar3)) + (1/6)*u_kbar4 ;
    vhatk_np1_RK4 = vhatk_old + (1/6)*v_kbar1 + ((1/3)*(v_kbar2 + v_kbar3)) + (1/6)*v_kbar4 ;
    whatk_np1_RK4 = whatk_old + (1/6)*w_kbar1 + ((1/3)*(w_kbar2 + w_kbar3)) + (1/6)*w_kbar4 ;
end