function post_product_ff = deAliased_product_psuedo_spect(Nx, Nx_pssp, vel_hat_A, vel_hat_B)
    preproduct_vel_hat_if = vel_hat_A .* vel_hat_B; % 3N/2 grid size 
    
    ff_preproduct_vel_hat_if = fftn(preproduct_vel_hat_if);

    %%
    
    post_product_ff = zeros(Nx,Nx,Nx);

    %positive quadrant
    post_product_ff(1:Nx/2 , 1:Nx/2 , 1:Nx/2) = ff_preproduct_vel_hat_if(1:Nx/2 , 1:Nx/2 , 1:Nx/2);

    %one neg 2 positive quadrant
    post_product_ff(Nx/2 + 1:end, 1:Nx/2 , 1:Nx/2) = ff_preproduct_vel_hat_if(Nx+1:end , 1:Nx/2 , 1:Nx/2);
    post_product_ff(1:Nx/2, Nx/2 + 1:end, 1:Nx/2) = ff_preproduct_vel_hat_if(1:Nx/2 , Nx+1:end , 1:Nx/2);
    post_product_ff(1:Nx/2 , 1:Nx/2 , Nx/2 + 1:end) = ff_preproduct_vel_hat_if(1:Nx/2 , 1:Nx/2 , Nx+1:end);

    %two neg 1 positive quadrant
    post_product_ff(Nx/2 + 1:end , Nx/2 + 1:end , 1:Nx/2) = ff_preproduct_vel_hat_if(Nx+1:end , Nx+1:end , 1:Nx/2);
    post_product_ff(1:Nx/2 , Nx/2 + 1:end , Nx/2 + 1:end) = ff_preproduct_vel_hat_if(1:Nx/2 , Nx+1:end , Nx+1:end);
    post_product_ff(Nx/2 + 1:end , 1:Nx/2 , Nx/2 + 1:end) = ff_preproduct_vel_hat_if(Nx+1:end , 1:Nx/2 , Nx+1:end);

    %all negative
    post_product_ff(Nx/2 + 1:end , Nx/2 + 1:end , Nx/2 + 1:end) = ff_preproduct_vel_hat_if(Nx+1:end , Nx+1:end , Nx+1:end);
    
end