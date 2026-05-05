function [u_pssp_if, v_pssp_if, w_pssp_if, Nx_pssp] = pre_product_psuedo_spect(Nx, uhatk, vhatk, whatk)

%%
%Psuedo-Spectral starting

Nx_pssp = (3/2) * Nx;

uhatk_pssp = zeros(Nx_pssp, Nx_pssp, Nx_pssp);
vhatk_pssp = uhatk_pssp;
whatk_pssp = uhatk_pssp;


%%
%zero padding for all 8 quadrants in 3D space.
% positive qudrant
uhatk_pssp(1:Nx/2 , 1:Nx/2 , 1:Nx/2) = uhatk(1:Nx/2 , 1:Nx/2 , 1:Nx/2);

%one neg 2 positive quadrant
uhatk_pssp(Nx+1:end , 1:Nx/2 , 1:Nx/2) = uhatk(Nx/2 + 1:end, 1:Nx/2 , 1:Nx/2);
uhatk_pssp(1:Nx/2 , Nx+1:end , 1:Nx/2) = uhatk(1:Nx/2 , Nx/2 + 1:end , 1:Nx/2);
uhatk_pssp(1:Nx/2 , 1:Nx/2 , Nx+1:end) = uhatk(1:Nx/2 , 1:Nx/2 , Nx/2 + 1:end);

%2 neg one pos
uhatk_pssp(Nx+1:end , Nx+1:end , 1:Nx/2) = uhatk(Nx/2 + 1:end, Nx/2 + 1:end , 1:Nx/2);
uhatk_pssp(Nx+1:end , 1:Nx/2 , Nx+1:end) = uhatk(Nx/2 + 1:end , 1:Nx/2 , Nx/2 + 1:end);
uhatk_pssp(1:Nx/2 , Nx+1:end , Nx+1:end) = uhatk(1:Nx/2 , Nx/2 + 1:end , Nx/2 + 1:end);

%all negative
uhatk_pssp(Nx+1:end , Nx+1:end , Nx+1:end) = uhatk(Nx/2 + 1:end , Nx/2 + 1:end , Nx/2 + 1:end);

%%

vhatk_pssp(1:Nx/2 , 1:Nx/2 , 1:Nx/2) = vhatk(1:Nx/2 , 1:Nx/2 , 1:Nx/2);

%one neg 2 positive quadrant
vhatk_pssp(Nx+1:end , 1:Nx/2 , 1:Nx/2) = vhatk(Nx/2 + 1:end, 1:Nx/2 , 1:Nx/2);
vhatk_pssp(1:Nx/2 , Nx+1:end , 1:Nx/2) = vhatk(1:Nx/2 , Nx/2 + 1:end , 1:Nx/2);
vhatk_pssp(1:Nx/2 , 1:Nx/2 , Nx+1:end) = vhatk(1:Nx/2 , 1:Nx/2 , Nx/2 + 1:end);

%2 neg one pos
vhatk_pssp(Nx+1:end , Nx+1:end , 1:Nx/2) = vhatk(Nx/2 + 1:end, Nx/2 + 1:end , 1:Nx/2);
vhatk_pssp(Nx+1:end , 1:Nx/2 , Nx+1:end) = vhatk(Nx/2 + 1:end , 1:Nx/2 , Nx/2 + 1:end);
vhatk_pssp(1:Nx/2 , Nx+1:end , Nx+1:end) = vhatk(1:Nx/2 , Nx/2 + 1:end , Nx/2 + 1:end);

%all negative
vhatk_pssp(Nx+1:end , Nx+1:end , Nx+1:end) = vhatk(Nx/2 + 1:end , Nx/2 + 1:end , Nx/2 + 1:end);

%%
whatk_pssp(1:Nx/2 , 1:Nx/2 , 1:Nx/2) = whatk(1:Nx/2 , 1:Nx/2 , 1:Nx/2);

%one neg 2 positive quadrant
whatk_pssp(Nx+1:end , 1:Nx/2 , 1:Nx/2) = whatk(Nx/2 + 1:end, 1:Nx/2 , 1:Nx/2);
whatk_pssp(1:Nx/2 , Nx+1:end , 1:Nx/2) = whatk(1:Nx/2 , Nx/2 + 1:end , 1:Nx/2);
whatk_pssp(1:Nx/2 , 1:Nx/2 , Nx+1:end) = whatk(1:Nx/2 , 1:Nx/2 , Nx/2 + 1:end);

%2 neg one pos
whatk_pssp(Nx+1:end , Nx+1:end , 1:Nx/2) = whatk(Nx/2 + 1:end, Nx/2 + 1:end , 1:Nx/2);
whatk_pssp(Nx+1:end , 1:Nx/2 , Nx+1:end) = whatk(Nx/2 + 1:end , 1:Nx/2 , Nx/2 + 1:end);
whatk_pssp(1:Nx/2 , Nx+1:end , Nx+1:end) = whatk(1:Nx/2 , Nx/2 + 1:end , Nx/2 + 1:end);

%all negative
whatk_pssp(Nx+1:end , Nx+1:end , Nx+1:end) = whatk(Nx/2 + 1:end , Nx/2 + 1:end , Nx/2 + 1:end);


%%
u_pssp_if = ifftn(uhatk_pssp);
v_pssp_if = ifftn(vhatk_pssp);
w_pssp_if = ifftn(whatk_pssp);
end