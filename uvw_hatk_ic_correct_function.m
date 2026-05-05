function [uhatk, vhatk, whatk, u_ifft_sc, v_ifft_sc, w_ifft_sc, COM, q_sq_rms_scaled , max_COM, max_rel_err , E] = uvw_hatk_ic_correct_function(Nx, l, m, n, ll, mm, nn, kx, ky, kz, u0, k0, k_mag)

E = 16 * sqrt(2/pi) .* ((u0^2)/k0) .* (k_mag./k0).^4 .* exp(-2 .* (k_mag./k0).^2);

%%

theta1 = 2 * pi * rand(size(k_mag));
theta2 = 2 * pi * rand(size(k_mag));
phi = 2 * pi * rand(size(k_mag));

%%

uhatk = zeros(Nx,Nx,Nx);
vhatk = zeros(Nx,Nx,Nx);
whatk = zeros(Nx,Nx,Nx);

alpha = zeros(size(k_mag));
beta  = zeros(size(k_mag));

k_mag_nz_1 = k_mag > 0; %by using boolean we exclude 0

A = sqrt( E(k_mag_nz_1) ./ (4*pi*(k_mag(k_mag_nz_1)).^2) );

%%

alpha(k_mag_nz_1) = A .* exp(1i*theta1(k_mag_nz_1)) .* cos(phi(k_mag_nz_1));
beta(k_mag_nz_1)  = A .* exp(1i*theta2(k_mag_nz_1)) .* sin(phi(k_mag_nz_1));

%%
uhatk_nz_0 = (k_mag > 0) & (kx == 0) & (ky == 0);

uhatk(uhatk_nz_0) = alpha(uhatk_nz_0);
vhatk(uhatk_nz_0) = beta(uhatk_nz_0);
whatk(uhatk_nz_0) = 0;

lhs = (abs(uhatk(uhatk_nz_0)).^2 + abs(vhatk(uhatk_nz_0)).^2).* 4 .* pi .* (k_mag(uhatk_nz_0)).^2;

rhs = E(uhatk_nz_0);

rel_err = abs(lhs - rhs) ./ max(1, abs(rhs));
max_rel_err = max(rel_err);

%%
uhatk_nz_1 = k_mag > 0 & (kx.^2 + ky.^2) > 0;

kx_nz1 = kx(uhatk_nz_1);
ky_nz1 = ky(uhatk_nz_1);
kz_nz1 = kz(uhatk_nz_1);
k_mag_nz_2 = k_mag(uhatk_nz_1);

uhatk(uhatk_nz_1) = (alpha(uhatk_nz_1) .* k_mag_nz_2 .* ky_nz1  + beta(uhatk_nz_1) .* kx_nz1 .* kz_nz1) ./ (k_mag_nz_2 .* sqrt(kx_nz1.^2 + ky_nz1.^2));

vhatk(uhatk_nz_1) = (beta(uhatk_nz_1) .* ky_nz1 .* kz_nz1  - alpha(uhatk_nz_1) .* k_mag_nz_2 .* kx_nz1) ./ (k_mag_nz_2 .* sqrt(kx_nz1.^2 + ky_nz1.^2));

%%
uhatk_nz_2 = k_mag > 0;

kx_nz2 = kx(uhatk_nz_2);
ky_nz2 = ky(uhatk_nz_2);
kz_nz3 = kz(uhatk_nz_2);
k_mag_nz_3 = k_mag(uhatk_nz_2);

whatk(uhatk_nz_2) = -(beta(uhatk_nz_2) .* sqrt(kx_nz2.^2 + ky_nz2.^2) ./ k_mag_nz_3);

%%
%setting 0 and Nyquist modes to 0

zero_handling = (ll == 0) & (mm == 0) & (nn == 0);
nyq_handling = (ll == -Nx/2) | (mm == -Nx/2) | (nn == -Nx/2);

zero_or_nyq = zero_handling | nyq_handling;

uhatk(zero_or_nyq) = 0;
vhatk(zero_or_nyq) = 0;
whatk(zero_or_nyq) = 0;

%%

indep_k = (nn > 0) | ((nn == 0) & (mm > 0)) | ((nn == 0) & (mm == 0) & (ll > 0));
for i = 1:Nx
    for j = 1:Nx
        for k = 1:Nx
            if indep_k(i,j,k) == 0
                continue %means skip the rest
            end

            mean_0 = (ll(i,j,k) == 0) & (mm(i,j,k) == 0) & (nn(i,j,k) == 0);
            oddball_nyq_x = (ll(i,j,k) == -Nx/2);
            oddball_nyq_y = (mm(i,j,k) == -Nx/2);
            oddball_nyq_z = (nn(i,j,k) == -Nx/2);
    
            if mean_0 || oddball_nyq_x || oddball_nyq_y || oddball_nyq_z
                continue
            end
            neg_l = find(l == - ll(i,j,k));
            neg_m = find(m == - mm(i,j,k));
            neg_n = find(n == - nn(i,j,k));
        
            uhatk(neg_l, neg_m, neg_n) = conj(uhatk(i,j,k));
            vhatk(neg_l, neg_m, neg_n) = conj(vhatk(i,j,k));
            whatk(neg_l, neg_m, neg_n) = conj(whatk(i,j,k));

        end
    end
end

%%
COM = (kx .* uhatk) + (ky .* vhatk) + (kz .* whatk); %satisifies the continuity

max_COM = max(abs(COM(:)));

%%

u_ifft1 = ifftn(uhatk);
v_ifft2 = ifftn(vhatk);
w_ifft3 = ifftn(whatk);

%%
q_sq_rms_1 = mean(u_ifft1(:).^2 + v_ifft2(:).^2 + w_ifft3(:).^2);

q_sq_target = 3 * (u0^2);

scale = sqrt(q_sq_target / q_sq_rms_1);

uhatk  = scale * uhatk;
vhatk  = scale * vhatk;
whatk  = scale * whatk;

u_ifft_sc = ifftn(uhatk);
v_ifft_sc = ifftn(vhatk);
w_ifft_sc = ifftn(whatk);

q_sq_rms_scaled = mean( real(u_ifft_sc(:)).^2 + real(v_ifft_sc(:)).^2 + real(w_ifft_sc(:)).^2 );

end