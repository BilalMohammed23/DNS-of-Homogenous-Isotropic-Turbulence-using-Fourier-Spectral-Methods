function [delta_t_max] = dt(max_vel,alpha,diff_gamma,h)

    delta_t_max_diff = (alpha * (h^2))/(3 * diff_gamma * pi^2);
    
    delta_t_max_conv = (alpha * h) / (sqrt(3) * pi * max_vel);

    delta_t_max = min(delta_t_max_diff, delta_t_max_conv);

end