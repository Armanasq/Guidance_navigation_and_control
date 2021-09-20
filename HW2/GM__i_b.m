function [GM_x__i_b] = GM__i_b(Rx_i_b)
r0 = 6378137.0;
mi_u = 3.986004418e14;               
j2 = 1.082627e-3;                                

r_x__i_b_norm = norm(Rx_i_b);                
t_temp = (Rx_i_b(3) / r_x__i_b_norm)^2;       
r_temp = [(1 - 5*(t_temp)) * Rx_i_b(1); (1 - 5*(t_temp)) * Rx_i_b(2); (3 - 5*(t_temp)) * Rx_i_b(3)];

GM_x__i_b=(-mi_u/r_x__i_b_norm^3)*(Rx_i_b+ (3/2)*j2*(r0^2/r_x__i_b_norm^2)*r_temp); 
end