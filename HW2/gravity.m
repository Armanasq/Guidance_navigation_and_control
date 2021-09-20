function [Gn_b_d] = gravity(L_b, h_b)
Wie = 72.921151467e-6;  
G0 = smg(L_b);	
mi_u = 3.986004418e14;  
flatting = 1/298.257223563; 
       
rp = 6356752.3142;      
r0 = 6378137.0;        

Gn_b_d=G0*(1-(2/r0)*(1+flatting+(Wie^2*r0^2*rp)/mi_u)*h_b + 3*h_b^2/(r0^2));
end