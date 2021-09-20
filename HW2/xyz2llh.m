function [L_b, lambda_b, h_b] = xyz2llh(r_e__e_b, it)
e = 0.0818191908425;           
R0 = 6378137;   
h_b = 0;    
RE = R0;

x = r_e__e_b(1);        
y = r_e__e_b(2);
z = r_e__e_b(3);
rr = sqrt(x^2 + y^2);
lambda_b = atan2(y, x);       

    

    for i=1:it   
        sin_L_b = z / ((1-e^2)*RE + h_b);
        L_b     = atan((z+e*e*RE*sin_L_b) / rr);
        RE      = R0 / sqrt(1-e^2*sin_L_b^2);	
        h_b     = rr / cos(L_b) - RE;        
    end
end
