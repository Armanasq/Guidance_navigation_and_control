function [r_e__e_b] = llh2xyz(L_b, lambda_b, h_b)

e = 0.0818191908425;          
r0 = 6378137;        

Re = r0 / sqrt(1 - e^2*(sin(L_b))^2);


x=(Re + h_b) * cos(L_b) * cos(lambda_b);
y=(Re + h_b) * cos(L_b) * sin(lambda_b);
z=((1 - e^2)*Re + h_b)  * sin(L_b);
r_e__e_b = [x; y; z];
end
