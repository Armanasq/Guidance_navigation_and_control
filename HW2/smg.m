function [G0] = smg(L_b)
e = 0.0818191908425;
G0 = 9.7803253359*(1+0.001931853*sin(L_b)^2)/sqrt(1-e^2*sin(L_b)^2); 
end