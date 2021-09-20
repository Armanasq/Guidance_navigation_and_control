% GN&C
% HW2
% A.Asgharpoor
% email: A.Asgharpoor@ut.ac.ir
% University of Tehran
% FNST
% ref:
%   Paul D Groves - Principles of GNSS, inertial, and multi-sensor integrated navigation systems 2nd ed-Artech House (2013)

clc
clear all
close all

disp('GN&C  HW2')
disp('Dr. Amiri')
disp('A.Asgharpoor          email: A.Asgharpoor@ut.ac.ir')
disp('FNST')
disp('===================================================================================')
fprintf('\n\n')

% Mt. Everest as:
%   Latitude (Lb) 27deg 59min 16sec N
%   Longitude (?b) 86deg 56min 40sec E
%   height (hb) 8850 meters 
%   (derived by GPS in 1999):

Wie = 72.921151467e-6;              % rad/sec (WGS84)
lat_ev=(27+59/60+16/3600)*pi/180;
lon_ev=(86+56/60+40/3600)*pi/180;
H_ev=8850;                          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q1) 
% a)
%
%  convert from geodetic curvilinear lat, lon, and height to ECEF rectangular x, y, and z coordinates

r_e__e_b=llh2xyz(lat_ev, lon_ev, H_ev);

disp('Q1)');
fprintf('\n')

disp('a)');
fprintf('\n')
fprintf('r_e__e_b = [%8.0f \n              %8.0f \n              %8.0f ]', r_e__e_b(1), r_e__e_b(2), r_e__e_b(3))
fprintf('\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q1) b
%
% convert ECEF x, y, and z coordinates to lat, lon, and height
%Tabulate the error in lat (in deg), lon (in deg) and height (in m) for 1 to 10 iterations.

iterations = 10;
H_b_error = zeros(1,iterations);
Lambda_b_error=zeros(1,iterations);
L_b_error=zeros(1,iterations);			


disp('b) Error')
fprintf('iterations      D_lat(°)         D_lon(°)         D_h (m)')
fprintf('\n')

for i=1:iterations
        [L_b, lambda_b, h_b] = xyz2llh(r_e__e_b, i);
            L_b_error(i)=(lat_ev-L_b)*180/pi;	
            Lambda_b_error(i)=(lon_ev-lambda_b)*180/pi;	
            H_b_error(i)=(H_ev - h_b);				
        fprintf('     %2.0f      %+12.9f      %+12.9f     %+12.9f', i, L_b_error(i), Lambda_b_error(i), H_b_error(i))
        fprintf('\n')
end


fprintf('\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q1) c)
% The acceleration due to gravity at the ellipsoid
% Somigliana model page 70 equ 2.134

G0=smg(lat_ev);
disp('c) The acceleration due to gravity at the ellipsoid')
fprintf('G0 = %12.8f (m/s^2)', G0)
fprintf('\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Q1) d)
% The magnitude of the centrifugal acceleration at the ellipsoid and at the peak

disp('d) The magnitude of the centrifugal acceleration');
% Earth angular velocity
Oe_i_e = [0, -1, 0;	1,  0, 0; 0,  0, 0] * Wie;

disp('At The Ellipsoid');
r_e__e_b = llh2xyz(lat_ev, lon_ev, 0);
cent_epl = -Oe_i_e*Oe_i_e*r_e__e_b; 

disp('  ECEF coords')
fprintf('       Vector = [ %+12.9f , %+12.9f , %+12.9f ] (m/s^2)', cent_epl(1), cent_epl(2), cent_epl(3))
fprintf('\n')

fprintf('       Mag = %+10.9f (m/s^2)\n\n', norm(cent_epl))
fprintf('\n\n')


disp('At The Peak:');
r_e__e_b = llh2xyz(lat_ev, lon_ev, H_ev);
cent_epl = -Oe_i_e * Oe_i_e * r_e__e_b;  
disp(' ECEF coords')
fprintf('       Vector  = [ %+12.9f , %+12.9f , %+12.9f ] (m/s^2)', cent_epl(1), cent_epl(2), cent_epl(3))
fprintf('\n')

fprintf('       Mag = %+12.9f (m/s^2)', norm(cent_epl))
fprintf('\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part e)
% The magnitude of the gravitational attraction at the ellipsoid and at the peak
% page 72 equ 2.141

disp('e)');
disp('The magnitude of the gravitational attraction at the ellipsoid:');
r_e__i_b=llh2xyz(lat_ev,lon_ev,0);
[GMe__i_b]=GM__i_b(r_e__i_b);
disp(' ECEF coords')
fprintf('       Vector  = [ %+10.9f , %+10.9f , %+10.9f ] (m/s^2)', GMe__i_b(1), GMe__i_b(2), GMe__i_b(3))
fprintf('\n')
fprintf('       Mag = %+10.9f (m/s^2)', norm(GMe__i_b))
fprintf('\n\n')

disp('The magnitude of the gravitational attraction at the peak:');
r_e__i_b = llh2xyz(lat_ev, lon_ev, H_ev);
[GMe__i_b]=GM__i_b(r_e__i_b);

disp(' ECEF coords')
fprintf('       Vector  = [ %+10.9f , %+10.9f , %+10.9f ] (m/s^2)', GMe__i_b(1), GMe__i_b(2), GMe__i_b(3))
fprintf('\n')
fprintf('       Mag = %+10.9f (m/s^2)', norm(GMe__i_b))
fprintf('\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('===================================================================================')

% Q2
% The acceleration due to gravity as a function of lat & height
% a)
% page71 equ 2.139

disp('Q2)');
disp('a)')
disp('The acceleration due to gravity at the peak')
Gn__bD = gravity(lat_ev, H_ev);
fprintf('       = %10.8f (m/s^2)', Gn__bD)
fprintf('\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% b)
disp('b)')
fprintf('       Delta_g = %+10.9f - %+10.9f  = %+10.9f (m/s^2)', Gn__bD, G0, Gn__bD - G0)
fprintf('\n')
MyWeight = 185.188; % 85kg = 185.188 Ibs
DeltaWeight = MyWeight * (Gn__bD - G0) / G0;
fprintf('\n')
fprintf('       Delta Weight ~%4.2f lbs less ', abs(DeltaWeight))
fprintf('\n')