clear all;
clc;

%% Specify time
dt = 0.32; T = 90*60*300;

t = 0:dt:T;
%% Initializations
% The inertia vector of satellite
    inertia = [6400,-76.4,-25.6;-76.4,4730,-40;-25.6,-40,8160]'; % in kg*m^2
% The initial position and velocity of satellite in Earth-Centered Earth-Fixed    
    ini_pos = [1029.7743e3;6699.3469e3;3.7896e3]'; % in m
    ini_vel = [-6.2119e3;0.9524e3;4.3946e3]'; % in m/s
% Initial angular velocity and angular acceleration of satellite
    cang = [0.01;0.01;0.01]'; % rad/s
    angacc = [0;0;0]'; % rad/s^2
% Constants
    G = 6.67428e-11; % Earth gravitational constant in m^3/kg*s^2
    M = 5.972e24; % Earth mass in kg
    Re = 6371.2e3; % Radius of earth in m
% Final Answers
    ang_vel = zeros(length(t),3); % To store the angular velocity of satellite after dt time
    torque = zeros(length(t),3); % To store the torque of satellite after dt time
% Minimum Principal Momentum
    J = 4726.01952; % in kg*m^2

%% Calculated Values
% Scalar linear velocity of satellite
    linvel = sqrt(dot(ini_vel, ini_vel)); % in m/s
% The altitude of satellite from earth's surface
    lla = ecef2lla(ini_pos); % ecef2lla() converts ecef coordinates to latitude, longitude and altitude
    alti = lla(3); % in meters
% Distance of satellite from center of earth 
    Rc = Re + alti; % in m
% Time period of satellite
    timePeriod = 2*pi/sqrt(dot(cang, cang)); % in s^-1 Since Time Period = 2pi/(Angular Velocity)
% Accelaration of satellite due to earths gravity
    ini_acc = -ini_pos * (G * M / Rc^2) / sqrt(dot(-ini_pos, -ini_pos)); % in m/s^2 Since acceleration has value of GM/R^2 and in direction opposite to position vector
% Position of satellite after dt time
    ini_pos1 = ini_pos + ini_vel * dt + 0.5 * ini_acc * dt^2; % in m Assuming constant accelaration for dt time Since d = u*t + 0.5*a*t^2
% Angle Of Inclination Of Orbit
    normal = cross(ini_pos, ini_pos1); % Normal to orbit plane
    normalDotK = normal(3) / sqrt(dot(normal, normal)); % Dot product of unit normal vector with k^
    angleOfInclinationOfOrbit = acos(normalDotK); % in rad
% Positive Scalar Gain of Bdot Law
    k = 4*pi*(1 + sin(angleOfInclinationOfOrbit)*J)/timePeriod; % in kg*m^2*s Since k = 4*pi*(1+sin(angle of inclination))*Jmin/Torb

%% Initializations for loop
pos = ini_pos;
vel = ini_vel;
acc = ini_acc;

%% Main loop
 for i=1:length(t)
     
    veli = vel;
    posi = pos;
    acci = acc;
    cangi = cang;
    
    lla = ecef2lla(posi);
    lat = lla(1); % Latitude
    long = lla(2); % Longitude
    alti = lla(3); % Altitude
    alti = min(alti, 6e5); % Since igrfmagm has a limit on altitude of 6e5
    
    [mag_field_vector1,hor_intensity,declination,inclination,total_intensity] = igrfmagm(alti,lat,long,decyear(2015,7,4),12); % igrfmagm used to calculate the magnetic feild of earth at particular position
    mag_field_vector2 = mag_field_vector1 * 1e-9; % in T Since the function returns the value in nT
    mag_field_vector2 = mag_field_vector2.'; % Taking transpose of the magnetic feild
    
    % Changing frame of reference of Magnetic Feild from GCI to Satellite Body
    if i == 1
%       Initialize with values = 0
        theta_x = 0; % in rad
        theta_y = 0; % in rad
        theta_z = 0; % in rad
    else 
        theta_x = theta_x + (cangi(1) * dt)+(0.5 * angacc(1) * dt^2); % in rad Since Change in angle = Angular Velocity * dt + 0.5 * Angular Accelaration * dt^2
        theta_y = theta_y + (cangi(2) * dt)+(0.5 * angacc(2) * dt^2); % in rad Since Change in angle = Angular Velocity * dt + 0.5 * Angular Accelaration * dt^2
        theta_z = theta_z + (cangi(3) * dt)+(0.5 * angacc(3) * dt^2); % in rad Since Change in angle = Angular Velocity * dt + 0.5 * Angular Accelaration * dt^2

    end

    rx = [1,0,0;0,cos(theta_x),-sin(theta_x);0,sin(theta_x),cos(theta_x)]; % Rotation Matrix Rx
    ry = [cos(theta_y),0,sin(theta_y);0,1,0;-sin(theta_y),0,cos(theta_y)]; % Rotation Matrix Ry
    rz = [cos(theta_z),-sin(theta_z),0;sin(theta_z),cos(theta_z),0;0,0,1]; % Rotation Matrix Rz
    mag_field_vector = rx * ry * rz * mag_field_vector2; % Since Rotation vector R = Rx*Ry*Rz
    
%% B-dot
%   Determinant of Magnetic Feild
        detb = sqrt(dot((mag_field_vector),(mag_field_vector)));
%   Magnetic Dipole Moment
        m = ((k)/detb*norm(mag_field_vector))*cross(cang, mag_field_vector); % in A*m^2 Since m = -k (w x b)/||B||
%   Torque
        vtorque = cross(m, mag_field_vector); % in N*m Since T = m x B
%   Angular Accelaration of Satellite
        angacc(1) = (vtorque(1) - ((inertia(3, :) - inertia(2, :)) * cang(2) * cang(3))) / inertia(1,:); % in rad/s^2 Since I1 x angacc1 + (I3 - I2)w2w3  = T1
        angacc(2) = (vtorque(2) - ((inertia(1, :) - inertia(3, :)) * cang(3) * cang(1))) / inertia(2,:); % in rad/s^2 Since I2 x angacc2 + (I1 - I3)w3w1  = T2
        angacc(3) = (vtorque(2) - ((inertia(2, :) - inertia(1, :)) * cang(1) * cang(2))) / inertia(3,:); % in rad/s^2 Since I3 x angacc3 + (I2 - I1)w1w2  = T3
%% Updating loop values

    cang = cangi + (angacc * dt); % calculates the new angular velocity
    pos = posi + (veli * dt) + (0.5 * acci * dt^2); % calculates the new position assuming constant acc for dt time
    vel = veli + (acci * dt); % calculates the new velocity
    acc = -pos * (G * M / Rc^2) / sqrt(dot(pos, pos)); % calculates the new acceleration
% Display the new angular velocity
    disp('new angular velocity');
    disp(cang);
    ang_vel(i,:) = cang;
    torque(i,:) = vtorque;
     
 end

%% Plots
% Plot for Angular velocity
    figure(1)
    subplot(3,1,1)
    plot(t/60,ang_vel(:,1));
    subplot(3,1,2)
    plot(t/60,ang_vel(:,2));
    subplot(3,1,3)
    plot(t/60,ang_vel(:,3));
% Plot for Torque
    figure(2)
    subplot(3,1,1)
    plot(t,torque(:,1));
    subplot(3,1,2)
    plot(t,torque(:,2));
    subplot(3,1,3)
    plot(t,torque(:,3));
