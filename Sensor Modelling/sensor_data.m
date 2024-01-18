%% Initial Set Up
close all; clear; clc;
[V, F] = spacecraftVFC;
R = [...
    0, 1, 0;...
    1, 0, 0;...
    0, 0, -1;...
];

V = V * R; % Transform vertices from NED to XYZ
handle = patch('Vertices', V, 'Faces', F,'Facecolor','yellow');
xlim([-100, 100]);
ylim([-100, 100]);
zlim([-100, 100]);
 title('Aircraft')
        xlabel('East')
        ylabel('North')
        zlabel('-Down')
view(3);

%% State Variables
pn=0;pe=0;pd=0; %Inertial Position of MAV along Inertial Frame
dpn=0;dpe=0;dpd=0; %Rate of Change of Inertial Position of MAV along Inertial Frame
%u=0;v=0;w=0; %Velocity of MAV measured along Body Frame
du=0;dv=0;dw=0; %Acceleration of MAV measured along Body Frame
p=0;q=0;r=0; %Angular rates of the MAV about the Body Frame 
dp=0;dq=0;dr=0; %Angular acceleration of the MAV about the Body Frame
phi=0;theta=0;psi=0; %Roll,Pitch and Yaw angles in degrees
phi=deg2rad(phi);theta=deg2rad(theta);psi=deg2rad(psi);
dphi=0;dtheta=0;dpsi=0; %Roll,Pitch and Yaw Rates

%% Simulation Settings
T=1; %Total time of simulation
dt=0.01; %Incremental Time Step
no=T/dt;

%% Sensor readings
%t,p,q,r,phi,theta,psi
true_values=zeros([no,7]);
sensor_values=zeros([no,7]);

%% Noise generator
noise=2*rand([1,6])-[1,1,1,1,1,1];
rms_value=[1,1,1,1,1,1];
noise=noise.*rms_value;

%% Inertial Properties - Aerosonde UAV
mass=13.5;
Jx=0.8244;
Jy=1.135;
Jz=1.759;
Jxz=0.1204;

%% Input Parameters
g=9.81;
S=0.55;
b=2.8956;
c=0.18994;
Sp=0.2027;
rho=1.2682;
k_mo=80;
Ktp=0;
Kw=0;
e=0.9;

CL0=0.28;
CD0=0.03;
Cm0=-0.02338;
Cla=3.45;
Cda=0.30;
Cma=-0.38;
Clq=0;
Cdq=0;
Cmq=-3.6;
Clde=-0.36;
Cdde=0;
Cmde=-0.25;
Cp=1;
Cdp=0.0437;
Cndr=-0.032;

Cy0=0;
Cl0=0;
Cn0=0;
Cyb=-0.98;
Clb=-0.12;
Cnb=0.25;
Cyp=0;
Clp=-0.26;
Cnp=0.022;
Cyr=0;
Clr=0.14;
Cnr=-0.35;
Cyda=0;
Clda=0.08;
Cnda=0.06;
Cydr=-0.17;
Cldr=0.105;

%% Mass Moment Co-efficients
G=Jx*Jz - Jxz*Jxz;
G1=(1/G)*Jxz*(Jx-Jy+Jz);
G2=(1/G)*(Jz*(Jz-Jy)+Jxz*Jxz);
G3=Jz/G;
G4=Jxz/G;
G5=(1/Jy)*(Jz-Jx);
G6=Jxz/Jy;
G7=(1/G)*((Jx-Jy)*Jx + Jxz*Jxz);
G8=Jx/G;

%% Computation
disp("Enter the ambient speed: ");
input_str = input('Va: ', 's');
Va= str2double(input_str);
disp("Enter the heading angle in degrees: ");
input_str = input('y: ', 's');
y= str2double(input_str);
y=y*pi/180;
disp("Enter the turn Radius: ");
input_str = input('R: ', 's');
turnR= str2double(input_str);
TRIM=states_calc(Va,y,turnR);
u=TRIM(5);
v=TRIM(6);
w=TRIM(7);
delta_a=TRIM(1);
delta_e=TRIM(2);
delta_r=TRIM(3);
delta_t=TRIM(4);

alpha=atan(w/u);
beta=asin(v/Va);
Cl=CL0+(Cla*alpha);
Cd=CD0+(Cda*alpha);
Cx=(-Cd*cos(alpha))+(Cl*sin(alpha));
Cxq=(-Cdq*cos(alpha))+(Clq*sin(alpha));
Cxde=(-Cdde*cos(alpha))+(Clde*sin(alpha));
Cz=(-Cd*sin(alpha))-(Cl*cos(alpha));
Czq=(-Cdq*sin(alpha))-(Clq*cos(alpha));
Czde=(-Cdde*sin(alpha))+(Clde*cos(alpha));

%Forces
a=Cx+(0.5*Cxq*c*q/Va)+(Cxde*delta_e);
j=Cy0+(Cyb*beta)+(0.5*Cyp*b*p/Va)+(0.5*Cyr*b*r/Va)+(Cyda*delta_a)+(Cydr*delta_r);
o=Cz+(0.5*Czq*c*q/Va)+(Czde*delta_e);
   
fx = (-mass*g*sin(theta)) + (0.5*rho*Va*Va*S*a) + (0.5*rho*Sp*Cp*((k_mo*delta_t*k_mo*delta_t)-(Va*Va)));
fy = (mass*g*cos(theta)*sin(phi))+(0.5*rho*Va*Va*S*j);
fz = (mass*g*cos(theta)*cos(phi))+ (0.5*rho*Va*Va*S*o);

%Moments
d=b*(Cl0+(Clb*beta)+(0.5*Clp*b*p/Va)+(0.5*Clr*b*r/Va)+(Clda*delta_a)+(Cldr*delta_r));
z=c*(Cm0+(Cma*alpha)+(0.5*Cmq*c*q/Va)+(Cmde*delta_e));
f=b*(Cn0+(Cnb*beta)+(0.5*Cnp*b*p/Va)+(0.5*Cnr*b*r/Va)+(Cnda*delta_a)+(Cndr*delta_r));
l=(0.5*rho*Va*Va*S*d) - (Ktp*Kw*delta_t*Kw*delta_t);
m=(0.5*rho*Va*Va*S*z);
n=(0.5*rho*Va*Va*S*f);

for i=1:no
    wv=[r*v - q*w;p*w - r*u; q*u - p*v];
    temp = wv +(1/mass).*[fx;fy;fz];
    du=temp(1);dv=temp(2);dw=temp(3);
    u=u+(du*dt); %Euler Forward Integration
    v=v+(dv*dt);
    w=w+(dw*dt);
    R1=[cos(theta)*cos(psi) (sin(phi)*sin(theta)*cos(psi))-(cos(phi)*sin(psi)) (cos(phi)*sin(theta)*cos(psi))+(sin(phi)*sin(psi)); cos(theta)*sin(psi) (sin(phi)*sin(theta)*sin(psi))+(cos(phi)*cos(psi)) (cos(phi)*sin(theta)*sin(psi))-(sin(phi)*cos(psi)); -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta)];
    temp=R1*[u;v;w];
    dpn=temp(1);dpe=temp(2);dpd=temp(3);
    pn=pn+(dpn*dt); %Euler forward Integration
    pe=pe+(dpe*dt);
    pd=pd+(dpd*dt);      
    
    temp=[(G1*p*q)-(G2*q*r);(G5*p*r)-G6*(p*p-r*r);(G7*p*q)-(G1*q*r)] + [(G3*l +G4*n);(m/Jy);(G4*l + G8*n)];
    dp=temp(1);dq=temp(2);dr=temp(3);
    p=p+(dp*dt); %Euler Forward Integration
    q=q+(dq*dt);
    r=r+(dr*dt);
    R2=[1 sin(phi)*tan(theta) cos(phi)*tan(theta);0 cos(phi) -sin(phi); 0 (sin(phi)/cos(theta)) (cos(phi)/cos(theta))];
    temp=R2*[p;q;r];
    dphi=temp(1);dtheta=temp(2);dpsi=temp(3);
    phi=phi+(dphi*dt); %Euler Forward Integration
    theta=theta+(dtheta*dt);
    psi=psi+(dpsi*dt);
    
    jitter=2*rand([1,6])-[1,1,1,1,1,1];
    true_values(i,1)=i*dt;
    sensor_values(i,1)=i*dt;
    true_values(i,2:7)=[p,q,r,phi,theta,psi];
    sensor_values(i,2:7)=[p,q,r,phi,theta,psi]+noise+jitter;

    V=V*R'; %XYZ to NED
    V = rotate(V', dphi, dtheta, dpsi)';
    V = translate(V', dpn, dpe, dpd)';
    V = V * R; % Transform vertices from NED to XYZ

    if exist('handle', 'var') && isvalid(handle)
        delete(handle);
    end
    handle = patch('Vertices', V, 'Faces', F,'Facecolor','yellow');
    view(3);

    % Adjust the limits based on the aircraft's position
    xlim([V(1,1) - 100, V(1,1) + 100]);
    ylim([V(1,2) - 100, V(1,2) + 100]);
    zlim([V(1,3) - 100, V(1,3) + 100]);

    pause(0.5);

    alpha=atan(w/u);
    beta=asin(v/Va);
    Cl=CL0+(Cla*alpha);
    Cd=CD0+(Cda*alpha);
    Cx=(-Cd*cos(alpha))+(Cl*sin(alpha));
    Cxq=(-Cdq*cos(alpha))+(Clq*sin(alpha));
    Cxde=(-Cdde*cos(alpha))+(Clde*sin(alpha));
    Cz=(-Cd*sin(alpha))-(Cl*cos(alpha));
    Czq=(-Cdq*sin(alpha))-(Clq*cos(alpha));
    Czde=(-Cdde*sin(alpha))+(Clde*cos(alpha));

    %Forces
    a=Cx+(0.5*Cxq*c*q/Va)+(Cxde*delta_e);
    j=Cy0+(Cyb*beta)+(0.5*Cyp*b*p/Va)+(0.5*Cyr*b*r/Va)+(Cyda*delta_a)+(Cydr*delta_r);
    o=Cz+(0.5*Czq*c*q/Va)+(Czde*delta_e);
   
    fx = (-mass*g*sin(theta)) + (0.5*rho*Va*Va*S*a) + (0.5*rho*Sp*Cp*((k_mo*delta_t*k_mo*delta_t)-(Va*Va)));
    fy = (mass*g*cos(theta)*sin(phi))+(0.5*rho*Va*Va*S*j);
    fz = (mass*g*cos(theta)*cos(phi))+ (0.5*rho*Va*Va*S*o);

    %Moments
    d=b*(Cl0+(Clb*beta)+(0.5*Clp*b*p/Va)+(0.5*Clr*b*r/Va)+(Clda*delta_a)+(Cldr*delta_r));
    z=c*(Cm0+(Cma*alpha)+(0.5*Cmq*c*q/Va)+(Cmde*delta_e));
    f=b*(Cn0+(Cnb*beta)+(0.5*Cnp*b*p/Va)+(0.5*Cnr*b*r/Va)+(Cnda*delta_a)+(Cndr*delta_r));
    l=(0.5*rho*Va*Va*S*d) - (Ktp*Kw*delta_t*Kw*delta_t);
    m=(0.5*rho*Va*Va*S*z);
    n=(0.5*rho*Va*Va*S*f);
end

disp('Exiting the program.');

%% Plots
figure
title("u,v,w");
plot(0.01:0.01:1,sensor_values(:,2),"Color",'b');
hold on
plot(0.01:0.01:1,true_values(:,2),"Color",'b','LineStyle','--');
hold on
plot(0.01:0.01:1,sensor_values(:,3),"Color",'g');
hold on
plot(0.01:0.01:1,true_values(:,3),"Color",'g','LineStyle','--');
hold on
plot(0.01:0.01:1,sensor_values(:,4),"Color",'r');
hold on
plot(0.01:0.01:1,true_values(:,4),"Color",'r','LineStyle','--');
%% Functions
function XYZ=rotate(XYZ, phi, theta, psi)
% Define rotation matrix
    R_roll = [...
        1, 0, 0;...
        0, cos(phi), -sin(phi);...
        0, sin(phi), cos(phi)];

    R_pitch = [...
        cos(theta), 0, sin(theta);...
        0, 1, 0;...
        -sin(theta), 0, cos(theta)];

    R_yaw = [...
        cos(psi), -sin(psi), 0;...
        sin(psi), cos(psi), 0;...
        0, 0, 1];

    R = R_roll * R_pitch * R_yaw;

    % Rotate vertices
    XYZ = R * XYZ;
end

function XYZ = translate(XYZ, pn, pe, pd)
    XYZ = XYZ + repmat([pn; pe; pd], 1, size(XYZ, 2));
end

function [V, F] = spacecraftVFC
    V = [...
        0 0 0;... % point 1
        -1.5 1.5 -1.5;... % point 2
        -1.5 -1.5 -1.5;... % point 3
        -1.5 -1.5 1.5;... % point 4
        -1.5 1.5 1.5;... % point 5
        -15 0 0;... % point 6
        -3 6.33 -1.33;... % point 7
        -7.6 6.33 -0.82;... % point 8
        -7.6 -6.33 -0.82;... % point 9
        -3 -6.33 -1.33;... % point 10
        -12 2.33 -0.33;... %point 11
        -15 2.33 0;... %point 12
        -15 -2.33 0;... %point 13
        -12 -2.33 -0.33;... point 14
        -12 0 -0.33;... %point 15
        -15 0 -4;... %point 16
        ];

    F = [...
        1,2,3,NaN;...
        1,4,5,NaN;...
        1,3,4,NaN;...
        1,2,5,NaN;...
        2,3,5,4;...
        2,3,6,NaN;...
        3,4,6,NaN;...
        2,5,6,NaN;...
        5,4,6,NaN;...
        7,8,9,10;
        11,12,13,14;
        6,15,16,NaN
        ];
end
