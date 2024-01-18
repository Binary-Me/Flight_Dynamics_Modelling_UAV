%% State Variables
pn=0;pe=0;pd=0; %Inertial Position of MAV along Inertial Frame
dpn=0;dpe=0;dpd=0; %Rate of Change of Inertial Position of MAV along Inertial Frame
u=0;v=0;w=0; %Velocity of MAV measured along Body Frame
du=0;dv=0;dw=0; %Acceleration of MAV measured along Body Frame
p=0;q=0;r=0; %Angular rates of the MAV about the Body Frame 
dp=0;dq=0;dr=0; %Angular acceleration of the MAV about the Body Frame
phi=0;theta=0;psi=0; %Roll,Pitch and Yaw angles in degrees
phi=deg2rad(phi);theta=deg2rad(theta);psi=deg2rad(psi);
dphi=0;dtheta=0;dpsi=0; %Roll,Pitch and Yaw Rates

%% Simulation Settings
dt=0.01;
no=10000; %No of sample inputs

%Physical Parameters - Aerosonde UAV
mass=25;
Jx=0.824;
Jy=1.135;
Jz=1.759;
Jxz=0.120;
g=9.81;
S=0.55;
b=2.9;
c=0.19;
rho=1.268;
e=0.9;

% Motor Parameters
n_cells=12;
V_i=3.3;
Vmax=n_cells*V_i;
Dp=0.508;
Kv=0.0659;
Kq=0.0659;
Rmo=0.042;
i0=1.5;
CQ2=-0.01664;
CQ1=0.004970;
CQ0=0.005230;
CT2=-0.1079;
CT1=-0.06044;
CT0=0.09357;

% Aerodynamic Co-efficients
CL0=0.23;
CD0=0.043;
Cm0=0.0135;
Cla=5.61;
Cda=0.030;
Cma=-2.74;
Clq=7.95;
Cdq=0;
Cmq=-38.21;
Clde=0.13;
Cdde=0.0135;
Cmde=-0.99;
M=50;
Cdp=0;

Cy0=0;
Cl0=0;
Cn0=0;
Cyb=-0.83;
Clb=-0.13;
Cnb=0.073;
Cyp=0;
Clp=-0.51;
Cnp=-0.069;
Cyr=0;
Clr=0.25;
Cnr=-0.095;
Cyda=0.075;
Clda=0.17;
Cnda=-0.011;
Cydr=0.19;
Cldr=0.0024;
Cndr=-0.069;

% Mass Moment Co-efficients
G=Jx*Jz - Jxz*Jxz;
G1=(1/G)*Jxz*(Jx-Jy+Jz);
G2=(1/G)*(Jz*(Jz-Jy)+Jxz*Jxz);
G3=Jz/G;
G4=Jxz/G;
G5=(1/Jy)*(Jz-Jx);
G6=Jxz/Jy;
G7=(1/G)*((Jx-Jy)*Jx + Jxz*Jxz);
G8=Jx/G;

%% Inputs
deltae=-1:0.01:0;
delta_t=0;
delta_a=0;
delta_r=0;
count=0;
theta_j=zeros([101,4]);
%%
Vav=[25,50,100,200];
for kota=1:4
    u = Vav(1,kota);
    for dum=1:101
        delta_e=deltae(dum);
        angles=zeros([3,no]);
        velocity_b=zeros([3,no]);
        moment=zeros([3,no]);
        force=zeros([3,no]);
        v = 0;
        w = 0;
        p = 0;
        q = 0;
        r = 0;
        for i=1:no
            Va=sqrt(u*u+v*v+w*w);
            alpha=atan(w/u);
            beta=asin(v/Va);
            Cl=CL0+(Cla*alpha);
            Cd=CD0+(Cda*alpha);
            Cx=(-Cd*cos(alpha))+(Cl*sin(alpha));
            Cxq=(-Cdq*cos(alpha))+(Clq*sin(alpha));
            Cxde=(-Cdde*cos(alpha))+(Clde*sin(alpha));
            Cz=(-Cd*sin(alpha))-(Cl*cos(alpha));
            Czq=(-Cdq*sin(alpha))-(Clq*cos(alpha));
            Czde=(-Cdde*sin(alpha))-(Clde*cos(alpha));
        
            %Forces and Moments due to Propulsion
            Vin=Vmax*delta_t;
            ap = (-rho*Dp^5/4*pi^2)*CQ0;
            bp = (rho*Dp^4/2*pi)*CQ1*Va + Kq*Kv/Rmo;
            cp= rho*Dp^3*CQ2*Va^2 - Kq*Vin/Rmo + Kq*i0;
            omega_p= (-bp + sqrt(bp^2 - 4*ap*cp))/(2*ap);
            Tp=(rho*Dp^4*CT0/(4*pi^2))*omega_p^2 + (rho*Dp^3*CT1*Va/(2*pi))*omega_p + (rho*Dp^2*CT2*Va^2);
            Qp=(rho*Dp^5*CQ0/(4*pi^2))*omega_p^2 + (rho*Dp^4*CQ1*Va/(2*pi))*omega_p + (rho*Dp^3*CQ2*Va^2);
        
            %Forces
            a=Cx+(0.5*Cxq*c*q/Va)+(Cxde*delta_e);
            j=Cy0+(Cyb*beta)+(0.5*Cyp*b*p/Va)+(0.5*Cyr*b*r/Va)+(Cyda*delta_a)+(Cydr*delta_r);
            o=Cz+(0.5*Czq*c*q/Va)+(Czde*delta_e);
           
            fx = (-mass*g*sin(theta)) + (0.5*rho*Va*Va*S*a) + Tp;
            fy = (mass*g*cos(theta)*sin(phi))+(0.5*rho*Va*Va*S*j);
            fz = (mass*g*cos(theta)*cos(phi))+ (0.5*rho*Va*Va*S*o);
        
            %Moments
            d=b*(Cl0+(Clb*beta)+(0.5*Clp*b*p/Va)+(0.5*Clr*b*r/Va)+(Clda*delta_a)+(Cldr*delta_r));
            z=c*(Cm0+(Cma*alpha)+(0.5*Cmq*c*q/Va)+(Cmde*delta_e));
            f=b*(Cn0+(Cnb*beta)+(0.5*Cnp*b*p/Va)+(0.5*Cnr*b*r/Va)+(Cnda*delta_a)+(Cndr*delta_r));
            l=(0.5*rho*Va*Va*S*d) - Qp;
            m=(0.5*rho*Va*Va*S*z);
            n=(0.5*rho*Va*Va*S*f);
        
            angles(1:3,i)=[phi;theta;psi];
            velocity_b(1:3,i)=[u;v;w];
            
            moment(1:3,i)=[l;m;n];
            force(1:3,i)=[fx;fy;fz];
           
            wv=[r*v - q*w;p*w - r*u; q*u - p*v];
            temp = wv +(1/mass).*[fx;fy;fz];
            du_dt=temp(1);dv_dt=temp(2);dw_dt=temp(3);
            u=u+(du_dt*dt); %Euler Forward Integration
            v=v+(dv_dt*dt);
            w=w+(dw_dt*dt);
            R1=[cos(theta)*cos(psi) (sin(phi)*sin(theta)*cos(psi))-(cos(phi)*sin(psi)) (cos(phi)*sin(theta)*cos(psi))+(sin(phi)*sin(psi)); cos(theta)*sin(psi) (sin(phi)*sin(theta)*sin(psi))+(cos(phi)*cos(psi)) (cos(phi)*sin(theta)*sin(psi))-(sin(phi)*cos(psi)); -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta)];
            temp=R1*[u;v;w];
            dpn_dt=temp(1);dpe_dt=temp(2);dpd_dt=temp(3);
            pn=pn+(dpn_dt*dt); %Euler forward Integration
            pe=pe+(dpe_dt*dt);
            pd=pd+(dpd_dt*dt);      
        
            temp=[(G1*p*q)-(G2*q*r);(G5*p*r)-G6*(p*p-r*r);(G7*p*q)-(G1*q*r)] + [(G3*l +G4*n);(m/Jy);(G4*l + G8*n)];
            dp_dt=temp(1);dq_dt=temp(2);dr_dt=temp(3);
            p=p+(dp_dt*dt); %Euler Forward Integration
            q=q+(dq_dt*dt);
            r=r+(dr_dt*dt);
            R2=[1 sin(phi)*tan(theta) cos(phi)*tan(theta);0 cos(phi) -sin(phi); 0 (sin(phi)/cos(theta)) (cos(phi)/cos(theta))];
            temp=R2*[p;q;r];
            dphi_dt=temp(1);dtheta_dt=temp(2);dpsi_dt=temp(3);
            phi=phi+(dphi_dt*dt); %Euler Forward Integration
            theta=theta+(dtheta_dt*dt);
            psi=psi+(dpsi_dt*dt);
        
        end
        theta_j(dum,kota)=angles(2,no);
    end
end
%%
figure
plot(-1:0.01:0,mod(rad2deg(real(theta_j(:,1)))+180,360)-180);
hold on
plot(-1:0.01:0,mod(rad2deg(real(theta_j(:,2)))+180,360)-180);
plot(-1:0.01:0,mod(rad2deg(real(theta_j(:,3)))+180,360)-180);
plot(-1:0.01:0,mod(rad2deg(real(theta_j(:,4)))+180,360)-180);
%plot(-1:0.01:1,mod(rad2deg(real(alpha(1,:)))+180,360)-180);
xlabel('\delta_e')
ylabel("\theta_s")
legend("Va=25m/s","Va=50m/s","Va=100m/s","Va=200m/s");
title("\theta_s vs \delta_e")

% figure

% % gamma=angles(2,:)-alpha;
% angles=real(angles);
% subplot(3, 1, 1);
% plot(1:no,mod(rad2deg(angles(1,:))+180,360)-180);
% title("phi");
% subplot(3, 1, 2);
% plot(1:no,mod(rad2deg(angles(2,:))+180,360)-180);
% title("theta");
% subplot(3, 1, 3);
% plot(1:no,mod(rad2deg(angles(3,:))+180,360)-180);
% title("psi");
% % 
% figure
% subplot(3, 1, 1);
% Va_d=sqrt(velocity_b(1,:).^2+velocity_b(2,:).^2+velocity_b(3,:).^2);
% plot(1:no,Va_d);
% title("Va");
% subplot(3, 1, 2);
% for i=1:no
%     alpha_m(i)=atan(velocity_b(3,i)/velocity_b(1,i));
% end
% plot(1:no,mod(rad2deg(alpha_m)+180,360)-180);
% title("\alpha");
% subplot(3, 1, 3);
% gamma=angles(2,:)-alpha;
% plot(1:no,mod(rad2deg(gamma)+180,360)-180);
% title("\gamma");

% plot(1:no,velocity_b(1,:));
% hold on
% plot(1:no,velocity_b(2,:));
% hold on
% plot(1:no,velocity_b(3,:));
% 
% legend('u','v','w');
% 
% title("Velocity in Body Frame");
% figure
% plot (1:no,Va_d);
% hold on
% plot(1:no,gamma.*(mass*g));

% figure
% subplot(3, 1, 1);
% plot(1:no,moment(1,:));
% title("l");
% subplot(3, 1, 2);
% plot(1:no,moment(2,:));
% title("m");
% subplot(3, 1, 3);
% plot(1:no,moment(3,:));
% title("n");

% figure
% plot(1:no,force(1,:));
% hold on
% plot(1:no,force(2,:));
% hold on
% plot(1:no,force(3,:));
% legend('fx','fy','fz');
% title("Forces in Body Frame");