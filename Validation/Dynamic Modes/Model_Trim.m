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
k_mo=80;
Ktp=0;
Kw=0;
Sp=0.2027;
Cp=1;

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
Va=25;
y=(pi/180)*0;
R=100;
trim=TRIM(Va,y,R);
alpha=trim(1);
beta=trim(2);
mu=trim(3);
delta_e=trim(4);
delta_t=trim(5);
delta_a=trim(6);
delta_r=trim(7);

u=Va*cos(alpha)*cos(beta);
v=Va*sin(beta);
w=Va*sin(alpha)*cos(beta);

stheta=sin(y)*cos(beta)*cos(alpha) + cos(y)*sin(mu)*sin(beta)*cos(alpha) + cos(y)*cos(mu)*sin(alpha);
theta=asin(stheta);
cphi=(-sin(y)*cos(beta)*sin(alpha) - cos(y)*sin(mu)*sin(beta)*sin(alpha) + cos(y)*cos(mu)*cos(alpha))/cos(theta);
sphi=(-sin(y)*sin(beta) + cos(y)*sin(mu)*cos(beta))/cos(theta);
phi=atan(sphi/cphi);


p=(Va*cos(y)/R)*(-sin(theta));
q=(Va*cos(y)/R)*sin(phi)*cos(theta);
r=(Va*cos(y)/R)*cos(phi)*cos(theta);
count=0;
angles=zeros([3,no]);
velocity_b=zeros([3,no]);
moment=zeros([3,no]);
force=zeros([3,no]);
loc=zeros([3,no]);
%%
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
    Tp=(0.5*rho*Sp*Cp*((k_mo*delta_t*k_mo*delta_t)-(Va*Va)));
    Qp=(Ktp*Kw*delta_t*Kw*delta_t);
   
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
    loc(1:3,i)=[pn;pe;pd];
    
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


figure
% alpha=atan(velocity_b(3,:)/velocity_b(1,:));
% gamma=angles(2,:)-alpha;
angles=real(angles);
subplot(3, 1, 1);
plot(1:no,mod(rad2deg(angles(1,:))+180,360)-180);
title("phi");
subplot(3, 1, 2);
plot(1:no,mod(rad2deg(angles(2,:))+180,360)-180);
title("theta");
subplot(3, 1, 3);
plot(1:no,mod(rad2deg(angles(3,:))+180,360)-180);
title("psi");
% 
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
% for i=1:no
%     beta_m(i)=atan(velocity_b(2,i)/Va_d(i));
% end
% plot(1:no,mod(rad2deg(beta_m)+180,360)-180);
% title("\beta");

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

figure
subplot(3, 1, 1);
plot(1:no,moment(1,:));
title("l");
subplot(3, 1, 2);
plot(1:no,moment(2,:));
title("m");
subplot(3, 1, 3);
plot(1:no,moment(3,:));
title("n");

figure
plot(1:no,loc(3,:));
title("pd");

% figure
% plot(1:no,force(1,:));
% hold on
% plot(1:no,force(2,:));
% hold on
% plot(1:no,force(3,:));
% legend('fx','fy','fz');
% title("Forces in Body Frame");

%%
function trim=TRIM(Va,y,R)
    alpha=0;
    beta=0;
    mu=(cos(y)*Va)^2/(9.8*R);
    del_t=0;
    del_e=0;
    del_a=0;
    del_r=0;
    
    u=Va*cos(alpha)*cos(beta);
    v=Va*sin(beta);
    w=Va*sin(alpha)*cos(beta);
    
    
    stheta=(sin(y)*cos(alpha)*cos(beta))+(cos(y)*cos(mu)*sin(alpha))+ (cos(y)*sin(mu)*sin(beta)*cos(alpha));
    theta=asin(stheta);
    cphi=((-sin(y)*sin(alpha)*cos(beta))+(cos(y)*cos(mu)*cos(alpha))-(cos(y)*sin(mu)*sin(beta)*sin(alpha)))/cos(theta);
    sphi=((-sin(y)*sin(beta))+(cos(y)*sin(mu)*cos(beta)))/cos(theta);
    phi=atan(sphi/cphi);
    
    
    p=(Va*cos(y)/R)*(-sin(theta));
    q=(Va*cos(y)/R)*sin(phi)*cos(theta);
    r=(Va*cos(y)/R)*cos(phi)*cos(theta);
    
    
    
    
    %%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%
    
    rho=1.2682;
    S=0.55;
    m=13.5;
    
    Jx=0.8244;
    Jy=1.135;
    Jz=1.759;
    Jxz=0.1204;
    
    b=2.8956;
    c=0.18994;
    Sprop=0.2027;
    Kmotor=80;
    KTp=0;
    Komg=0;
    e=0.9;
    
    CLo=0.28;
    CDo=0.03;
    Cmo=-0.02338;
    CLalpha=3.45;
    CDalpha=0.30;
    Cmalpha=-0.38;
    CLq=0;
    CDq=0;
    Cmq=-3.6;
    CLdele=-0.36;
    CDdele=0;
    Cmdele=-0.5;
    Cprop=1;
    M=50;
    alphao=0.4712;
    epsilon=0.1592;
    CDp=0.0437;
    Cndelr=-0.032;
    
    CYo=0;
    Clo=0;
    Cno=0;
    CYbeta=-0.98;
    Clbeta=-0.12;
    Cnbeta=0.25;
    CYp=0;
    Clp=-0.26;
    Cnp=0.022;
    CYr=0;
    Clr=0.14;
    Cnr=-0.35;
    CYdela=0;
    Cldela=0.08;
    Cndela=0.06;
    CYdelr=-0.17;
    Cldelr=0.105;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% constants dependent of alpha
    
    CL=CLo+(CLalpha*alpha);
    CD=CDo+(CDalpha*alpha);
        
    CX=(-CD*cos(alpha))+(CL*sin(alpha));
    CXq=(-CDq*cos(alpha))+(CLq*sin(alpha));
    CXdele=(-CDdele*cos(alpha))+(CLdele*sin(alpha));
    CZ=(-CD*sin(alpha))-(CL*cos(alpha));
    CZq=(-CDq*sin(alpha))-(CLq*cos(alpha));
    CZdele=(-CDdele*sin(alpha))-(CLdele*cos(alpha));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% geometric constants
    
    % gamma 1-8
    
    g=(Jx*Jz)-(Jxz^2);
    g1=Jxz*(Jx-Jy+Jz)/g;
    g2=((Jz*(Jz-Jy))+Jxz^2)/g;
    g3=Jz/g;
    g4=Jxz/g;
    g5=(Jz-Jx)/Jy;
    g6=Jxz/Jy;
    g7=((Jx*(Jx-Jy))+Jxz^2)/g;
    g8=Jx/g;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Cpo=(g3*Clo)+(g4*Cno);
    Cpbeta=(g3*Clbeta)+(g4*Cnbeta);
    Cpp=(g3*Clp)+(g4*Cnp);
    Cpr=(g3*Clr)+(g4*Cnr);
    Cpdela=(g3*Cldela)+(g4*Cndela);
    Cpdelr=(g3*Cldelr)+(g4*Cndelr);
    Cro=(g4*Clo)+(g8*Cno);
    Crbeta=(g4*Clbeta)+(g8*Cnbeta);
    Crp=(g4*Clp)+(g8*Cnp);
    Crr=(g4*Clr)+(g8*Cnr);
    Crdela=(g4*Cldela)+(g8*Cndela);
    Crdelr=(g4*Cldelr)+(g8*Cndelr);
    
    % X=[u,v,w,phi,theta,psi,p,q,r];
    
    
    X_guess=[alpha;beta;mu;del_e;del_t;del_a;del_r];
    options=optimoptions('fsolve','TolX',1e-6);
    trim=fsolve(@(X)sct(X,Va,y,R),X_guess,options);
end

function residue=sct(X,Va,y,R)  % sct~ steady coordinated turn



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    alpha=X(1);
    beta=X(2);
    mu=X(3);
    del_e=X(4);
    del_t=X(5);
    del_a=X(6);
    del_r=X(7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rho=1.2682;
    S=0.55;
    m=13.5;
    
    Jx=0.8244;
    Jy=1.135;
    Jz=1.759;
    Jxz=0.1204;

    b=2.8956;
    c=0.18994;
    Sprop=0.2027;
    Kmotor=80;
    KTp=0;
    Komg=0;
    e=0.9;

    CLo=0.28;
    CDo=0.03;
    Cmo=-0.02338;
    CLalpha=3.45;
    CDalpha=0.30;
    Cmalpha=-0.38;
    CLq=0;
    CDq=0;
    Cmq=-3.6;
    CLdele=-0.36;
    CDdele=0;
    Cmdele=-0.5;
    Cprop=1;
    M=50;
    alphao=0.4712;
    epsilon=0.1592;
    CDp=0.0437;
    Cndelr=-0.032;

    CYo=0;
    Clo=0;
    Cno=0;
    CYbeta=-0.98;
    Clbeta=-0.12;
    Cnbeta=0.25;
    CYp=0;
    Clp=-0.26;
    Cnp=0.022;
    CYr=0;
    Clr=0.14;
    Cnr=-0.35;
    CYdela=0;
    Cldela=0.08;
    Cndela=0.06;
    CYdelr=-0.17;
    Cldelr=0.105;

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    g=(Jx*Jz)-(Jxz^2);
    g1=Jxz*(Jx-Jy+Jz)/g;
    g2=((Jz*(Jz-Jy))+Jxz^2)/g;
    g3=Jz/g;
    g4=Jxz/g;
    g5=(Jz-Jx)/Jy;
    g6=Jxz/Jy;
    g7=((Jx*(Jx-Jy))+Jxz^2)/g;
    g8=Jx/g;
    g9=1/Jy;
    
    
    
    CL=CLo+(CLalpha*alpha);
    CD=CDo+(CDalpha*alpha);

    CX=(-CD*cos(alpha))+(CL*sin(alpha));
    CXq=(-CDq*cos(alpha))+(CLq*sin(alpha));
    CXdele=(-CDdele*cos(alpha))+(CLdele*sin(alpha));
    CZ=(-CD*sin(alpha))-(CL*cos(alpha));
    CZq=(-CDq*sin(alpha))-(CLq*cos(alpha));
    CZdele=(-CDdele*sin(alpha))-(CLdele*cos(alpha));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% geometric constants

    % gamma 1-8

  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Cpo=(g3*Clo)+(g4*Cno);
    Cpbeta=(g3*Clbeta)+(g4*Cnbeta);
    Cpp=(g3*Clp)+(g4*Cnp);
    Cpr=(g3*Clr)+(g4*Cnr);
    Cpdela=(g3*Cldela)+(g4*Cndela);
    Cpdelr=(g3*Cldelr)+(g4*Cndelr);
    Cro=(g4*Clo)+(g8*Cno);
    Crbeta=(g4*Clbeta)+(g8*Cnbeta);
    Crp=(g4*Clp)+(g8*Cnp);
    Crr=(g4*Clr)+(g8*Cnr);
    Crdela=(g4*Cldela)+(g8*Cndela);
    Crdelr=(g4*Cldelr)+(g8*Cndelr);
    
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    stheta=(sin(y)*cos(alpha)*cos(beta))+(cos(y)*cos(mu)*sin(alpha))+ (cos(y)*sin(mu)*sin(beta)*cos(alpha));
    theta=asin(stheta);
    cphi=((-sin(y)*sin(alpha)*cos(beta))+(cos(y)*cos(mu)*cos(alpha))-(cos(y)*sin(mu)*sin(beta)*sin(alpha)))/cos(theta);
    sphi=((-sin(y)*sin(beta))+(cos(y)*sin(mu)*cos(beta)))/cos(theta);
    phi=atan(sphi/cphi);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E11=((Va*cos(y)/R)^2)*[-g1*cos(theta)*sin(theta)*sin(phi)-g2*cos(phi)*sin(phi)*(cos(theta))^2;
        -g5*cos(theta)*sin(theta)*cos(phi)-g6*((sin(theta))^2-(cos(theta)*cos(phi))^2);
        -g7*cos(theta)*sin(theta)*sin(phi)-g1*cos(phi)*sin(phi)*(cos(theta))^2];
    
    E12=0.5*rho*S*Va^2*[g3,0,g4;0,g9,0;g4,0,g8];
    
    E13=[-KTp*(Komg*del_t)^2;0;0];
    
    E14=[b,0,0;0,c,0;0,0,b];
    
    E15=[Cpo+Cpbeta*beta+((Cpp*sin(theta)+Cpr*cos(theta)*cos(phi))*b*cos(y)*0.5/R)+Cpdela*del_a+Cpdelr*del_r;
        Cmo+Cmalpha*alpha+(Cmq*c*cos(theta)*sin(phi)*cos(y)*0.5/R)+Cmdele*del_e;
        Cro+Crbeta*beta+((-Crp*sin(theta)+Crr*cos(theta)*cos(phi))*b*cos(y)*0.5/R)+Crdela*del_a+Crdelr*del_r];
    
    E1=E11+E12*(E13+(E14*E15));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    E21=-(((Va*cos(y))^2)/R)*[-cos(mu)*sin(beta)*cos(alpha)+sin(mu)*sin(alpha);
        cos(mu)*cos(beta);
        -cos(mu)*sin(beta)*sin(alpha)-sin(mu)*cos(alpha)];
    
    E22=9.8*[-sin(theta);cos(theta)*sin(phi);cos(theta)*cos(phi)];
    
    E23=(rho*Sprop*Cprop*0.5/m)*[(Kmotor*del_t)^2-(Va*cos(alpha)*cos(beta))^2;0;0];
    
    E24=((rho*S*Va^2)/(2*m))*eye(3);%*[-cos(alpha),0,sin(alpha);0,1,0;-sin(alpha),0,-cos(alpha)];
    
    %CL=CLo+(CLalpha*alpha);
    
    p=-sin(theta)*Va*cos(y)/R;
    q=sin(phi)*cos(theta)*Va*cos(y)/R;
    r=cos(phi)*cos(theta)*Va*cos(y)/R;
    
    E25=[CX+(CXq*q*c*0.5/Va)+(CXdele*del_e);
        CYo+(CYbeta*beta)+((CYp*p+CYr*r)*(b*0.5/Va))+(CYdela*del_a)+(CYdelr*del_r);
        CZ+(CZ*alpha)+(CZq*q*c*0.5/Va)+(CZdele*del_e)];

    E2=E21+E22+E23+(E24*E25);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    E3=(CYbeta*beta)+((CYp*p+CYr*r)*(b*0.5/Va))+(CYdela*del_a)+(CYdelr*del_r);
    
    residue=[E1;E2;E3];
    
end