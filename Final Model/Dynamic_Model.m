function Outputs= Dynamic_Model(Va,y,turnR)
    
    Struct=Model_Parameters();
    TRIM=States_calc(Va,y,turnR,Struct);

    %Variables
    pn=0;pe=0;pd=0; %Inertial Position of MAV along Inertial Frame
    dpn=0;dpe=0;dpd=0; %Rate of Change of Inertial Position of MAV along Inertial Frame
    %u=0;v=0;w=0; %Velocity of MAV measured along Body Frame
    du=0;dv=0;dw=0; %Acceleration of MAV measured along Body Frame
    p=0;q=0;r=0; %Angular rates of the MAV about the Body Frame 
    dp=0;dq=0;dr=0; %Angular acceleration of the MAV about the Body Frame
    phi=0;theta=0;psi=0; %Roll,Pitch and Yaw angles in degrees
    phi=deg2rad(phi);theta=deg2rad(theta);psi=deg2rad(psi);
    dphi=0;dtheta=0;dpsi=0; %Roll,Pitch and Yaw Rates

    u=TRIM(5);
    v=TRIM(6);
    w=TRIM(7);
    delta_a=TRIM(1);
    delta_e=TRIM(2);
    delta_r=TRIM(3);
    delta_t=TRIM(4);
    alpha=atan(w/u);
    beta=asin(v/Va);

    %Calculating Dependant Parameters
    Cl=Struct.CL0+(Struct.Cla*alpha);
    Cd=Struct.CD0+(Struct.Cda*alpha);
    Cx=(-Cd*cos(alpha))+(Cl*sin(alpha));
    Cxq=(-Struct.Cdq*cos(alpha))+(Struct.Clq*sin(alpha));
    Cxde=(-Struct.Cdde*cos(alpha))+(Struct.Clde*sin(alpha));
    Cz=(-Cd*sin(alpha))-(Cl*cos(alpha));
    Czq=(-Struct.Cdq*sin(alpha))-(Struct.Clq*cos(alpha));
    Czde=(-Struct.Cdde*sin(alpha))+(Struct.Clde*cos(alpha));

    %Forces
    a=Cx+(0.5*Cxq*Struct.c*q/Va)+(Cxde*delta_e);
    j=Struct.Cy0+(Struct.Cyb*beta)+(0.5*Struct.Cyp*Struct.b*p/Va)+(0.5*Struct.Cyr*Struct.b*r/Va)+(Struct.Cyda*delta_a)+(Struct.Cydr*delta_r);
    o=Cz+(0.5*Czq*Struct.c*q/Va)+(Czde*delta_e);
   
    fx = (-Struct.mass*Struct.g*sin(theta)) + (0.5*Struct.rho*Va*Va*Struct.S*a) + (0.5*Struct.rho*Struct.Sp*Struct.Cp*((Struct.k_mo*delta_t*Struct.k_mo*delta_t)-(Va*Va)));
    fy = (Struct.mass*Struct.g*cos(theta)*sin(phi))+(0.5*Struct.rho*Va*Va*Struct.S*j);
    fz = (Struct.mass*Struct.g*cos(theta)*cos(phi))+ (0.5*Struct.rho*Va*Va*Struct.S*o);

    %Moments
    d=Struct.b*(Struct.Cl0+(Struct.Clb*beta)+(0.5*Struct.Clp*Struct.b*p/Va)+(0.5*Struct.Clr*Struct.b*r/Va)+(Struct.Clda*delta_a)+(Struct.Cldr*delta_r));
    z=Struct.c*(Struct.Cm0+(Struct.Cma*alpha)+(0.5*Struct.Cmq*Struct.c*q/Va)+(Struct.Cmde*delta_e));
    f=Struct.b*(Struct.Cn0+(Struct.Cnb*beta)+(0.5*Struct.Cnp*Struct.b*p/Va)+(0.5*Struct.Cnr*Struct.b*r/Va)+(Struct.Cnda*delta_a)+(Struct.Cndr*delta_r));
    l=(0.5*Struct.rho*Va*Va*Struct.S*d) - (Struct.Ktp*Struct.Kw*delta_t*Struct.Kw*delta_t);
    m=(0.5*Struct.rho*Va*Va*Struct.S*z);
    n=(0.5*Struct.rho*Va*Va*Struct.S*f);

    T=Time_inputs();
    no=T.Time/T.dt;

    rates_ned=zeros([3,no+1]);
    rates_ned(1:3,1)=[dpn;dpe;dpd];

    rates_euler=zeros([3,no+1]);
    rates_euler(1:3,1)=[dphi;dtheta;dpsi];

    omega_actual=zeros([3,no+1]);
    omega_actual(1:3,1)=[p;q;r];

    euler_angles_actual=zeros([3,no+1]);
    euler_angles_actual(1:3,1)=[phi;theta;psi];

    for i=1:no
        wv=[r*v - q*w;p*w - r*u; q*u - p*v];
        temp = wv +(1/Struct.mass).*[fx;fy;fz];
        du=temp(1);dv=temp(2);dw=temp(3);
        u=u+(du*T.dt); %Euler Forward Integration
        v=v+(dv*T.dt);
        w=w+(dw*T.dt);
        R1=[cos(theta)*cos(psi) (sin(phi)*sin(theta)*cos(psi))-(cos(phi)*sin(psi)) (cos(phi)*sin(theta)*cos(psi))+(sin(phi)*sin(psi)); cos(theta)*sin(psi) (sin(phi)*sin(theta)*sin(psi))+(cos(phi)*cos(psi)) (cos(phi)*sin(theta)*sin(psi))-(sin(phi)*cos(psi)); -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta)];
        temp=R1*[u;v;w];
        dpn=temp(1);dpe=temp(2);dpd=temp(3);
        pn=pn+(dpn*T.dt); %Euler forward Integration
        pe=pe+(dpe*T.dt);
        pd=pd+(dpd*T.dt);      
    
        temp=[(Struct.G1*p*q)-(Struct.G2*q*r);(Struct.G5*p*r)-Struct.G6*(p*p-r*r);(Struct.G7*p*q)-(Struct.G1*q*r)] + [(Struct.G3*l +Struct.G4*n);(m/Struct.Jy);(Struct.G4*l + Struct.G8*n)];
        dp=temp(1);dq=temp(2);dr=temp(3);
        p=p+(dp*T.dt); %Euler Forward Integration
        q=q+(dq*T.dt);
        r=r+(dr*T.dt);
        R2=[1 sin(phi)*tan(theta) cos(phi)*tan(theta);0 cos(phi) -sin(phi); 0 (sin(phi)/cos(theta)) (cos(phi)/cos(theta))];
        temp=R2*[p;q;r];
        dphi=temp(1);dtheta=temp(2);dpsi=temp(3);
        phi=phi+(dphi*T.dt); %Euler Forward Integration
        theta=theta+(dtheta*T.dt);
        psi=psi+(dpsi*T.dt);

        omega_actual(1:3,i+1)=[p;q;r];
        euler_angles_actual(1:3,i+1)=[phi;theta;psi];
        rates_euler(1:3,i+1)=[dphi*T.dt;dtheta*T.dt;dpsi*T.dt];
        rates_ned(1:3,i+1)=[dpd*T.dt;dpe*T.dt;dpd*T.dt];

        %Updating net forces and moments acting after small displacement
        alpha=atan(w/u);
        beta=asin(v/Va);

        %Calculating Dependant Parameters
        Cl=Struct.CL0+(Struct.Cla*alpha);
        Cd=Struct.CD0+(Struct.Cda*alpha);
        Cx=(-Cd*cos(alpha))+(Cl*sin(alpha));
        Cxq=(-Struct.Cdq*cos(alpha))+(Struct.Clq*sin(alpha));
        Cxde=(-Struct.Cdde*cos(alpha))+(Struct.Clde*sin(alpha));
        Cz=(-Cd*sin(alpha))-(Cl*cos(alpha));
        Czq=(-Struct.Cdq*sin(alpha))-(Struct.Clq*cos(alpha));
        Czde=(-Struct.Cdde*sin(alpha))+(Struct.Clde*cos(alpha));

        %Forces
        a=Cx+(0.5*Cxq*Struct.c*q/Va)+(Cxde*delta_e);
        j=Struct.Cy0+(Struct.Cyb*beta)+(0.5*Struct.Cyp*Struct.b*p/Va)+(0.5*Struct.Cyr*Struct.b*r/Va)+(Struct.Cyda*delta_a)+(Struct.Cydr*delta_r);
        o=Cz+(0.5*Czq*Struct.c*q/Va)+(Czde*delta_e);
   
        fx = (-Struct.mass*Struct.g*sin(theta)) + (0.5*Struct.rho*Va*Va*Struct.S*a) + (0.5*Struct.rho*Struct.Sp*Struct.Cp*((Struct.k_mo*delta_t*Struct.k_mo*delta_t)-(Va*Va)));
        fy = (Struct.mass*Struct.g*cos(theta)*sin(phi))+(0.5*Struct.rho*Va*Va*Struct.S*j);
        fz = (Struct.mass*Struct.g*cos(theta)*cos(phi))+ (0.5*Struct.rho*Va*Va*Struct.S*o);

        %Moments
        d=Struct.b*(Struct.Cl0+(Struct.Clb*beta)+(0.5*Struct.Clp*Struct.b*p/Va)+(0.5*Struct.Clr*Struct.b*r/Va)+(Struct.Clda*delta_a)+(Struct.Cldr*delta_r));
        z=Struct.c*(Struct.Cm0+(Struct.Cma*alpha)+(0.5*Struct.Cmq*Struct.c*q/Va)+(Struct.Cmde*delta_e));
        f=Struct.b*(Struct.Cn0+(Struct.Cnb*beta)+(0.5*Struct.Cnp*Struct.b*p/Va)+(0.5*Struct.Cnr*Struct.b*r/Va)+(Struct.Cnda*delta_a)+(Struct.Cndr*delta_r));
        l=(0.5*Struct.rho*Va*Va*Struct.S*d) - (Struct.Ktp*Struct.Kw*delta_t*Struct.Kw*delta_t);
        m=(0.5*Struct.rho*Va*Va*Struct.S*z);
        n=(0.5*Struct.rho*Va*Va*Struct.S*f);
    end
    Outputs.omega=omega_actual;
    Outputs.euler=euler_angles_actual;
    Outputs.translate=rates_ned;
    Outputs.rotate=rates_euler;
end




