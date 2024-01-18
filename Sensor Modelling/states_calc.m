function TRIM = states_calc(Va,y,R)
    alpha=0;
    beta=0;
    mu=(cos(y)*Va)^2/(9.8*R);
    deltat=0;
    deltae=0;
    deltaa=0;
    deltar=0;
    u=Va*cos(alpha)*cos(beta);
    v=Va*sin(beta);
    w=Va*sin(alpha)*cos(beta);
    stheta=(sin(y)*cos(alpha))+(cos(y)*cos(mu)*sin(alpha));
    theta=asin(stheta);
    cphi=((-sin(y)*sin(alpha))+(cos(y)*cos(mu)*cos(alpha)))/cos(theta);
    sphi=cos(y)*sin(mu)/cos(theta); 
    phi=acos(cphi);
    p=(Va*cos(y)/R)*(-sin(theta));
    q=(Va*cos(y)/R)*sin(phi)*cos(theta);
    r=(Va*cos(y)/R)*cos(phi)*cos(theta);

    %Inertial Properties - Aerosonde UAV
    mass=13.5;
    Jx=0.8244;
    Jy=1.135;
    Jz=1.759;
    Jxz=0.1204;

    % Parameters
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

    Cl=CL0+(Cla*alpha);
    Cd=CD0+(Cda*alpha);
    Cx=(-Cd*cos(alpha))+(Cl*sin(alpha));
    Cxq=(-Cdq*cos(alpha))+(Clq*sin(alpha));
    Cxde=(-Cdde*cos(alpha))+(Clde*sin(alpha));
    Cz=(-Cd*sin(alpha))-(Cl*cos(alpha));
    Czq=(-Cdq*sin(alpha))-(Clq*cos(alpha));
    Czde=(-Cdde*sin(alpha))+(Clde*cos(alpha));
    Cp0= G3*Cl0 + G4*Cn0;
    Cpb= G3*Clb + G4*Cnb;
    Cpp= G3*Clp + G4*Cnp;
    Cpr= G3*Clr + G4*Clr;
    Cpda= G3*Clda + G4*Cnda;
    Cpdr= G3*Cldr + G4*Cndr;
    Cr0= G4*Cl0 + G8*Cn0;
    Crb= G4*Clb + G8*Cnb;
    Crp= G4*Clp + G8*Cnp;
    Crr= G4*Clr + G8*Clr;
    Crda= G4*Clda + G8*Cnda;
    Crdr= G4*Cldr + G8*Cndr;

    %Model Co-efficients
    Yv = (((rho * S * b * v) / (4 * mass * Va)) * ...
    (Cyp * p + Cyr * r)) + ((rho * S * v / mass) * ...
    (Cy0 + Cyb * beta + Cyda * deltaa + Cydr * deltar)) +  ...
    ((rho * S * Cyb / (2 * mass)) * sqrt(u^2 + w^2));
    Yp = w + (rho * Va * S*b / (4 * mass)) * Cyp;
    Yr = -u + (rho * Va * S*b / (4 * mass)) * Cyr;
    Ydeltaa = ((rho * Va^2 * S) / (2 * mass)) * Cyda;
    Ydeltar = ((rho * Va^2 * S) / (2 * mass)) * Cydr;
    Lv = (((rho * S * b^2 * v) / (4 * Va)) * (Cpp * p + Cpr * r)) + ((rho * S * b * v) * (Cp0 + Cpb * beta + Cpda * deltaa + Cpdr * deltar))+ ((rho * S * b * Cpb) / (2 * sqrt(u^2 + w^2)));
    Lp = G1*q + ((rho * Va * S * b^2 * Cpp) /4);
    Lr = -G2*q + ((rho * Va * S * b^2 * Cpr) / 4);
    Ldeltaa = (rho * Va^2 * S * b * Cpda) /2;
    Ldeltar = (rho * Va^2 * S * b * Cpdr) /2;
    Nv = (((rho * S * b^2 * v) / (4 * Va)) * ...
    (Crp * p + Crr * r)) + ((rho * S * b * v) * ...
    (Cr0 + Crb * beta + Crda * deltaa + Crdr * deltar)) +  ...
    ((rho * S * b * Crb / 2) * sqrt(u^2 + w^2));
    Np = G7*q + ((rho * Va * S * b^2 * Crp) /4);
    Nr = -G1*q + ((rho * Va * S * b^2 * Crr) / 4);
    Ndeltaa = (rho * Va^2 * S * b * Crda) /2;
    Ndeltar = (rho * Va^2 * S * b * Crdr) /2;
    Xu = ((u * rho * S /mass) * (-Cda+ Cla*alpha+ Cxde * deltae)) ...
    - (rho * S * w * Cx / (2 * mass)) + ((rho * S * c * Cxq * u * q) / (4 * mass * Va)) ...
    - (rho * Sp * Cp * u/mass);
    Xw = -q + ((w * rho * S / mass) * (-Cda+ Cla*alpha + Cxde * deltae)) ...
    + (rho * S * c * Cxq * w * q / (4 * mass * Va)) ...
    + (rho * S * Cla * u / (2 * mass)) - (rho * Sp * Cp* w / mass);
    Xq = -w + ((rho * Va * S * Cxq * c) / (4 * mass));
    Xdeltae = (rho * Va^2 * S * Cxde) / (2 * mass);
    Xdeltat = (rho * Sp * Cp * k_mo^2 * deltat) / mass;
    Zu = q+ ((u * rho * S/mass)* (-Cla -Cda*alpha+ Czde * deltae)) ...
    - (rho * S * -Cda * w / (2 * mass)) ...
    + (u* rho * S * Czq * c * q / (4 * mass * Va));
    Zw = ((w * rho * S/mass)* (-Cla -Cda*alpha + Czde * deltae)) ...
    - (rho * S * Cz * u / (2 * mass)) ...
    + (w* rho * S * Czq * c * q / (4 * mass * Va));
    Zq = u + (rho * Va * S * Czq * c) / (4 * mass);
    Zdeltae = (rho * Va^2 * S * Czde) / (2 * mass);
    Mu = ((u * rho * S* c/Jy)* (Cm0 + Cma * alpha+ Cmde * deltae)) ...
    - (rho * S * c * Cma * w / (2 * Jy)) ...
    + (u* rho * S * Cmq * c^2 * q / (4 * Jy * Va));
    Mw = ((w * rho * S*c/mass)* (Cm0 + Cma * alpha+ Cmde * deltae)) ...
    - (rho * S * c * Cma * u / (2 * Jy)) ...
    + (w* rho * S * Cmq * c^2 * q / (4 * Jy * Va));
    Mq = (rho * Va * S * Cmq * c^2) / (4 * Jy);
    Mdeltae = (rho * Va^2 * S * c * Cmde) / (2 * Jy);

    %State Space Model
    % State vector
    X = [beta; p; r; phi;u; alpha; q; theta];
    % Control inputs
    U = [deltaa; deltar;deltae; deltat];

    % Coefficient matrices
    A = [Yv, Yp/(Va*cos(beta)), Yr/(Va*cos(beta)), (g*cos(theta)*cos(phi)/(Va*cos(beta))),0,0,0,0;
     Lv*Va*cos(beta), Lp, Lr, 0, 0,0,0,0;
     Nv*Va*cos(beta), Np, Nr, 0, 0,0,0,0;
     0, 1, cos(phi)*tan(theta), (q*cos(phi)*tan(theta)-r*sin(phi)*tan(theta)),0,0,0,0;
     0,0,0,0,Xu, Xw * Va * cos(alpha), Xq , -g * cos(theta);
     0,0,0,0,Zu/(Va * cos(alpha)),Zw, Zq/(Va * cos(alpha)), -g * sin(theta)/(Va * cos(alpha));
     0,0,0,0,Mu, Mw * Va * cos(alpha), Mq, 0;
     0,0,0,0,0, 0, 1, 0];

    B = [Ydeltaa/(Va*cos(beta)), Ydeltar/(Va*cos(beta)),0,0;
     Ldeltaa, Ldeltar,0,0;
     Ndeltaa, Ndeltar,0,0;
     0, 0,0,0;
     0,0, Xdeltae, Xdeltat;
     0,0,Zdeltae/(Va * cos(alpha)), 0;
     0,0,Mdeltae, 0;
     0,0,0, 0];

    % State derivative (dx/dt)
    Xdot = A * X + B * U;

    % Calculating Trim States - Using Gradient Descent
    % Initialize state variables
    X = [beta; p; r; phi; u; alpha; q; theta]; % Initial guess for states
    learning_rate = 0.01; % Learning rate
    num_iterations = 1000; % Number of iterations

    % Define desired trim conditions for X_dot (you may replace this with your criteria)
    desired_X_dot = zeros(size(X));

    % Gradient descent iterations
    for iteration = 1:num_iterations
        % Compute the cost function
        X_dot_predicted = A * X + B * U; % Predicted X_dot
        cost = norm(Xdot - X_dot_predicted)^2; % Squared error
    
        % Compute the gradient of the cost function with respect to X
        gradient = -2 * A' * (Xdot - X_dot_predicted); % Gradient
    
        % Update the state variables using gradient descent
        X = X - learning_rate * gradient;
    
        % Check for convergence (you can define your own convergence criteria)
        if cost < 0.001
            break;
        end
    end

    % The final values of the state variables are the trimmed states
    trimmed_states = X;

    % Calculating Trimmed States
    beta_trim=trimmed_states(1);
    p_trim=trimmed_states(2);
    r_trim=trimmed_states(3);
    phi_trim=trimmed_states(4);
    u_trim=trimmed_states(5);
    alpha_trim=trimmed_states(6);
    q_trim=trimmed_states(7);
    theta_trim=trimmed_states(8);
    v_trim=Va*sin(beta_trim);
    w_trim=Va*sin(alpha_trim)*cos(beta_trim);
    dphi_trim=0;
    dtheta_trim=0;
    dpsi_trim=Va*cos(y)/R;

    % Elevator
    temp=((Jxz*(p_trim*p_trim-r_trim*r_trim)) + ((Jx-Jz)*p_trim*r_trim))/(0.5*rho*Va*Va*c*S);
    deltae_trim=(temp - Cm0 - (Cma*alpha_trim) - (0.5*Cmq*c*q_trim/Va))/Cmde;

    % Throttle
    temp=(2*mass*(-r_trim*v_trim + q_trim*w_trim + g*sin(theta_trim))) - ((rho*Va*Va*S)*(Cx*alpha_trim+(0.5*Cxq*(alpha_trim)*c*q_trim/Va)+(Cxde*(alpha_trim)*deltae_trim)));
    deltat_trim=sqrt((temp/(rho*Sp*Cp*k_mo*k_mo)) + ((Va*Va)/(k_mo*k_mo)));

    C_mat=[Cpda Cpdr; Crda Crdr];
    temp1=(-G1*p_trim*q_trim + G2*q_trim*r_trim)/(0.5*rho*Va*Va*S*b);
    temp2=(-G7*p_trim*q_trim + G1*q_trim*r_trim)/(0.5*rho*Va*Va*S*b);
    temp=inv(C_mat)*[temp1-Cp0-Cpb*beta_trim-(0.5*Cpp*b*p_trim/Va)-(0.5*Cpr*b*r_trim/Va);temp2-Cr0-Crb*beta_trim-(0.5*Crp*b*p_trim/Va)-(0.5*Crr*b*r_trim/Va)];
    deltaa_trim=temp(1);
    deltar_trim=temp(2);
    TRIM(1)=deltaa_trim;
    TRIM(2)=deltae_trim;
    TRIM(3)=deltar_trim;
    TRIM(4)=deltat_trim;
    TRIM(5)=u_trim;
    TRIM(6)=v_trim;
    TRIM(7)=w_trim;
end

