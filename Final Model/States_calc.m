function TRIM = States_calc(Va,y,R,Struct)
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

    Cl=Struct.CL0+(Struct.Cla*alpha);
    Cd=Struct.CD0+(Struct.Cda*alpha);
    Cx=(-Cd*cos(alpha))+(Cl*sin(alpha));
    Cxq=(-Struct.Cdq*cos(alpha))+(Struct.Clq*sin(alpha));
    Cxde=(-Struct.Cdde*cos(alpha))+(Struct.Clde*sin(alpha));
    Cz=(-Cd*sin(alpha))-(Cl*cos(alpha));
    Czq=(-Struct.Cdq*sin(alpha))-(Struct.Clq*cos(alpha));
    Czde=(-Struct.Cdde*sin(alpha))+(Struct.Clde*cos(alpha));
    
    Cp0= Struct.G3*Struct.Cl0 + Struct.G4*Struct.Cn0;
    Cpb= Struct.G3*Struct.Clb + Struct.G4*Struct.Cnb;
    Cpp= Struct.G3*Struct.Clp + Struct.G4*Struct.Cnp;
    Cpr= Struct.G3*Struct.Clr + Struct.G4*Struct.Cnr;
    Cpda= Struct.G3*Struct.Clda + Struct.G4*Struct.Cnda;
    Cpdr= Struct.G3*Struct.Cldr + Struct.G4*Struct.Cndr;
    Cr0= Struct.G4*Struct.Cl0 + Struct.G8*Struct.Cn0;
    Crb= Struct.G4*Struct.Clb + Struct.G8*Struct.Cnb;
    Crp= Struct.G4*Struct.Clp + Struct.G8*Struct.Cnp;
    Crr= Struct.G4*Struct.Clr + Struct.G8*Struct.Cnr;
    Crda= Struct.G4*Struct.Clda + Struct.G8*Struct.Cnda;
    Crdr= Struct.G4*Struct.Cldr + Struct.G8*Struct.Cndr;

    %Model Co-efficients
    Yv = (((Struct.rho * Struct.S * Struct.b * v) / (4 * Struct.mass * Va)) * ...
    (Struct.Cyp * p + Struct.Cyr * r)) + ((Struct.rho * Struct.S * v / Struct.mass) * ...
    (Struct.Cy0 + Struct.Cyb * beta + Struct.Cyda * deltaa + Struct.Cydr * deltar)) +  ...
    ((Struct.rho * Struct.S * v * Struct.Cyb / (2 * Struct.mass)) * sqrt(u^2 + w^2));
    Yp = w + (Struct.rho * Va * Struct.S*Struct.b / (4 * Struct.mass)) * Struct.Cyp;
    Yr = -u + (Struct.rho * Va * Struct.S*Struct.b / (4 * Struct.mass)) * Struct.Cyr;
    Ydeltaa = ((Struct.rho * Va^2 * Struct.S) / (2 * Struct.mass)) * Struct.Cyda;
    Ydeltar = ((Struct.rho * Va^2 * Struct.S) / (2 * Struct.mass)) * Struct.Cydr;
    Lv = (((Struct.rho * Struct.S * Struct.b^2 * v) / (4 * Va)) * (Cpp * p + Cpr * r)) + ((Struct.rho * Struct.S * Struct.b * v) * (Cp0 + Cpb * beta + Cpda * deltaa + Cpdr * deltar))+ ((Struct.rho * Struct.S * Struct.b * Cpb) / (2 * sqrt(u^2 + w^2)));
    Lp = Struct.G1*q + ((Struct.rho * Va * Struct.S * Struct.b^2 * Cpp) /4);
    Lr = -Struct.G2*q + ((Struct.rho * Va * Struct.S * Struct.b^2 * Cpr) / 4);
    Ldeltaa = (Struct.rho * Va^2 * Struct.S * Struct.b * Cpda) /2;
    Ldeltar = (Struct.rho * Va^2 * Struct.S * Struct.b * Cpdr) /2;
    Nv = (((Struct.rho * Struct.S * Struct.b^2 * v) / (4 * Va)) * ...
    (Crp * p + Crr * r)) + ((Struct.rho * Struct.S * Struct.b * v) * ...
    (Cr0 + Crb * beta + Crda * deltaa + Crdr * deltar)) +  ...
    ((Struct.rho * Struct.S * Struct.b * Crb / 2) * sqrt(u^2 + w^2));
    Np = Struct.G7*q + ((Struct.rho * Va * Struct.S * Struct.b^2 * Crp) /4);
    Nr = -Struct.G1*q + ((Struct.rho * Va * Struct.S * Struct.b^2 * Crr) / 4);
    Ndeltaa = (Struct.rho * Va^2 * Struct.S * Struct.b * Crda) /2;
    Ndeltar = (Struct.rho * Va^2 * Struct.S * Struct.b * Crdr) /2;
    Xu = ((u * Struct.rho * Struct.S /Struct.mass) * (-Struct.Cda+ Struct.Cla*alpha+ Cxde * deltae)) ...
    - (Struct.rho * Struct.S * w * Cx / (2 * Struct.mass)) + ((Struct.rho * Struct.S * Struct.c * Cxq * u * q) / (4 * Struct.mass * Va)) ...
    - (Struct.rho * Struct.Sp * Struct.Cp * u/Struct.mass);
    Xw = -q + ((w * Struct.rho * Struct.S / Struct.mass) * (-Struct.Cda+ Struct.Cla*alpha + Cxde * deltae)) ...
    + (Struct.rho * Struct.S * Struct.c * Cxq * w * q / (4 * Struct.mass * Va)) ...
    + (Struct.rho * Struct.S * Struct.Cla * u / (2 * Struct.mass)) - (Struct.rho * Struct.Sp * Struct.Cp* w / Struct.mass);
    Xq = -w + ((Struct.rho * Va * Struct.S * Cxq * Struct.c) / (4 * Struct.mass));
    Xdeltae = (Struct.rho * Va^2 * Struct.S * Cxde) / (2 * Struct.mass);
    Xdeltat = (Struct.rho * Struct.Sp * Struct.Cp * Struct.k_mo^2 * deltat) / Struct.mass;
    Zu = q+ ((u * Struct.rho * Struct.S/Struct.mass)* (-Struct.Cla -Struct.Cda*alpha+ Czde * deltae)) ...
    - (Struct.rho * Struct.S * -Struct.Cda * w / (2 * Struct.mass)) ...
    + (u* Struct.rho * Struct.S * Czq * Struct.c * q / (4 * Struct.mass * Va));
    Zw = ((w * Struct.rho * Struct.S/Struct.mass)* (-Struct.Cla -Struct.Cda*alpha + Czde * deltae)) ...
    - (Struct.rho * Struct.S * Cz * u / (2 * Struct.mass)) ...
    + (w* Struct.rho * Struct.S * Czq * Struct.c * q / (4 * Struct.mass * Va));
    Zq = u + (Struct.rho * Va * Struct.S * Czq * Struct.c) / (4 * Struct.mass);
    Zdeltae = (Struct.rho * Va^2 * Struct.S * Czde) / (2 * Struct.mass);
    Mu = ((u * Struct.rho * Struct.S* Struct.c/Struct.Jy)* (Struct.Cm0 + Struct.Cma * alpha+ Struct.Cmde * deltae)) ...
    - (Struct.rho * Struct.S * Struct.c * Struct.Cma * w / (2 * Struct.Jy)) ...
    + (u* Struct.rho * Struct.S * Struct.Cmq * Struct.c^2 * q / (4 * Struct.Jy * Va));
    Mw = ((w * Struct.rho * Struct.S*Struct.c/Struct.mass)* (Struct.Cm0 + Struct.Cma * alpha+ Struct.Cmde * deltae)) ...
    - (Struct.rho * Struct.S * Struct.c * Struct.Cma * u / (2 * Struct.Jy)) ...
    + (w* Struct.rho * Struct.S * Struct.Cmq * Struct.c^2 * q / (4 * Struct.Jy * Va));
    Mq = (Struct.rho * Va * Struct.S * Struct.Cmq * Struct.c^2) / (4 * Struct.Jy);
    Mdeltae = (Struct.rho * Va^2 * Struct.S * Struct.c * Struct.Cmde) / (2 * Struct.Jy);


    %State Space Model
    % State vector
    X = [beta; p; r; phi;u; alpha; q; theta];
    % Control inputs
    U = [deltaa; deltar;deltae; deltat];

    % Coefficient matrices
    A = [Yv, Yp/(Va*cos(beta)), Yr/(Va*cos(beta)), (Struct.g*cos(theta)*cos(phi)/(Va*cos(beta))),0,0,0,0;
     Lv*Va*cos(beta), Lp, Lr, 0, 0,0,0,0;
     Nv*Va*cos(beta), Np, Nr, 0, 0,0,0,0;
     0, 1, cos(phi)*tan(theta), (q*cos(phi)*tan(theta)-r*sin(phi)*tan(theta)),0,0,0,0;
     0,0,0,0,Xu, Xw * Va * cos(alpha), Xq , -Struct.g * cos(theta);
     0,0,0,0,Zu/(Va * cos(alpha)),Zw, Zq/(Va * cos(alpha)), -Struct.g * sin(theta)/(Va * cos(alpha));
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
    temp=((Struct.Jxz*(p_trim*p_trim-r_trim*r_trim)) + ((Struct.Jx-Struct.Jz)*p_trim*r_trim))/(0.5*Struct.rho*Va*Va*Struct.c*Struct.S);
    deltae_trim=(temp - Struct.Cm0 - (Struct.Cma*alpha_trim) - (0.5*Struct.Cmq*Struct.c*q_trim/Va))/Struct.Cmde;

    % Throttle
    temp=(2*Struct.mass*(-r_trim*v_trim + q_trim*w_trim + Struct.g*sin(theta_trim))) - ((Struct.rho*Va*Va*Struct.S)*(Cx*alpha_trim+(0.5*Cxq*(alpha_trim)*Struct.c*q_trim/Va)+(Cxde*(alpha_trim)*deltae_trim)));
    deltat_trim=sqrt((temp/(Struct.rho*Struct.Sp*Struct.Cp*Struct.k_mo*Struct.k_mo)) + ((Va*Va)/(Struct.k_mo*Struct.k_mo)));

    C_mat=[Cpda Cpdr; Crda Crdr];
    temp1=(-Struct.G1*p_trim*q_trim + Struct.G2*q_trim*r_trim)/(0.5*Struct.rho*Va*Va*Struct.S*Struct.b);
    temp2=(-Struct.G7*p_trim*q_trim + Struct.G1*q_trim*r_trim)/(0.5*Struct.rho*Va*Va*Struct.S*Struct.b);
    temp=C_mat\[temp1-Cp0-Cpb*beta_trim-(0.5*Cpp*Struct.b*p_trim/Va)-(0.5*Cpr*Struct.b*r_trim/Va);temp2-Cr0-Crb*beta_trim-(0.5*Crp*Struct.b*p_trim/Va)-(0.5*Crr*Struct.b*r_trim/Va)];
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


