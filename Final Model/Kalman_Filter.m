function Q_kalman = Kalman_Filter(w,sensor)
    Q_kalman=zeros(4,length(w));
    Q_kalman(:,1)=[1;0;0;0];
    r=zeros(4,length(w));
    r(:,1)=[1;0;0;0];
    R=zeros(3,3,length(w));

    variance=Sensor_variance_inputs();
    E=diag(sqrt(2).*[variance.w1,variance.w2,variance.w3]);
    O=diag(sqrt(2).*[variance.h1,variance.h2,variance.h3,variance.k1,variance.k2,variance.k3]);

    T=Time_inputs();
    total_time=T.Time;
    dt=T.dt;

    Home=Home_Cordinates();
    h0=Home.h;
    k0=Home.k;

    for i=1:(total_time/dt)
        phi=mag(w(:,i))*dt/2;
        temp=quatprod(Q_kalman(:,i),[cos(phi); uni(w(:,i)).*sin(phi)]);
        A=eye(3) - (crossprodv(w(:,i))).*dt;
        R_temp=A*R(:,:,i)*A' + E*dt^2./4; %Predict Step
        s=qinverse(temp);
        z_temp=[quatprod(quatprod(s,[0;h0]),temp); quatprod(quatprod(s,[0;k0]),temp)];
        u_temp=[z_temp(2:4);z_temp(6:8)];
        C=2.*[crossprodv(u_temp(1:3));crossprodv(u_temp(4:6))];
        L=R_temp*C'/(C*R_temp*C'+O);  %Update Step
        R(:,:,i+1)=R_temp-L*C*R_temp;
        r(2:4,i+1)=L*[sensor.h(:,i+1)-u_temp(1:3); sensor.k(:,i+1)-u_temp(4:6)];
        r(1,i+1)=sqrt(1-(mag(r(2:4,i+1)))^2);
        Q_kalman(:,i+1)=quatprod(temp, r(:,i+1));
    end
end

