no=length(T);
dt=0.01;
Q_kalman=zeros(4,no);
Q_kalman(:,1)=[1;0;0;0];
r=zeros(4,no);
r(:,1)=[1;0;0;0];
R=zeros(3,3,no);

E=0.01.*eye(3);
O=0.01.*eye(6);

h0=[1;0;0];
k0=[0;1;0];

w=w_sensor;

for i=1:(no-1)
    phi=mag(w(:,i))*dt/2;
    temp=Norm(quatprod(Q_kalman(:,i),[cos(phi); uni(w(:,i)).*sin(phi)]));
    A=eye(3) - (crossprodv(w(:,i))).*dt;
    R_temp=A*R(:,:,i)*A' + E*dt^2./4; %Predict Step
    s=qinverse(temp);
    z_temp=[Norm(quatprod(Norm(quatprod(s,[0;h0])),temp)); Norm(quatprod(Norm(quatprod(s,[0;k0])),temp))];
    u_temp=[z_temp(2:4);z_temp(6:8)];
    C=2.*[crossprodv(u_temp(1:3));crossprodv(u_temp(4:6))];
    L=R_temp*C'/(C*R_temp*C'+O);  %Update Step
    R(:,:,i+1)=R_temp-L*C*R_temp;
    r(2:4,i+1)=L*[sensor_h(:,i+1)-u_temp(1:3); sensor_k(:,i+1)-u_temp(4:6)];
    r(1,i+1)=sqrt(1-(mag(r(2:4,i+1)))^2);
    Q_kalman(:,i+1)=Norm(quatprod(temp, r(:,i+1)));
end

euler_est=zeros([3,no]);
for i=1:length(Q_kalman)
    q0 = real(Q_kalman(1,i));
    q1 = real(Q_kalman(2,i));
    q2 = real(Q_kalman(3,i));
    q3 = real(Q_kalman(4,i));
    
    % Calculate roll (phi)
    phi = atan2(2*(q0*q1 + q2*q3), 1 - 2*(q1^2 + q2^2));

    % Calculate pitch (theta)
    theta = asin(2*(q0*q2 - q3*q1));

    % Calculate yaw (psi)
    psi = atan2(2*(q0*q3 + q1*q2), 1 - 2*(q2^2 + q3^2));

    % Output Euler angles in radians
    euler_est(:,i) = [phi; theta; psi];
end

%%
figure
angles=real(angles);
subplot(3, 1, 1);
plot(T,mod(rad2deg(angles(1,:))+180,360)-180);
hold on
plot(T,mod(rad2deg(euler_est(1,:))+180,360)-180);
legend("\phi",'\phi_{est}');
xlabel("Time (in sec)");
ylabel("\phi in degrees");

subplot(3, 1, 2);
plot(T,mod(rad2deg(angles(2,:))+180,360)-180);
hold on
plot(T,mod(rad2deg(euler_est(2,:))+180,360)-180);
legend("\theta",'\theta_{est}');
xlabel("Time (in sec)");
ylabel("\theta in degrees");

subplot(3, 1, 3);
plot(T,mod(rad2deg(angles(3,:))+180,360)-180);
hold on
plot(T,mod(rad2deg(euler_est(3,:))+180,360)-180);
legend("\psi",'\psi_{est}');
xlabel("Time (in sec)");
ylabel("\psi in degrees");

sgtitle("AKF - Euler Angles - Steady Coordinated Turn")

%%
function qinv = qinverse(q)
    qinv(1,1)=q(1,1);
    qinv(2,1)=-q(2,1);
    qinv(3,1)=-q(3,1);
    qinv(4,1)=-q(4,1);
end

function V = crossprodv(u)
    V(1,1)=0;
    V(1,2)=-u(3);
    V(1,3)=u(2);
    V(2,1)=u(3);
    V(2,2)=0;
    V(2,3)=-u(1);
    V(3,1)=-u(2);
    V(3,2)=u(1);
    V(3,3)=0;
end

function a = mag(b)
    norm=0;
    for i=1:3
        norm=norm+(b(i)^2);
    end
    norm=sqrt(norm);
    a=norm;
end


function a = uni(b)
    norm=0;
    for i=1:3
        norm=norm+(b(i)^2);
    end
    norm=sqrt(norm);
    a=b./norm;
end

function prod = quatprod(p,q)
    p0=p(1)*q(1)-p(2)*q(2)-p(3)*q(3)-p(4)*q(4);
    p1=p(1)*q(2)+p(2)*q(1)+p(3)*q(4)-p(4)*q(3);
    p2=p(1)*q(3)-p(2)*q(4)+p(3)*q(1)+p(4)*q(2);
    p3=p(1)*q(4)+p(2)*q(3)-p(3)*q(2)+p(4)*q(1);
    prod=[p0;p1;p2;p3];
end

function n=Norm(q)
    sum=sqrt(q(1)^2+q(2)^2+q(3)^2+q(4)^2);
    n=q./sum;
end