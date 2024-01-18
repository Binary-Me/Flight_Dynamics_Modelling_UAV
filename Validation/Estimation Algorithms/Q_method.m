%% Daven ports q method

% properties of satellite
h0 = [1 0 0]';
k0 = [0 1 0]';

q_from_Davenport=zeros([4,no]);
euler_est=zeros([3,no]);

%Davenports
for i = 1:no
    [q_from_Davenport(:, i)] = DavenportsQ(sensor_h(:, i), sensor_k(:, i), h0, k0);
end

%for i = 2:1000
    %[q_from_Davenport(:, i)] = correct(q_from_Davenport(:, i), q_from_Davenport(:, i - 1));
%end
%%
for i=1:length(q_from_Davenport)
    q0 = real(q_from_Davenport(1,i));
    q1 = real(q_from_Davenport(2,i));
    q2 = real(q_from_Davenport(3,i));
    q3 = real(q_from_Davenport(4,i));
    
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
plot(T,mod(rad2deg(euler_est(1,:))+180,360)-180,':');
legend("\phi",'\phi_{est}');
xlabel("Time (in sec)");
ylabel("\phi in degrees");

subplot(3, 1, 2);
plot(T,mod(rad2deg(angles(2,:))+180,360)-180);
hold on
plot(T,mod(rad2deg(euler_est(2,:))+180,360)-180,':');
legend("\theta",'\theta_{est}');
xlabel("Time (in sec)");
ylabel("\theta in degrees");

subplot(3, 1, 3);
plot(T,mod(rad2deg(angles(3,:))+180,360)-180);
hold on
plot(T,mod(rad2deg(euler_est(3,:))+180,360)-180,":");
legend("\psi",'\psi_{est}');
xlabel("Time (in sec)");
ylabel("\psi in degrees");

sgtitle("Davenport - Euler Angles - Steady Coordinated Turn")
%%
function [Q] = DavenportsQ(u, v, g, h)
    a = 1/2;
    b = 1/2;
    

    D = a * (Cr(u)) * (Cnr(g)) + b * (Cr(v)) * (Cnr(h));
    [q, e] = eig(D);

    max = 0;
    for i = 1:4
        if (e(i, i) > max)
            max = e(i, i);
            pos = i;
        end
    end

    Q(1, 1) = q(1, i);
    Q(2, 1) = q(2, i);
    Q(3, 1) = q(3, i);
    Q(4, 1) = q(4, i);

end

function [outputArg1] = correct(inputArg1, inputArg2)

    for i = 1:4
        if (abs(inputArg1(i, 1) - inputArg2(i, 1)) > 0.3)
            outputArg1 = -1 * inputArg1;
            break;
        else
            outputArg1 = inputArg1;
        end
    end
end

function [CR] = Cr(u)

    CR = [0 u'; -u mat(u)];

end

function [CR] = Cnr(u)
    CR = [0 -u'; u mat(u)];
end

function [S] = mat(u)
    S = [0 -u(3, 1) u(2, 1); u(3, 1) 0 -u(1, 1); -u(2, 1) u(1, 1) 0];
end

