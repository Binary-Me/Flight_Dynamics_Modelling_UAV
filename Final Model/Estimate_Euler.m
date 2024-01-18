function euler = Estimate_Euler(w,sensor)
    Q_kalman = Kalman_Filter(w,sensor);
    euler=zeros([3,length(w)]);
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
        euler(:,i) = [phi; theta; psi];
    end
end

