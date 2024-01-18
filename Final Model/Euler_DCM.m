function DCM = Euler_DCM(euler)
DCM=zeros([3*length(euler),3]);
k=1;
for i=1:length(euler)
    roll=euler(1,i);
    pitch=euler(2,i);
    yaw=euler(3,i);
    R_roll = [1 0 0; 0 cos(roll) -sin(roll); 0 sin(roll) cos(roll)];
    R_pitch = [cos(pitch) 0 sin(pitch); 0 1 0; -sin(pitch) 0 cos(pitch)];
    R_yaw = [cos(yaw) -sin(yaw) 0; sin(yaw) cos(yaw) 0; 0 0 1];
    temp=R_yaw * R_pitch * R_roll;
    normalized_dcm = zeros(3, 3);
    for l = 1:3
        column_vector = temp(:, l);
        normalized_dcm(:, l) = column_vector / norm(column_vector);
    end
    DCM(k:k+2,1:3)=normalized_dcm;
    k=k+3;
end



