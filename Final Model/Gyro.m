function gyro= Gyro(omega)
    variance=Sensor_variance_inputs();
    constant_noise=[variance.w1; variance.w2; variance.w3].*(2*rand([3,1])-[1;1;1]);
    for i=1:length(omega)
        jitter=[variance.w1;variance.w2;variance.w3].*(2*rand([3,1])-[1;1;1]);
        gyro.omega(:,i)=omega(:,i)+constant_noise+jitter;
    end
end