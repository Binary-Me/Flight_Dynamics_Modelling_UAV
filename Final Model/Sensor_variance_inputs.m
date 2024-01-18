function variance = Sensor_variance_inputs()
    %Variance in sensor measurement
    %Tracker
    variance.h1=0.0001;
    variance.h2=0.0001;
    variance.h3=0.0001;
    variance.k1=0.0001;
    variance.k2=0.0001;
    variance.k3=0.0001;
    %Gyroscope
    variance.w1=0.001;
    variance.w2=0.001;
    variance.w3=0.001;
end

