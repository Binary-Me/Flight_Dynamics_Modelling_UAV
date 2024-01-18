function [actual,sensor] = Tracker(DCM)
    Home=Home_Cordinates();
    actual.h=zeros([3,length(DCM)/3]);
    actual.k=zeros([3,length(DCM)/3]);
    j=1;
    for i=1:length(DCM)/3
        actual.h(:,i)=DCM(j:j+2,1:3)*Home.h;
        actual.k(:,i)=DCM(j:j+2,1:3)*Home.k;
        j=j+3;
    end

    sensor.h=zeros([3,length(DCM)/3]);
    sensor.k=zeros([3,length(DCM)/3]);

    variance=Sensor_variance_inputs();
    constant_noise_k=[variance.k1; variance.k2; variance.k3].*(2*rand([3,1])-[1;1;1]);
    constant_noise_h=[variance.h1; variance.h2; variance.h3].*(2*rand([3,1])-[1;1;1]);

    for i=1:length(DCM)/3
        jitter_h=[variance.h1; variance.h2; variance.h3].*(2*rand([3,1])-[1;1;1]);
        jitter_k=[variance.k1; variance.k2; variance.k3].*(2*rand([3,1])-[1;1;1]);
        sensor.h(:,i)=actual.h(:,i)+constant_noise_h+jitter_h;
        sensor.k(:,i)=actual.k(:,i)+constant_noise_k+jitter_k;
    end


end

