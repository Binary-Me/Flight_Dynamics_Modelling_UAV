%%Interface
disp("Enter the ambient speed: ");
input_str = input('Va: ', 's');
Va= str2double(input_str);
disp("Enter the heading angle in degrees: ");
input_str = input('y: ', 's');
y= str2double(input_str);
y=y*pi/180;
disp("Enter the turn Radius: ");
input_str = input('R: ', 's');
turnR= str2double(input_str);

Output=Dynamic_Model(Va,y,turnR);
Animate(Output.rotate,Output.translate);
gyro=Gyro(Output.omega);
DCM_actual=Euler_DCM(Output.euler);
[actual,sensor] = Tracker(DCM_actual);
estimate_euler=Estimate_Euler(gyro.omega,sensor);
%%
angle_est=mod(rad2deg(real(estimate_euler(1,:))),360);
angle_est(angle_est > 180) = angle_est(angle_est > 180) - 360;
angle_act=mod(rad2deg(real(Output.euler(1,:))),360);
angle_act(angle_act > 180) = angle_act(angle_act > 180) - 360;
figure
plot(0:0.001:200,-(angle_est),'r');
hold on
plot(0:0.001:200,(angle_act),'b');
legend("estimated","actual");
title("Phi vs T")
xlabel(" Time in sec")
ylabel("Phi in degrees")