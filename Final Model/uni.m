%Obtaining a Unit Vector
function a = uni(b)
norm=0;
for i=1:3
    norm=norm+(b(i)^2);
end
norm=sqrt(norm);
a=b./norm;
end