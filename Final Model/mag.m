%Obtaining Magnitude of a Vector
function a = mag(b)
norm=0;
for i=1:3
    norm=norm+(b(i)^2);
end
norm=sqrt(norm);
a=norm;
end