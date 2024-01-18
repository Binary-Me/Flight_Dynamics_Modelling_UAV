function q=Matrix_to_Quat(C)
    q(1)=0.25*sqrt(1+C(1,1)+C(2,2)+C(3,3));
    q(2)=(C(3,2) - C(2,3))/4/q(1);
    q(3)=(C(1,3) - C(3,1))/4/q(1);
    q(4)=(C(2,1) - C(1,2))/4/q(1);
end