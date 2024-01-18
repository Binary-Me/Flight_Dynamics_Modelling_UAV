function n=Norm(q)
    sum=sqrt(q(1)^2+q(2)^2+q(3)^2+q(4)^2);
    n=q./sum;
end