function qinv = qinverse(q)
    qinv(1,1)=q(1,1);
    qinv(2,1)=-q(2,1);
    qinv(3,1)=-q(3,1);
    qinv(4,1)=-q(4,1);
end