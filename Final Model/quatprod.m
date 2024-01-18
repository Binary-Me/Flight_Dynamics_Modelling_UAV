%To find the product of two quaternions
function prod = quatprod(p,q)
    p0=p(1)*q(1)-p(2)*q(2)-p(3)*q(3)-p(4)*q(4);
    p1=p(1)*q(2)+p(2)*q(1)+p(3)*q(4)-p(4)*q(3);
    p2=p(1)*q(3)-p(2)*q(4)+p(3)*q(1)+p(4)*q(2);
    p3=p(1)*q(4)+p(2)*q(3)-p(3)*q(2)+p(4)*q(1);
    prod=[p0;p1;p2;p3];
end