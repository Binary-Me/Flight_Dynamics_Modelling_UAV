mu_values=zeros([1,31,4]);
R_values=[100;500;1000;1500];
for k=1:4
    R=R_values(k);
    j=1;
    for i=5:1:35
        % Trim inputs
        Va=i;
        y=(pi/180)*0;
        
        % Trim states initial guess
        
        alpha=0;
        beta=0;
        mu=(cos(y)*Va)^2/(9.8*R);
        del_t=0;
        del_e=0;
        del_a=0;
        del_r=0;
        
        X_guess=[alpha;beta;mu;del_e;del_t;del_a;del_r];
        options=optimoptions('fsolve','TolX',1e-6);
        trm_state=fsolve(@(X)sct(X,Va,y,R),X_guess,options);
        
        alpha=trm_state(1)*180/pi;
        beta=trm_state(2)*180/pi;
        mu=mod(rad2deg(trm_state(3)),360);
        del_e=trm_state(4);
        del_t=trm_state(5);
        del_a=trm_state(6);
        del_r=trm_state(7);
        
        T=table(alpha,beta,mu,del_e,del_t,del_a,del_r);
        mu_values(1,j,k)=mu;
        j=j+1;
    end
end
%%
figure
plot(5:1:35,mu_values(:,:,1));
hold on
plot(5:1:35,mu_values(:,:,2));
plot(5:1:35,mu_values(:,:,3));
plot(5:1:35,mu_values(:,:,4));
xlabel('V_a')
ylabel("\mu")
legend("R=100m","R=500m","R=1000m","R=2000m");
title("\mu vs V_a")

%%

function residue=sct(X,Va,y,R)  % sct~ steady coordinated turn

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    alpha=X(1);
    beta=X(2);
    mu=X(3);
    del_e=X(4);
    del_t=X(5);
    del_a=X(6);
    del_r=X(7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rho=1.2682;
    S=0.55;
    m=13.5;
    
    Jx=0.8244;
    Jy=1.135;
    Jz=1.759;
    Jxz=0.1204;

    b=2.8956;
    c=0.18994;
    Sprop=0.2027;
    Kmotor=80;
    KTp=0;
    Komg=0;
    e=0.9;

    CLo=0.28;
    CDo=0.03;
    Cmo=-0.02338;
    CLalpha=3.45;
    CDalpha=0.30;
    Cmalpha=-0.38;
    CLq=0;
    CDq=0;
    Cmq=-3.6;
    CLdele=-0.36;
    CDdele=0;
    Cmdele=-0.5;
    Cprop=1;
    M=50;
    alphao=0.4712;
    epsilon=0.1592;
    CDp=0.0437;
    Cndelr=-0.032;

    CYo=0;
    Clo=0;
    Cno=0;
    CYbeta=-0.98;
    Clbeta=-0.12;
    Cnbeta=0.25;
    CYp=0;
    Clp=-0.26;
    Cnp=0.022;
    CYr=0;
    Clr=0.14;
    Cnr=-0.35;
    CYdela=0;
    Cldela=0.08;
    Cndela=0.06;
    CYdelr=-0.17;
    Cldelr=0.105;

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    g=(Jx*Jz)-(Jxz^2);
    g1=Jxz*(Jx-Jy+Jz)/g;
    g2=((Jz*(Jz-Jy))+Jxz^2)/g;
    g3=Jz/g;
    g4=Jxz/g;
    g5=(Jz-Jx)/Jy;
    g6=Jxz/Jy;
    g7=((Jx*(Jx-Jy))+Jxz^2)/g;
    g8=Jx/g;
    g9=1/Jy;
    
    
    
    CL=CLo+(CLalpha*alpha);
    CD=CDo+(CDalpha*alpha);

    CX=(-CD*cos(alpha))+(CL*sin(alpha));
    CXq=(-CDq*cos(alpha))+(CLq*sin(alpha));
    CXdele=(-CDdele*cos(alpha))+(CLdele*sin(alpha));
    CZ=(-CD*sin(alpha))-(CL*cos(alpha));
    CZq=(-CDq*sin(alpha))-(CLq*cos(alpha));
    CZdele=(-CDdele*sin(alpha))-(CLdele*cos(alpha));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% geometric constants

    % gamma 1-8

  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Cpo=(g3*Clo)+(g4*Cno);
    Cpbeta=(g3*Clbeta)+(g4*Cnbeta);
    Cpp=(g3*Clp)+(g4*Cnp);
    Cpr=(g3*Clr)+(g4*Cnr);
    Cpdela=(g3*Cldela)+(g4*Cndela);
    Cpdelr=(g3*Cldelr)+(g4*Cndelr);
    Cro=(g4*Clo)+(g8*Cno);
    Crbeta=(g4*Clbeta)+(g8*Cnbeta);
    Crp=(g4*Clp)+(g8*Cnp);
    Crr=(g4*Clr)+(g8*Cnr);
    Crdela=(g4*Cldela)+(g8*Cndela);
    Crdelr=(g4*Cldelr)+(g8*Cndelr);
    
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    stheta=(sin(y)*cos(alpha)*cos(beta))+(cos(y)*cos(mu)*sin(alpha))+ (cos(y)*sin(mu)*sin(beta)*cos(alpha));
    theta=asin(stheta);
    cphi=((-sin(y)*sin(alpha)*cos(beta))+(cos(y)*cos(mu)*cos(alpha))-(cos(y)*sin(mu)*sin(beta)*sin(alpha)))/cos(theta);
    sphi=((-sin(y)*sin(beta))+(cos(y)*sin(mu)*cos(beta)))/cos(theta);
    phi=atan(sphi/cphi);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E11=((Va*cos(y)/R)^2)*[-g1*cos(theta)*sin(theta)*sin(phi)-g2*cos(phi)*sin(phi)*(cos(theta))^2;
        -g5*cos(theta)*sin(theta)*cos(phi)-g6*((sin(theta))^2-(cos(theta)*cos(phi))^2);
        -g7*cos(theta)*sin(theta)*sin(phi)-g1*cos(phi)*sin(phi)*(cos(theta))^2];
    
    E12=0.5*rho*S*Va^2*[g3,0,g4;0,g9,0;g4,0,g8];
    
    E13=[-KTp*(Komg*del_t)^2;0;0];
    
    E14=[b,0,0;0,c,0;0,0,b];
    
    E15=[Cpo+Cpbeta*beta+((Cpp*sin(theta)+Cpr*cos(theta)*cos(phi))*b*cos(y)*0.5/R)+Cpdela*del_a+Cpdelr*del_r;
        Cmo+Cmalpha*alpha+(Cmq*c*cos(theta)*sin(phi)*cos(y)*0.5/R)+Cmdele*del_e;
        Cro+Crbeta*beta+((-Crp*sin(theta)+Crr*cos(theta)*cos(phi))*b*cos(y)*0.5/R)+Crdela*del_a+Crdelr*del_r];
    
    E1=E11+E12*(E13+(E14*E15));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    E21=-(((Va*cos(y))^2)/R)*[-cos(mu)*sin(beta)*cos(alpha)+sin(mu)*sin(alpha);
        cos(mu)*cos(beta);
        -cos(mu)*sin(beta)*sin(alpha)-sin(mu)*cos(alpha)];
    
    E22=9.8*[-sin(theta);cos(theta)*sin(phi);cos(theta)*cos(phi)];
    
    E23=(rho*Sprop*Cprop*0.5/m)*[(Kmotor*del_t)^2-(Va*cos(alpha)*cos(beta))^2;0;0];
    
    E24=((rho*S*Va^2)/(2*m))*eye(3);%*[-cos(alpha),0,sin(alpha);0,1,0;-sin(alpha),0,-cos(alpha)];
    
    %CL=CLo+(CLalpha*alpha);
    
    p=-sin(theta)*Va*cos(y)/R;
    q=sin(phi)*cos(theta)*Va*cos(y)/R;
    r=cos(phi)*cos(theta)*Va*cos(y)/R;
    
    E25=[CX+(CXq*q*c*0.5/Va)+(CXdele*del_e);
        CYo+(CYbeta*beta)+((CYp*p+CYr*r)*(b*0.5/Va))+(CYdela*del_a)+(CYdelr*del_r);
        CZ+(CZ*alpha)+(CZq*q*c*0.5/Va)+(CZdele*del_e)];

    E2=E21+E22+E23+(E24*E25);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    E3=(CYbeta*beta)+((CYp*p+CYr*r)*(b*0.5/Va))+(CYdela*del_a)+(CYdelr*del_r);
    
    residue=[E1;E2;E3];
    
end




