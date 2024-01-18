function Struct= Model_Parameters()
    % Inertial Properties - Aerosonde UAV
    Struct.mass=13.5;
    Struct.Jx=0.8244;
    Struct.Jy=1.135;
    Struct.Jz=1.759;
    Struct.Jxz=0.1204;

    % Input Parameters
    Struct.g=9.81;
    Struct.S=0.55;
    Struct.b=2.8956;
    Struct.c=0.18994;
    Struct.Sp=0.2027;
    Struct.rho=1.2682;
    Struct.k_mo=80;
    Struct.Ktp=0;
    Struct.Kw=0;
    Struct.e=0.9;

    Struct.CL0=0.28;
    Struct.CD0=0.03;
    Struct.Cm0=-0.02338;
    Struct.Cla=3.45;
    Struct.Cda=0.30;
    Struct.Cma=-0.38;
    Struct.Clq=0;
    Struct.Cdq=0;
    Struct.Cmq=-3.6;
    Struct.Clde=-0.36;
    Struct.Cdde=0;
    Struct.Cmde=-0.25;
    Struct.Cp=1;
    Struct.Cdp=0.0437;
    Struct.Cndr=-0.032;

    Struct.Cy0=0;
    Struct.Cl0=0;
    Struct.Cn0=0;
    Struct.Cyb=-0.98;
    Struct.Clb=-0.12;
    Struct.Cnb=0.25;
    Struct.Cyp=0;
    Struct.Clp=-0.26;
    Struct.Cnp=0.022;
    Struct.Cyr=0;
    Struct.Clr=0.14;
    Struct.Cnr=-0.35;
    Struct.Cyda=0;
    Struct.Clda=0.08;
    Struct.Cnda=0.06;
    Struct.Cydr=-0.17;
    Struct.Cldr=0.105;

    % Mass Moment Co-efficients
    Struct.G=Struct.Jx*Struct.Jz - Struct.Jxz*Struct.Jxz;
    Struct.G1=(1/Struct.G)*Struct.Jxz*(Struct.Jx-Struct.Jy+Struct.Jz);
    Struct.G2=(1/Struct.G)*(Struct.Jz*(Struct.Jz-Struct.Jy)+Struct.Jxz*Struct.Jxz);
    Struct.G3=Struct.Jz/Struct.G;
    Struct.G4=Struct.Jxz/Struct.G;
    Struct.G5=(1/Struct.Jy)*(Struct.Jz-Struct.Jx);
    Struct.G6=Struct.Jxz/Struct.Jy;
    Struct.G7=(1/Struct.G)*((Struct.Jx-Struct.Jy)*Struct.Jx + Struct.Jxz*Struct.Jxz);
    Struct.G8=Struct.Jx/Struct.G;

end

