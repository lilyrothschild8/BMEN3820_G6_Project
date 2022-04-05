function dndt = G6_1ODE(t,N,flag,I_max,k1,k2,k3,k4,kp,ke,a,dc)
    %G6_1ODE defines the ODEs for the CD8 T Cell Model
    %call using G6_1ODE(t,N,flag,I_max,k1,k2,k3,k4,p,ke,phi,qs)), 
    %where 'flag' denotes the location where parameters are set
    
        %INPUT parameters:
        %k1 is the infection spread rate
        %I_max is the maximum number of infected cells 
        %k2 is the killing rate of infected cells by CD8 T cells
        %k3 antigen driven CDT8 cell proliferation rate 
        %k4 antigen driven CDT8 cell proliferation rate 
        %kp antigen driven proliferation threshold
        %ke antigen driven suppresion threshold
        %a is the rate of growth of cytokine pathology
        %dc is the rate of loss of cytokine pathology
        %a is the rate growth cytokine pathology 
        %dc is the rate of loss of cytokine pathology 
        
    %OUTPUT:
    %dndt is the population array
    %dndt(1) infected cells = I 
    %dndt(2) activated CD8 T Cells = E
    %dndt(3) pathology = P
    
dndt = zeros(size(N));
dndt(1)=k1.*N(1).*(1-(N(1))/I_max)-k2*N(1).*N(2); 
dndt(2)=k3*(N(1).*N(2))./(kp+N(1))-k4*(N(1).*N(2))/(ke+N(1));
dndt(3)=a*N(1).*N(2)-dc*N(3);
end
