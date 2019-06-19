R_ref= 0.0422 ;
R_atm= R_ref *(3/1000+1) ; %lets start with atmospheric d34S of 3 permil, and later do a spin up.
COS32=500; % This is the amount in the atmosphere. for simplicity I work here with concetrations in ppt instead of mass
COS34=COS32*R_atm ; % This is the amount of COS with 34S
epsilon_U=-5; %assuming fractionation of -5 permil in uptake by plants
%R_P=R_ref*(20/1000+1); % Assuming 20 permil for production from the ocean
%R_P=R_ref*(3/1000+1); % Assuming 3 permil for anthropogenic production
R_P=R_ref*(8/1000+1); % Assuming 8 permil for total production from both sources
P_COS32=[50,110,130,140,150,200,250,300,250,200,150,50]' ; %assuming arbitrary production for each mont
U_COS32=[25,50,60,70,100,300,500,400,200,150,100,25]' ;  %assuming arbitrary uptake for each mont
Dt=0.1 ; % time step

P_COS34=P_COS32 * R_P ;
for y=1:10 % years loop
    for m=1:12 % months loop
        R_atm=COS34/COS32 ;
        U_COS34 = U_COS32(m) *(1+epsilon_U /1000)*R_atm ;
        COS32=COS32 + P_COS32(m)*Dt- U_COS32(m)*Dt ; % COS 32 iteration
        COS34=COS34 + P_COS34(m)*Dt- U_COS34*Dt ; % COS 34 iteration
        d34S(y,m)=((COS34/COS32)/R_ref-1)*1000 ; %calculate delta-34S
        COS(y,m)=COS32 ; % save the COS concetration for later 
    end
end
subplot(2,1,1) 
plot(1:12,COS(10,:))
subplot(2,1,2)
plot(1:12,d34S(10,:))



