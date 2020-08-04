function HBV_infect_model
% HBV infection kinetics model
% ODEs contained in HBV_infect_eqs_acute.m and HBV_infect_eqs_chronic.m
% Seth Thor
% Seth Jurgens

tspan = [0 1440];% min

%y0 = [0 0 0]; %Vs Vene Vcyt
%y0 = [0 0 0 0]; %Vs Vene Vcyt cccDNA 
%y0 = [0 0 0 0 0]; %Vs Vene Vcyt cccDNA mRNAsub 
%y0 = [0 0 0 0 0 0]; %Vs Vene Vcyt cccDNA mRNAsub pmRNA
%y0 = [0 0 0 0 0 0 0]; %Vs Vene Vcyt cccDNA mRNAsub pmRNA RNA+
%y0 = [0 0 0 0 0 0 0 0]; %Vs Vene Vcyt cccDNA mRNAsub pmRNA RNA+ DNA+
%y0 = [0 0 0 0 0 0 0 0 0]; %Vs Vene Vcyt cccDNA mRNAsub pmRNA RNA+ DNA+ mRNAcyt
%y0 = [0 0 0 0 0 0 0 0 0 0]; %Vs Vene Vcyt cccDNA mRNAsub pmRNA RNA+ DNA+ mRNAcyt Ant
%y0 = [0 0 0 0 0 0 0 0 0 0 0]; %Vs Vene Vcyt cccDNA mRNAsub pmRNA RNA+ DNA+ mRNAcyt Ant CP
y0=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %chronic 
%y0=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100]; %acute
[t,y] = ode45('HBV_infect_eqs_chronic',tspan,y0);
%[t,y] = ode45('HBV_infect_eqs_acute',tspan,y0);

f1 = figure;
plot(t,y(:,1)) %Vs
hold on
plot(t,y(:,2),'--') %Vene
plot(t,y(:,3),'r') %Vint
plot(t,y(:,16),'r--')%Vnew --> accumulating
%plot(t,y(:,17),'g') %Vex
hold off
xlabel('Time (min)'), ylabel('Count')
legend('V_s','V_e_n_e','V_i_n_t','V_n_e_w','Vex')

f2 = figure;
plot(t,y(:,4)) %Raw viral genetic material
hold on
plot(t,y(:,5),'--') %cccDNA --> accumulating slowly, degrade too slow
plot(t,y(:,6),'r') %mRNAsub
plot(t,y(:,13),'r--') %mRNAcyt --> accumulating degrade too slow
plot(t,y(:,7),'g') %pmRNAsub
plot(t,y(:,8),'g--') %pmRNA --> accumulating slowly
hold off
xlabel('Time (min)'), ylabel('Count')
legend('Gen','cccDNA','mRNA_s_u_b','mRNA_c_y_t','pmRNA_s_u_b','pmRNA_c_y_t')

f3 = figure;
plot(t,y(:,9),'k') %POL and core particle from pmRNA
hold on
plot(t,y(:,10)) %RNA+
plot(t,y(:,11),'r') %DNA-
plot(t,y(:,12),'r--') %DNA+
plot(t,y(:,14),'g')%Ant
plot(t,y(:,15),'g--')%Pre-CP --> accumulating
hold off
xlabel('Time (min)'), ylabel('Count')
legend('POLcore','RNA+','DNA-','DNA+','Ant','PreCP')

end

