function F = HBV_infect_eqs_acute(t,y)
%This is the function file contains the equations for the HBV_infect_model
%code that runs the infection kinetics

%Acute Infection --> Vex accumulates
%Rate constants ===========================================================
k_a = 0.003; %0.288; %1/min; %0.0048 %surface attachment rate (1/s) ????
k_e = 0.006; %1/min; %6.5E-4; %endocytosis rate (1/s) 
nk_fus = 0.02; %1/min; %0.0023; %endosome fusion rate (1/s)
k_tran = nk_fus; %1/min; %0.0023; %transfer to lysosome (1/s)
k_nuc = 1/11; %1/min; %1/(60+600); %virus nuclear import (1/s)
%k_alpha = k_-nuc; %1/min; %nuclear import of RNA+ (1/s)
k_mRNAtranscript = 0.333; %0.333; %1/min; 1/180; %mRNAsub transcription rate (1/s) ****
k_pmRNAtranscript = 0.333; %pmRNAsub transcription rate *****
k_revtranscript = 0.333; %reverse transcription rate (1/s)
k_translate = 0.143; %1/min; 1/420; %translation rate (1/s)
k_RNAtransport = 1/11; %1 %RNA nuclear export 10 min, 1 min to reach ribosomes
k_beta = 1/25; %1/min; 1/90 %pmRNA to RNA+ ?????
k_viralassemble = 0.1; %1/min; 1/600; %viral assembly
k_secrete = 1; %1/min; 1/60; %secretion rate
k_out = 0.216; %1/min; 0.0036; %viral export rate (1/s)
k_degrade = 0.001; %1000 min; 1E-2; %degradation rate (degradation too slow)
k_repair = 1/10; %repair genetic material
k_pack = 1/10; %packaging rate

%Equations ================================================================

Vs = y(1);
Vene = y(2);
Vint = y(3);
Gen = y(4);
cccDNA = y(5);
mRNAsub = y(6);
pmRNAsub = y(7);
pmRNAcyt = y(8);
POLcore = y(9);
RNAplus = y(10);
DNAminus = y(11);
DNAplus = y(12);
mRNAcyt = y(13);
Ant = y(14);
PreCP = y(15);
Vnew = y(16);
Vex = y(17);


F(1,1) = k_a*Vex - k_e*Vs; %Surface virus
F(2,1) = k_e*Vs - nk_fus*Vene - k_tran*Vene; %Endosomal virus
F(3,1) = nk_fus*Vene - k_nuc*Vint; % Internalized virus

F(4,1) = k_nuc*(Vint + Vnew) - k_repair*Gen - k_degrade*Gen; %Raw viral genetic material
F(5,1) = k_repair*Gen + k_RNAtransport*RNAplus - k_degrade*cccDNA; %cccDNA DEGRADATION TOO SLOW
F(6,1) = k_mRNAtranscript*cccDNA - k_RNAtransport*mRNAsub - k_degrade*mRNAsub; %mRNAsubgenomic
F(7,1) = k_pmRNAtranscript*cccDNA - k_RNAtransport*pmRNAsub - k_degrade*pmRNAsub; %pmRNAsubgenomic
F(8,1) = k_RNAtransport*pmRNAsub - k_pack*pmRNAcyt - k_degrade*pmRNAcyt; %pmRNAcyt
F(9,1) = k_translate*pmRNAcyt - k_pack*POLcore - k_degrade*POLcore; %POL and core particle from pmRNA

F(10,1) = k_beta*pmRNAcyt - k_RNAtransport*RNAplus - k_revtranscript*RNAplus - k_degrade*RNAplus; %RNA+
F(11,1) = k_revtranscript*RNAplus - k_repair*DNAminus - k_degrade*DNAminus; %DNA-
F(12,1) = k_repair*DNAminus - k_viralassemble*DNAplus*Ant - k_degrade*DNAplus; %DNA+
F(13,1) = k_RNAtransport*mRNAsub - k_degrade*mRNAcyt; %mRNAcyt
F(14,1) = k_translate*mRNAcyt - k_viralassemble*Ant - k_secrete*Ant - k_degrade*Ant; %surface antigen 
F(15,1) = k_translate*mRNAcyt - k_secrete*PreCP - k_degrade*PreCP; %Pre core protein
F(16,1) = k_viralassemble*DNAplus*Ant - k_out*Vnew; %new virus %PRODUCTION RATE IS TOO FAST
F(17,1) = -k_a*Vex + k_out*Vnew; %Extracellular virus

end

