total_dil <- function(Cf,Mf = 2.8E4,Uw = 1E5,Pshed = 1, Ptotal) { #Concentration of target - gc/ mg-dw, Mass rate of feces production (wet) - 28,000 mg-dry/day, solid fraction of feces - 0.2, daily water use - 100,000 mL/day.
  (Cf*Mf*Pshed)/(Uw*Ptotal)
}

l_dil <- function(K,TSS = 0.28,Pshedders, Ptotal,ffecal,Cfeces){ #0.28 mg/mL - calculated from Rose + 100L/day water use
  
  (TSS*ffecal*Cfeces)/(Ptotal*(1+K*TSS))
}

s_dil <- function(K,TSS = 0.28,Pshedders, Ptotal,ffecal,Cfeces){
  
  (K*(TSS*ffecal*(Pshedders/Ptotal)*Cfeces))/(1+K*TSS)
}