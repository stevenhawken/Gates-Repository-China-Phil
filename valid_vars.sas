

*******************************************************************;
*data step statements for variables required by GA validation macro;
*******************************************************************;
 
%macro valid_vars;

sga10=.;sga3=.;iugr=.;
bw=bweight_raw;gest=floor(gestage);
if BW ne . and gest ne . and sex ne "" then do; 
if (sex="M" and GEST=24 and BW<498) then sga10=1;
  else if (sex="M" and GEST=25 and BW<569) then sga10=1;
  else if (sex="M" and GEST=26 and BW<649) then sga10=1;
  else if (sex="M" and GEST=27 and BW<738) then sga10=1;
  else if (sex="M" and GEST=28 and BW<837) then sga10=1;
  else if (sex="M" and GEST=29 and BW<948) then sga10=1;
  else if (sex="M" and GEST=30 and BW<1071) then sga10=1;
  else if (sex="M" and GEST=31 and BW<1208) then sga10=1;
  else if (sex="M" and GEST=32 and BW<1359) then sga10=1;
  else if (sex="M" and GEST=33 and BW<1438) then sga10=1;
  else if (sex="M" and GEST=34 and BW<1709) then sga10=1;
  else if (sex="M" and GEST=35 and BW<1957) then sga10=1;
  else if (sex="M" and GEST=36 and BW<2182) then sga10=1;
  else if (sex="M" and GEST=37 and BW<2386) then sga10=1;
  else if (sex="M" and GEST=38 and BW<2570) then sga10=1;
  else if (sex="M" and GEST=39 and BW<2735) then sga10=1;
  else if (sex="M" and GEST=40 and BW<2882) then sga10=1;
  else if (sex="M" and GEST=41 and BW<3011) then sga10=1;
  else if (sex="M" and GEST=42 and BW<3124) then sga10=1;
  else if (sex="F" and GEST =24 and BW<470) then sga10=1;
  else if (sex="F" and GEST =25 and BW<537) then sga10=1;
  else if (sex="F" and GEST =26 and BW<613) then sga10=1;
  else if (sex="F" and GEST =27 and BW<697) then sga10=1;
  else if (sex="F" and GEST =28 and BW<791) then sga10=1;
  else if (sex="F" and GEST =29 and BW<896) then sga10=1;
  else if (sex="F" and GEST =30 and BW<1012) then sga10=1;
  else if (sex="F" and GEST =31 and BW<1141) then sga10=1;
  else if (sex="F" and GEST =32 and BW<1284) then sga10=1;
  else if (sex="F" and GEST =33 and BW<1409) then sga10=1;
  else if (sex="F" and GEST =34 and BW<1679) then sga10=1;
  else if (sex="F" and GEST =35 and BW<1923) then sga10=1;
  else if (sex="F" and GEST =36 and BW<2141) then sga10=1;
  else if (sex="F" and GEST =37 and BW<2334) then sga10=1;
  else if (sex="F" and GEST =38 and BW<2504) then sga10=1;
  else if (sex="F" and GEST =39 and BW<2652) then sga10=1;
  else if (sex="F" and GEST =40 and BW<2780) then sga10=1;
  else if (sex="F" and GEST =41 and BW<2887) then sga10=1;
  else if (sex="F" and GEST =42 and BW<2976) then sga10=1;
  else sga10=0;

if (sex="M" and GEST=24 and BW<443  ) then sga3=1;
  else if (sex="M" and GEST=25 and BW<507) then sga3=1;
  else if (sex="M" and GEST=26 and BW<578) then sga3=1;
  else if (sex="M" and GEST=27 and BW<657) then sga3=1;
  else if (sex="M" and GEST=28 and BW<746) then sga3=1;
  else if (sex="M" and GEST=29 and BW<845) then sga3=1;
  else if (sex="M" and GEST=30 and BW<954) then sga3=1;
  else if (sex="M" and GEST=31 and BW<1076) then sga3=1;
  else if (sex="M" and GEST=32 and BW<1184) then sga3=1;
  else if (sex="M" and GEST=33 and BW<1211) then sga3=1;
  else if (sex="M" and GEST=34 and BW<1456) then sga3=1;
  else if (sex="M" and GEST=35 and BW<1705) then sga3=1;
  else if (sex="M" and GEST=36 and BW<1931) then sga3=1;
  else if (sex="M" and GEST=37 and BW<2136) then sga3=1;
  else if (sex="M" and GEST=38 and BW<2321) then sga3=1;
  else if (sex="M" and GEST=39 and BW<2487) then sga3=1;
  else if (sex="M" and GEST=40 and BW<2634) then sga3=1;
  else if (sex="M" and GEST=41 and BW<2764) then sga3=1;
  else if (sex="M" and GEST=42 and BW<2878) then sga3=1;
  else if (sex="F" and GEST =24 and BW<419) then sga3=1;
  else if (sex="F" and GEST =25 and BW<479) then sga3=1;
  else if (sex="F" and GEST =26 and BW<546) then sga3=1;
  else if (sex="F" and GEST =27 and BW<621) then sga3=1;
  else if (sex="F" and GEST =28 and BW<704) then sga3=1;
  else if (sex="F" and GEST =29 and BW<798) then sga3=1;
  else if (sex="F" and GEST =30 and BW<901) then sga3=1;
  else if (sex="F" and GEST =31 and BW<1016) then sga3=1;
  else if (sex="F" and GEST =32 and BW<1144) then sga3=1;
  else if (sex="F" and GEST =33 and BW<1198) then sga3=1;
  else if (sex="F" and GEST =34 and BW<1464) then sga3=1;
  else if (sex="F" and GEST =35 and BW<1705) then sga3=1;
  else if (sex="F" and GEST =36 and BW<1919) then sga3=1;
  else if (sex="F" and GEST =37 and BW<2110) then sga3=1;
  else if (sex="F" and GEST =38 and BW<2278) then sga3=1;
  else if (sex="F" and GEST =39 and BW<2423) then sga3=1;
  else if (sex="F" and GEST =40 and BW<2548) then sga3=1;
  else if (sex="F" and GEST =41 and BW<2654) then sga3=1;
  else if (sex="F" and GEST =42 and BW<2740) then sga3=1;
 else sga3=0;
end;
IUGR=sga3; 
if BW>=3000 then bw3k=0;
else if 0<BW<3000 then bw3k=1;

if gestage>=37 then  preterm=1;
else if 0<gestage<37 then preterm=2;

if 0<gestage<28 then premie=1;
else if 28=<gestage<32 then premie=2;
else if 32=<gestage<37 then premie=3;
else if gestage>=37 then premie=4;
gest_wk=floor(gestage);
if BW>=2500 then lbw=1;
else if 0<BW<2500 then lbw=2;
format sex multiple;
%mend valid_vars;