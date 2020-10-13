/******************************************************************
* Program name: GA_validation_prep_data.sas       

******************************************************************/
proc format;
value prem
1='1. <28 wks   '
2='2. 28-<32 wks'
3='3. 32-<37 wks'
4='4. >=37 wks  '
;
value term
1='1. >=37 wks'
2='2. <37 wks '
;
value prem34wks
0='1. >=34 wks '
1='2. <34 wks'
;
value prem37wks
0='1. >=37 wks '
1='2. <37 wks'
;
value prem34_37wks
1='1. <34 wks '
2='2. >=34 to <37 wks'
3='3. >=37 wks'
;
value lbw
1='1. 2500+'
2='2. <2500'
;
run;

%let today = %sysfunc(today(), date9.);

%macro win(data=,covars=,id=person_id,out=,varlist=,trans=,transvarlist=,tukey=,weekmonth=none,std=y,type=autoscale,where=,plots=y);

proc datasets noprint ;
delete quantiles quantiles2 one two no_zeros;
quit;

data one(keep= &covars &id &weekmonth &varlist) two(drop= &varlist);
set &data;
none=1;
array vars &varlist;
do over vars;
if .<vars<0 then vars=0; /* some values below instrument detection limit  are -ve, set to 0*/
end; 
run;

data no_zeros;
set one;
array vars &varlist;
do over vars;
	if vars=0 then vars=.; /* dataset to determine extreme quantiles for non-zero values for Tukey fence extreme value processing*/
end; 
run;

proc sort data=one;
by &weekmonth;
run;
proc sort data=no_zeros;
by &weekmonth;
run;

proc summary data=no_zeros qmethod=p2;
by &weekmonth;
var &varlist;
output out=quantiles2(drop=_TYPE_ _FREQ_) qrange= q1= q3= 
/ autoname;run;

data one;
merge one quantiles2;
by &weekmonth;
%if %length(&tukey)>0 %then %do;
 %let word_cnt=%sysfunc(countw(&varlist));
%do i=1 %to &word_cnt;
	%let var=%scan(&varlist,&i, %str( ));
	if &var ne . then do;
		&var._low_fence= &var._q1 - ( &tukey * &var._qrange ); 
		&var._up_fence=  &var._q3 + ( &tukey * &var._qrange ); 
		if 0 < &var < &var._low_fence then &var = &var._low_fence;
		else if &var > &var._up_fence then &var = &var._up_fence;
	end;
      drop  &var._q1 &var._q3 &var._qrange &var._low_fence &var._up_fence;  
%end;
%end;
run;

%if %length(&trans)>0 %then %do;
data one;
set one;
 %let word_cnt=%sysfunc(countw(&transvarlist));
%do i=1 %to &word_cnt;
	%let var=%scan(&transvarlist,&i, %str( ));
	if &var ne . then do;
		&var= &trans ( &var );
	end; 
%end;

run;
%end;

%if &std=y %then %do;

	proc summary data=one;
	&where
	by &weekmonth;
	var &varlist;
	output out=quantiles(drop=_TYPE_ _FREQ_) mean= median= stddev= min= max= p10= p90=  / autoname;
	run;

	data one;
	merge one quantiles;
	by &weekmonth;
	run;


	data one;
	set one;
	 %let word_cnt=%sysfunc(countw(&varlist));
	%do i=1 %to &word_cnt;
		%let var=%scan(&varlist,&i, %str( ));

%if &type=autoscale %then %do;
if &var ne . then do; 
	&var = ( &var - &var._mean ) /&var._stddev; 
end;
%end; 

%else %if &type=pareto %then %do; 
if &var ne . then do;
	&var = ( &var - &var._mean ) /sqrt(&var._stddev);
end;
%end; 

%else %if &type=range %then %do; 
if &var ne . then do;
	&var = ( &var - &var._mean ) /(&var._p90 - &var._p10); 
end;
%end; 

%else %if &type=centre %then %do;
if &var ne . then do; 
	&var = ( &var - &var._mean );
end;
%end; 

%else %if &type=centre_med %then %do; 
if &var ne . then do;
	&var = ( &var - &var._median ); 
end;
%end; 

%else %if &type=level %then %do;
if &var ne . then do; 
	&var = ( &var - &var._mean )/ &var._mean ;
end; 
%end;

%else %if &type=level_med %then %do; 
if &var ne . then do;
	&var = ( &var - &var._median )/ &var._median ; 
end;
%end;

drop &var._mean &var._stddev &var._p10 &var._p90 &var._median &var._min &var._max;

%end;

	run;
%end;

%if &plots=y %then %do;
ods _all_ close;
options papersize="ISO A4" orientation=portrait;
ods html;
ods graphics / width=4cm height=3cm;
* start gridded layout with four columns ;
ods layout gridded columns=4 ;
%let word_cnt=%sysfunc(countw(&varlist));
%do i=1 %to &word_cnt;
%let var=%scan(&varlist,&i, %str( ));
ods region;
proc univariate data=one;
var &var ;
histogram/nmidpoints=200;
ods select histogram;
run;
%end;
ods layout end;
ods graphics off;
ods html close;

%end;

proc sort data=one;
by &id ;
run;

proc sort data=two;
by &id;
run;

data &out;
merge two one(drop= &covars);
by &id ;
run;
ods graphics;
ods html;

%mend win;

*Sample Call;

%win
(data=zambia_heel,
covars=gestage,
id=study_id,
out=zambia_heel_pto,
varlist=bweight ASA Ala Arg BIO C0 C2 C3 C4 C5 C6 C8 C10 C12 C14 C16 
C18 C10_1 C12_1 C14OH C14_1 C14_2 C16OH C16_1OH C18OH C18_1 C18_2 
C18_1OH C3DC C4DC C4OH C5DC C5OH C5_1 C6DC C8_1 Cit GALT Gly IRT 
LEU MET Orn PHE TSH TYR Val _17OHP fa_ratio,
trans=log,
transvarlist=/*bweight*/ ASA Ala Arg BIO C0 C2 C3 C4 C5 C6 C8 C10 C12 C14 C16 
C18 C10_1 C12_1 C14OH C14_1 C14_2 C16OH C16_1OH C18OH C18_1 C18_2 
C18_1OH C3DC C4DC C4OH C5DC C5OH C5_1 C6DC C8_1 Cit GALT Gly IRT 
LEU MET Orn PHE TSH TYR Val _17OHP /*fa_ratio*/,
tukey=3,
weekmonth=none,
std=y,
type=pareto,
where=,
plots=y);
