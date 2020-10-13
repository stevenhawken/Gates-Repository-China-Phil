

%let today = %sysfunc(today(), date9.);/*macro to add todays date to filenames*/

/*macro to derive Harrell's percentile knot locations for restricted cubic splines*/

%macro spline_pctls(data,var,suff,pctls);
proc datasets noprint;
delete percentiles;
run;

proc univariate data=&data noprint;
   var &var;
   output out=percentiles pctlpts=&pctls pctlpre=p_; /* specify the percentiles */
run;
%global pctls_&var.&suff; 
data _null_;
set percentiles;
call symput("pctls_&var.&suff",catx(',', of _numeric_));   /* put all values into a comma-separated list */
run;
%mend spline_pctls;