
double_X = table2array (X_indicators_GWP);
double_Y = table2array (Y_GWP_100);
[indices_GWP_indicators,d,t,e,w,bm,dm,tm,em,wm]=betaKS3(double_X ,double_Y);

% Write the table to a CSV file
csvwrite('Users/GF20PZ/OneDrive - Aalborg Universitet/Dokumenter/AAU/FIsh Farm model/From July 2022/borgonovo_indices/borgonovo_GWP100_indicators.csv',indices_GWP_indicators)







double_X = table2array (X_FEP_parameters);
double_Y = table2array (Y_FEP);
[indices_FEP_parameters,d,t,e,w,bm,dm,tm,em,wm]=betaKS3(double_X ,double_Y);

% Write the table to a CSV file
csvwrite('Users/GF20PZ/OneDrive - Aalborg Universitet/Dokumenter/AAU/FIsh Farm model/From July 2022/borgonovo_indices/borgonovo_FEP_parameters.csv',indices_FEP_parameters)


double_X = table2array (X_FEP_indicators);
double_Y = table2array (Y_FEP);
[indices_FEP_indicators,d,t,e,w,bm,dm,tm,em,wm]=betaKS3(double_X ,double_Y);

% Write the table to a CSV file
csvwrite('Users/GF20PZ/OneDrive - Aalborg Universitet/Dokumenter/AAU/FIsh Farm model/From July 2022/borgonovo_indices/borgonovo_FEP_indicators.csv',indices_FEP_indicators)



double_X = table2array (X_TETPinf_parameters);
double_Y = table2array (Y_TETPinf);
[indices_TETPinf_parameters,d,t,e,w,bm,dm,tm,em,wm]=betaKS3(double_X ,double_Y);

% Write the table to a CSV file
csvwrite('Users/GF20PZ/OneDrive - Aalborg Universitet/Dokumenter/AAU/FIsh Farm model/From July 2022/borgonovo_indices/borgonovo_TETPinf_parameters.csv',indices_TETPinf_parameters)



double_X = table2array (X_TETPinf_indicators);
double_Y = table2array (Y_TETPinf);
[indices_TETPinf_indicators,d,t,e,w,bm,dm,tm,em,wm]=betaKS3(double_X ,double_Y);

% Write the table to a CSV file
csvwrite('Users/GF20PZ/OneDrive - Aalborg Universitet/Dokumenter/AAU/FIsh Farm model/From July 2022/borgonovo_indices/borgonovo_TETPinf_indicators.csv',indices_TETPinf_indicators)





double_X = table2array (X_eutrophication_parameters);
double_Y = table2array (Y_eutrophication);
[indices_eutrophication_parameters,d,t,e,w,bm,dm,tm,em,wm]=betaKS3(double_X ,double_Y);

% Write the table to a CSV file
csvwrite('Users/GF20PZ/OneDrive - Aalborg Universitet/Dokumenter/AAU/FIsh Farm model/From July 2022/borgonovo_indices/borgonovo_eutrophication_parameters.csv',indices_eutrophication_parameters)



double_X = table2array (X_eutrophication_indicators);
double_Y = table2array (Y_eutrophication);
[indices_eutrophication_indicators,d,t,e,w,bm,dm,tm,em,wm]=betaKS3(double_X ,double_Y);

% Write the table to a CSV file
csvwrite('Users/GF20PZ/OneDrive - Aalborg Universitet/Dokumenter/AAU/FIsh Farm model/From July 2022/borgonovo_indices/borgonovo_eutrophication_indicators.csv',indices_eutrophication_indicators)




FETPinf

double_X = table2array (X_FETPinf_indicators);
double_Y = table2array (Y_FETPinf);
[indices_FETPinf_indicators,d,t,e,w,bm,dm,tm,em,wm]=betaKS3(double_X ,double_Y);

% Write the table to a CSV file
csvwrite('Users/GF20PZ/OneDrive - Aalborg Universitet/Dokumenter/AAU/FIsh Farm model/From July 2022/borgonovo_indices/borgonovo_FETPinf_indicators.csv',indices_FETPinf_indicators)




double_X = table2array (X_FETPinf_parameters);
double_Y = table2array (Y_FETPinf);
[indices_FETPinf_parameters,d,t,e,w,bm,dm,tm,em,wm]=betaKS3(double_X ,double_Y);

% Write the table to a CSV file
csvwrite('Users/GF20PZ/OneDrive - Aalborg Universitet/Dokumenter/AAU/FIsh Farm model/From July 2022/borgonovo_indices/borgonovo_FETPinf_parameters.csv',indices_FETPinf_parameters)

