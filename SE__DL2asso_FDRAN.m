clear all;  %  important %nakagami-m   coverage probability analitical expression
            %jlf            

% % ============ SPATIAL MODEL ========= %nagani-m channeel

lambda_d =50; %(NUMBER OF RECEIVING VEHICULAR NODES PER KM)


lambda_v = 20; %NUMBER OF TIER 1 NODES PER KM SQ. 




sigam2dbm=-96; %这里不太对，应该是-96左右。但是变化太小了实在
sigma2 = 10^(sigam2dbm/10-3)

% % ============ TRANSMIT PARAMETERS ========= %
pd_dbm =23; % Transmit power of tier 1 node in dBm
pd = 10^(.1*(pd_dbm-30));
pv_dbm = 20; % Transmit power of tier 2 node in dBm
pv = 10^(.1*(pv_dbm-30));






nn=1
for i=2.1:0.2:6
% % ============ PROPAGATION PARAMETERS =================================================== %
alphad =i; % PATH_LOSS EXPONENT 
alphau = i;  %tier 2




% ~~~ SHADOWING PARAMETERS ~~~~ %
ln_mu1 = 0;% Mean of log-normal shadowing gain in dB for TIER 1
ln_sig1 = 2; % Std deviation of log-normal shadowing gain in dB for TIER 1



% % % ======================================================================================


sir_threshold_dB =10;%-50:5:50;
sir_threshold = 10.^(.1*sir_threshold_dB);

%---------------END of PARAMETERS ------------------------------------------------------




sf_c = exp( 1/alphad*log(10)/10*ln_mu1 + (1/alphad*log(10)/10)^2*ln_sig1^2/2 );   %tier 1  1/alpha是一维的


%替代理论后的密度
sh_lambda_1 = sf_c*lambda_d; %1 
%sh_lambda_2 = sf_u0*lambda_2; %20
sh_lambda_a = sf_c*(lambda_v);%21


%%%%%%%输入参数
lambda = sh_lambda_1 ;
alpha = alphad ;
p=pd;




%pdfrm = @(xM) 2*sh_lambda_1*xM.*exp(-sh_lambda_1*xM.^2);

gamma=@(n) factorial( n-1 );

pdf_rn= @(x,n) (2*lambda*x)^n*exp(-2*lambda*x)/(x*gamma(n));
pdf_r_1_to_n=@(x,n) (2*lambda)^n*exp(-2*lambda*x);

 
% j=@(x, alpha,p) sir_threshold*x^alpha;%这里把p去掉试试

%单个

lapI_1 =  @(j, x_d,pa)   exp(-2*lambda*integral (@(x)  (1-1/(j*pa*x^(-alpha) +1))  ,x_d, Inf,'ArrayValued', true)       );

% %2个拉普拉斯
% lapI_12 =  @(j1,j2, x_1,x_2)   exp(-2*lambda*integral (@(x)  (1-1/( (j1*x^(-alpha) +1) *(j2*x^(-alpha) +1)         ) )  ,x_2, Inf,'ArrayValued', true)       );
% 
% cov_dl_1=@( threshold,pa)   integral(@(x)  pdf_rn(x,1) *exp(-threshold/pa*sigma2*x^alpha)*lapI_1 ( threshold*x^alpha, x ),0, Inf,'ArrayValued', true   );
% 
% 
% 
% cov_dl_2=@( threshold,pa)   integral(@(x)   pdf_rn(x,2) *exp(-threshold/pa*sigma2*x^alpha)*lapI_1 ( threshold*x^alpha, x ) ,0, Inf,'ArrayValued', true   );

cov_dl_12=@( pa) integral(@(x_2)  ...
                                                    integral( @(x_1 )   pdf_r_1_to_n(x_2,2)*...
                                                                         integral( @(threshold )   exp(-( exp(threshold)-1 )*sigma2/(pa*(  x_1^(-alpha)+x_2^(-alpha)     )       ))*...
                                                    lapI_1( ( exp(threshold)-1 )/(pa*(x_1^(-alpha)+x_2^(-alpha))),x_2,pa ), 0, inf,'ArrayValued', true ), 0, x_2,'ArrayValued', true )...
                                                      ,0, Inf,'ArrayValued', true );

pa=p;
% aa=cov_dl_1(sir_threshold,pa)
% bb=cov_dl_2(sir_threshold,pa)

ab=cov_dl_12(pa)

se2DL(nn) =ab;
nn=nn+1;
end
plot(se2DL)

