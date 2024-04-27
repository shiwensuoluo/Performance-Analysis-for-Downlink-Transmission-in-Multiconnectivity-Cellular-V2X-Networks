clear;
%%%use monte carlo to campute the coverage and rate
%参数
d = 50;  % OBSERVATION WINDOW RADIUS
m_asso =2;
n=1;
for ii = 10
% % ============ SPATIAL MODEL ========= %nagani-m channeel

lambda_d = 50; %NUMBER OF TIER 1 NODES PER KM SQ. 原来是0.5

lambda_v = 20; %(NUMBER OF RECEIVING VEHICULAR NODES PER KM)





% % ============ TRANSMIT PARAMETERS ========= %
p1_dbm = 23; % Transmit power of tier 1 node in dBm
p1 = 10^(.1*(p1_dbm-30));


pv_dbm = 20;
pv = 10^(.1*(pv_dbm-30));

sigam2dbm=-96 ;
sigma2 = 10^(sigam2dbm/10-3);

nn=1
for al=2.1:0.3:6
    
  
alphad =al; % PATH_LOSS EXPONENT 
alpha = al;  %tier 2



% ~~~ SHADOWING PARAMETERS ~~~~ %
ln_mu1 = 0;% Mean of log-normal shadowing gain in dB for TIER 1
ln_sig1 = 2; % Std deviation of log-normal shadowing gain in dB for TIER 1



% % % ======================================================================================

num_iterations = 10;%原文是1e4

sir_threshold_dB =ii;
sir_threshold = 10^(.1*sir_threshold_dB);


rate_threshold = 10; % in Mbps
%-------------------------------------------------------------------
%======================Simulation=============================

sir_measured = zeros(1, num_iterations);



%%仿真布置
numdbs = poissrnd(lambda_d*(d),1,num_iterations); % Number of TIER 1 nodes in observation window for all iterations   服从lambade分布的  

for j1=1:num_iterations

    num_DBS = numdbs(j1); % Number of TIER 1 nodes
    
    %生成下行基站位置
%     loc_DBS=poissrnd(lambda_v,[1,num_DBS])/max(poissrnd(lambda_v,[1,num_DBS]));
%     location_DBS = [-d/2+d*loc_DBS];
    location_DBS = [-d/2+d*rand(1,num_DBS)];
    [~,index_loc_DBS] = sort (abs(location_DBS) );
    
    %生成车辆位置
    num_v = poissrnd(lambda_v*d);
   location_V=[-d/2+d*rand(1,num_v)];
    
%      loc_v=poissrnd(lambda_v,[1,num_v])/max(poissrnd(lambda_v,[1,num_v]));
%      location_V=[-d/2+d*loc_v];
   
    [~,index_loc_veh] = sort (abs(location_V) );
    
   

    % ============ TIER 2 generation ===================
   


    
    % ==================== ============== ============= ======== ========= ========== =
    % ================== Load & SIR Computation =======================
  
   %val是值 ind是索引位置
   


                     
                  
        for iv = 1:num_v
            
            
            shad_d = 10.^(.1*(ln_mu1 + ln_sig1*randn(num_DBS,1)))';   %
            
       %%%%shadow     
         
           
  
      %%%%%%%%%%channel parameter  
          chnl_d =exprnd(1,1,num_DBS);
        
    
            p_max_dbs2v = p1*shad_d.*(abs(location_V(iv)-location_DBS)).^(-alphad); %找距离车辆最近的基站
            [~,max_to_min_DL] = sort( p_max_dbs2v);%%找这个一大列最大值。   下行  %最大功率的SBS和MBS
                         
                         
                         
            loc_m_asso = max_to_min_DL(num_DBS-m_asso+1:num_DBS);%因为是从小到大排，所以调整一下
            
            
            p_sinr_dbs2v=p1*shad_d.*chnl_d.*(abs(location_V(iv)-location_DBS)).^(-alphad);
            
            p_sinr_dl = sum(p_sinr_dbs2v(loc_m_asso) );%加了信道增益的分子
            
            loc_sinr=max_to_min_DL(1:num_DBS-m_asso);
            
            sinr=p_sinr_dl/(sum(p_sinr_dbs2v(loc_sinr))  +sigma2);%计算sinr
                
           if sinr>sir_threshold
               cov_pr_v(iv)=1;
           else
               cov_pr_v(iv)=0;
                 
           end

        end  
        

    
 cov_pr_d(j1)= sum(cov_pr_v) / num_v

    
    
end

sim_cov_2d_alpha(nn) = mean(cov_pr_d)
nn=nn+1;
end
end
plot(sim_cov_2d_alpha)









