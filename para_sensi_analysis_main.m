%% parameter sensitivity analysis  7 parameters  including ddddd
clc
clear
theta_ori = [-15.5 0.5 1.5 0.8 5 1 0.5];

for kk = 1:4
    kk
    theta = theta_ori;
    theta(kk) = theta_ori(kk)+0.1;
    

    
    [integral_value,max_virus,max_symptom,max_elisa,final_elisa] = datagen_sensitivity_analysis(theta);
    
    theta(kk) = theta_ori(kk);

   
    [integral_value_new,max_virus_new,max_symptom_new,max_elisa_new,final_elisa_new] = datagen_sensitivity_analysis(theta);
    
    final_sen_results(kk,1) = ((integral_value_new - integral_value)/integral_value);% /(0.2/0.9);
    final_sen_results(kk,2) = ((max_virus_new - max_virus)/max_virus);% /(0.2/0.9);
    final_sen_results(kk,3) = ((max_symptom_new - max_symptom)/max_symptom);% /(0.2/0.9);
    final_sen_results(kk,4) = ((max_elisa_new - max_elisa)/max_elisa); % /(0.2/0.9);
    final_sen_results(kk,5) = ((final_elisa_new - final_elisa)/final_elisa);% /(0.2/0.9);
end

for kk = 5:7
    kk
    theta = theta_ori;
    theta(kk) = theta_ori(kk)*0.9;
    

    
    [integral_value,max_virus,max_symptom,max_elisa,final_elisa] = datagen_sensitivity_analysis(theta);
    
    theta(kk) = theta_ori(kk)*1.1;

   
    [integral_value_new,max_virus_new,max_symptom_new,max_elisa_new,final_elisa_new] = datagen_sensitivity_analysis(theta);
    
    final_sen_results(kk,1) = ((integral_value_new - integral_value)/integral_value)/(0.2/0.9);
    final_sen_results(kk,2) = ((max_virus_new - max_virus)/max_virus) /(0.2/0.9);
    final_sen_results(kk,3) = ((max_symptom_new - max_symptom)/max_symptom)/(0.2/0.9);
    final_sen_results(kk,4) = ((max_elisa_new - max_elisa)/max_elisa)/(0.2/0.9);
    final_sen_results(kk,5) = ((final_elisa_new - final_elisa)/final_elisa)/(0.2/0.9);
end
    
    

