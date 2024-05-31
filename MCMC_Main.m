%% main
clc
clear
test_speed_one_run;

exp_data(:,1) = timepoints;
exp_data(:,2) = person_ll;
exp_data(:,3) = person_ll_virus;

params0 = [-15.5;0.5;4;1;0.4];
ll_init = ll_function(exp_data,params0);
paraset = [];
paraset = [paraset params0];
loglikelihood_set = [];
loglikelihood_set = [loglikelihood_set ll_init];
all_guesses = [];
all_guesses = [all_guesses params0];

run_number = 10000;
for i = 2:run_number
%     if rem(i,100) == 0
%         i/run_number
%     end
    paramtest = [random('Uniform',-15.5,-15.5);random('Uniform',0.5,0.5);random('Uniform',3,7);random('Uniform',1,1);random('Uniform',0.2,0.8)];
    ll_test = ll_function(exp_data,paramtest);
    % Metroplis-Hastings
    if ll_test >= loglikelihood_set(end) + 1*log(random('Uniform',0,1))
        loglikelihood_set = [loglikelihood_set ll_test];
        paraset = [paraset paramtest];
        all_guesses = [all_guesses paramtest];
    else
        loglikelihood_set = [loglikelihood_set loglikelihood_set(end)];
        paraset = [paraset paraset(:,end)];
        all_guesses = [all_guesses paramtest];
    end
end

    

