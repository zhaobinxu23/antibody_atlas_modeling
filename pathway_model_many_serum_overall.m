function F=pathway_model_many_serum_overall(t,y,para)
for i = 1: 3
    y(i) = max(0,y(i));
end

%% y(1) represent antibody  y(2) serum  y(3) = complex
F(1,1) = -para(1)*y(1)*y(2) + para(2)*y(3);
F(2,1) = -para(1)*y(1)*y(2) + para(2)*y(3);
F(3,1) = para(1)*y(1)*y(2) - para(2)*y(3);



end




%created by the program testexcel_IL
