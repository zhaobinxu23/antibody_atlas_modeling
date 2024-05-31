function F=pathway_model_single_antibody_immune_res(t,y,para)

for i = 1: 5
    y(i) = max(0,y(i));
end

F(1,1) = para(10)*y(3)-para(1)*y(1)*y(2)+para(3)*y(3)+para(6)*y(5)+para(8)*y(5)-para(7)*y(1)*y(4)-para(2)*y(1);
F(2,1) = para(10)*y(3)-para(1)*y(1)*y(2)+para(4)*y(2);
F(3,1) = -para(5)*y(3)+para(1)*y(1)*y(2)-para(10)*y(3);
F(4,1) = para(9)-para(7)*y(1)*y(4)+para(8)*y(5);
F(5,1) = para(7)*y(1)*y(4)-para(8)*y(5)-para(5)*y(5);

end




%created by the program testexcel_IL
