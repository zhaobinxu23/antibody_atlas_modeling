function F=pathway_model_many_antibody_immune_res_new(t,y,para,para_new)
for i = 1: 302
    y(i) = max(0,y(i));
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j,1) = -para(i)*y(10*(i-1)+j)*y(101)+ para_new(j)*y(10*(i-1)+j+101)-para(11)*y(10*(i-1)+j) + para(12)*y(10*(i-1)+j+101)-para(16)*y(10*(i-1)+j)*y(202) + para(17)*y(10*(i-1)+j+202)+para(15)*y(10*(i-1)+j+202);

        F_new(10*(i-1)+j) = -para(i)*y(10*(i-1)+j)*y(101)+ para_new(j)*y(10*(i-1)+j+101);
        F_new_self(10*(i-1)+j) = -para(16)*y(10*(i-1)+j)*y(202) + para(17)*y(10*(i-1)+j+202);
    end
end
F(101,1) = para(13)*y(101)+sum(F_new);

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+101,1) = para(i)*y(10*(i-1)+j)*y(101) - para_new(j)*y(10*(i-1)+j+101) - para(14)*y(10*(i-1)+j+101);
    end
end


% F(202,1) = 0; %% without environmental antigen

 F(202,1) = para(18)+ sum(F_new_self); %% with environmental antigen

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+202,1) = para(16)*y(10*(i-1)+j)*y(202) - para(17)*y(10*(i-1)+j+202) - para(14)*y(10*(i-1)+j+202);
    end
end

% y(1) = max(0,y(1));
% y(2) = max(0,y(2));
% y(3) = max(0,y(3));
% y(4) = max(0,y(4));
% y(5) = max(0,y(5));
% 
% F(1,1)=-para(1)*y(1)*y(3)+para(6)*y(4)-para(4)*y(1)+para(5)*y(4);
% 
% F(2,1)=-para(2)*y(2)*y(3)+para(8)*y(5)-para(4)*y(2)+para(5)*y(5);
% 
% F(3,1)=-para(1)*y(1)*y(3)+para(6)*y(4)-para(2)*y(2)*y(3)+para(8)*y(5)+para(3)*y(3);
% 
% F(4,1) = para(1)*y(1)*y(3) - para(6)*y(4) -para(7)*y(4);
% 
% F(5,1) = para(2)*y(2)*y(3)- para(8)*y(5)-para(7)*y(5);

end




%created by the program testexcel_IL
