function F=pathway_model_many_antibody_immune_IgG_trans(t,y,para,para_new)
for i = 1: 603
    y(i) = max(0,y(i));
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j,1) = -para(i)*y(10*(i-1)+j)*y(201)+ para_new(j)*y(10*(i-1)+j+201)-para(11)*y(10*(i-1)+j) + para(12)*y(10*(i-1)+j+201)+para(22)*y(10*(i-1)+j+100)*y(10*(i-1)+j+301)/(1e9+y(10*(i-1)+j+301))-para(16)*y(10*(i-1)+j)*y(402) + para(17)*y(10*(i-1)+j+403)+para(15)*y(10*(i-1)+j+403);
        F_new(10*(i-1)+j) = -para(i)*y(10*(i-1)+j)*y(201)+ para_new(j)*y(10*(i-1)+j+201);
        F_new_self(10*(i-1)+j) = -para(16)*y(10*(i-1)+j)*y(402) + para(17)*y(10*(i-1)+j+403);
    end
end


for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+100,1) = -para(i)*y(10*(i-1)+j+100)*y(201) + para_new(j)*y(10*(i-1)+j+301) -para(21)*y(10*(i-1)+j+100) + para(12)*y(10*(i-1)+j+301)-para(22)*y(10*(i-1)+j+100)*y(10*(i-1)+j+301)/(1e9+y(10*(i-1)+j+301))-para(19)*y(10*(i-1)+j+100)*y(403) + para(17)*y(10*(i-1)+j+503)+para(15)*y(10*(i-1)+j+503);
        F_new_IgM(10*(i-1)+j+100) = -para(i)*y(10*(i-1)+j+100)*y(201) + para_new(j)*y(10*(i-1)+j+301);
        F_new_self_IgM(10*(i-1)+j+100) = -para(19)*y(10*(i-1)+j+100)*y(403) + para(17)*y(10*(i-1)+j+503);
    end
end

F(201,1) = para(13)*y(201)+sum(F_new)+sum(F_new_IgM);


for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+201,1) = para(i)*y(10*(i-1)+j)*y(201) - para_new(j)*y(10*(i-1)+j+201) - para(14)*y(10*(i-1)+j+201);
        F(10*(i-1)+j+301,1) = para(i)*y(10*(i-1)+j+100)*y(201) - para_new(j)*y(10*(i-1)+j+301) - para(14)*y(10*(i-1)+j+301);
    end
end

F(402,1) = para(18)+ sum(F_new_self); %% with environmental antigen
F(403,1) = para(20)+ sum(F_new_self_IgM); %% with environmental antigen

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+403,1) = para(16)*y(10*(i-1)+j)*y(402) - para(17)*y(10*(i-1)+j+403) - para(14)*y(10*(i-1)+j+403);
        F(10*(i-1)+j+503,1) = para(19)*y(10*(i-1)+j+100)*y(403) - para(17)*y(10*(i-1)+j+503) - para(14)*y(10*(i-1)+j+503);
    end
end

end

