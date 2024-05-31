Main_immune_imprinting_many_antibodies;
clear x0

for i = 1:100
data_new(i,:) = interp1(t,y(:,i),(0:1:500));
end


for k = 1:101
 for i = 1:100
    xx0(i) = data_new(i,k);
 end  
 x0(1) = sum(xx0);
 x0(2) = 0;
 x0(3) = 0;
 para(1) = 10^(-13.5);
 para(2) = 10^(-0.5);
 
[tt, zz]=ode15s(@pathway_model_many_serum_overall,[0 50],x0,[],para);

percent(k) = interp1(tt,zz(:,1),50)/x0(1);
 
end


for k = 1:101
    
   
theta(1) = -15.5;
for i = 1:100
    x0(i) = data_new(i,k)*percent(k);
end

x0(101) = 1e16;%% virus

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+101) = 0;
    end
end

x0(202) = 0;%% environmental antigen concentration

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+202) = 0;
    end
end


para_2(1) = 10^(theta(1)-4.5); 
para_2(2) = 10^(theta(1)-3.5); 
para_2(3) = 10^(theta(1)-2.5); 
para_2(4) = 10^(theta(1)-1.5); 
para_2(5) = 10^(theta(1)-0.5);  
para_2(6) = 10^(theta(1)+0.5); 
para_2(7) = 10^(theta(1)+1.5);  
para_2(8) = 10^(theta(1)+2.5); 
para_2(9) = 10^(theta(1)+3.5);  
para_2(10) = 10^(theta(1)+4.5); 

para_new(1) = 1e-3; 
para_new(2) = 1e-2;
para_new(3) = 1e-1; 
para_new(4) = 1e0;
para_new(5) = 1e1; 
para_new(6) = 1e2;
para_new(7) = 1e3; 
para_new(8) = 1e4;
para_new(9) = 1e5; 
para_new(10) = 1e6;


para_2(11) = 0;%% 0.01
para_2(12) = 0;%% 5
para_2(13) = 0;%% 1
para_2(14) = 0;%% 0.5
para_2(15) = 1;
para_2(16) = 0.01*10.5/0.5*1e-16;
para_2(17) = 10;
para_2(18) = 0;

[ttt, zzz]=ode15s(@pathway_model_many_antibody_immune_res_new,[0 10],x0,[],para_2,para_new);

out_IgG(k) = 1e16-interp1(ttt,zzz(:,101),10);
end

plot(out_IgG);
hold on

clear x0
%% 
for k = 1:101
 for i = 1:100
    xx0(i) = data_new(i,k);
 end  
 x0(1) = sum(xx0);
 x0(2) = 1e15;
 x0(3) = 0;
 para(1) = 10^(-13.5);
 para(2) = 10^(-0.5);
 
[tt, zz]=ode15s(@pathway_model_many_serum_overall,[0 50],x0,[],para);

percent(k) = interp1(tt,zz(:,1),50)/x0(1);

 
    
end


for k = 1:101
    
   
theta(1) = -15.5;
for i = 1:100
    x0(i) = data_new(i,k)*percent(k);
end

x0(101) = 1e16;%% virus

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+101) = 0;
    end
end

x0(202) = 0;%% environmental antigen concentration

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+202) = 0;
    end
end

para_2(1) = 10^(theta(1)-4.5); 
para_2(2) = 10^(theta(1)-3.5); 
para_2(3) = 10^(theta(1)-2.5); 
para_2(4) = 10^(theta(1)-1.5); 
para_2(5) = 10^(theta(1)-0.5);  
para_2(6) = 10^(theta(1)+0.5); 
para_2(7) = 10^(theta(1)+1.5);  
para_2(8) = 10^(theta(1)+2.5); 
para_2(9) = 10^(theta(1)+3.5);  
para_2(10) = 10^(theta(1)+4.5); 

para_new(1) = 1e-3; 
para_new(2) = 1e-2;
para_new(3) = 1e-1; 
para_new(4) = 1e0;
para_new(5) = 1e1; 
para_new(6) = 1e2;
para_new(7) = 1e3; 
para_new(8) = 1e4;
para_new(9) = 1e5; 
para_new(10) = 1e6;


para_2(11) = 0;%% 0.01
para_2(12) = 0;%% 5
para_2(13) = 0;%% 1
para_2(14) = 0;%% 0.5
para_2(15) = 1;
para_2(16) = 0.01*10.5/0.5*1e-16;
para_2(17) = 10;
para_2(18) = 0;


[ttt, zzz]=ode15s(@pathway_model_many_antibody_immune_res_new,[0 10],x0,[],para_2,para_new);

out_IgG(k) = 1e16-interp1(ttt,zzz(:,101),10);
end

plot(out_IgG);
hold on
clear x0
%%
for k = 1:101
 for i = 1:100
    xx0(i) = data_new(i,k);
 end  
 x0(1) = sum(xx0);
 x0(2) = 2e15;
 x0(3) = 0;
 para(1) = 10^(-13.5);
 para(2) = 10^(-0.5);
 
[tt, zz]=ode15s(@pathway_model_many_serum_overall,[0 50],x0,[],para);

percent(k) = interp1(tt,zz(:,1),50)/x0(1);

 
    
end


for k = 1:101
    
   
theta(1) = -15.5;
for i = 1:100
    x0(i) = data_new(i,k)*percent(k);
end

x0(101) = 1e16;%% virus

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+101) = 0;
    end
end

x0(202) = 0;%% environmental antigen concentration

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+202) = 0;
    end
end


para_2(1) = 10^(theta(1)-4.5); 
para_2(2) = 10^(theta(1)-3.5); 
para_2(3) = 10^(theta(1)-2.5); 
para_2(4) = 10^(theta(1)-1.5); 
para_2(5) = 10^(theta(1)-0.5);  
para_2(6) = 10^(theta(1)+0.5); 
para_2(7) = 10^(theta(1)+1.5);  
para_2(8) = 10^(theta(1)+2.5); 
para_2(9) = 10^(theta(1)+3.5);  
para_2(10) = 10^(theta(1)+4.5); 

para_new(1) = 1e-3; 
para_new(2) = 1e-2;
para_new(3) = 1e-1; 
para_new(4) = 1e0;
para_new(5) = 1e1; 
para_new(6) = 1e2;
para_new(7) = 1e3; 
para_new(8) = 1e4;
para_new(9) = 1e5; 
para_new(10) = 1e6;


para_2(11) = 0;%% 0.01
para_2(12) = 0;%% 5
para_2(13) = 0;%% 1
para_2(14) = 0;%% 0.5
para_2(15) = 1;
para_2(16) = 0.01*10.5/0.5*1e-16;
para_2(17) = 10;
para_2(18) = 0;


[ttt, zzz]=ode15s(@pathway_model_many_antibody_immune_res_new,[0 10],x0,[],para_2,para_new);

out_IgG(k) = 1e16-interp1(ttt,zzz(:,101),10);
end

plot(out_IgG);
hold on
clear x0

%%
for k = 1:101
 for i = 1:100
    xx0(i) = data_new(i,k);
 end  
 x0(1) = sum(xx0);
 x0(2) = 5e15;
 x0(3) = 0;
 para(1) = 10^(-13.5);
 para(2) = 10^(-0.5);
 
[tt, zz]=ode15s(@pathway_model_many_serum_overall,[0 50],x0,[],para);

percent(k) = interp1(tt,zz(:,1),50)/x0(1);
    
end


for k = 1:101
    
   
theta(1) = -15.5;
for i = 1:100
    x0(i) = data_new(i,k)*percent(k);
end

x0(101) = 1e16;%% virus

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+101) = 0;
    end
end

x0(202) = 0;%% environmental antigen concentration

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+202) = 0;
    end
end


para_2(1) = 10^(theta(1)-4.5); 
para_2(2) = 10^(theta(1)-3.5); 
para_2(3) = 10^(theta(1)-2.5); 
para_2(4) = 10^(theta(1)-1.5); 
para_2(5) = 10^(theta(1)-0.5);  
para_2(6) = 10^(theta(1)+0.5); 
para_2(7) = 10^(theta(1)+1.5);  
para_2(8) = 10^(theta(1)+2.5); 
para_2(9) = 10^(theta(1)+3.5);  
para_2(10) = 10^(theta(1)+4.5); 

para_new(1) = 1e-3; 
para_new(2) = 1e-2;
para_new(3) = 1e-1; 
para_new(4) = 1e0;
para_new(5) = 1e1; 
para_new(6) = 1e2;
para_new(7) = 1e3; 
para_new(8) = 1e4;
para_new(9) = 1e5; 
para_new(10) = 1e6;


para_2(11) = 0;%% 0.01
para_2(12) = 0;%% 5
para_2(13) = 0;%% 1
para_2(14) = 0;%% 0.5
para_2(15) = 1;
para_2(16) = 0.01*10.5/0.5*1e-16;
para_2(17) = 10;
para_2(18) = 0;


[ttt, zzz]=ode15s(@pathway_model_many_antibody_immune_res_new,[0 10],x0,[],para_2,para_new);

out_IgG(k) = 1e16-interp1(ttt,zzz(:,101),10);
end
plot(out_IgG);
hold on
