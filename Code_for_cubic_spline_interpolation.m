format short;
syms x1 x2 x3 x4 x5 x6 x7;
x = [x1,x2,x3,x4,x5,x6,x7];
syms f f1 f2;
R = [0.02832*0.5,0.02738*0.5,0.02689*0.5,0.02665*0.5,0.02668*0.5,0.02702*0.5,0.02777*0.5,0.02882*0.5];
r = R.*2;
yr = [0.5,2,3,4,5,7,10,30];
num = yr.*2;
dyr = []; 
for t = 2:length(yr)
   dyr(t-1) = (yr(t) - yr(t-1))/0.5; 
end
num(1) = [];
dnum = [];
for i = 2:length(num)
   dnum(i-1) = num(i) - num(i-1);
end

%solve x1 ---- 2yr spot rate
f1 = (1/(1+R(1))+1/(1+0.5*(2*r(1)+x1)/3)^2+1/(1+0.5*(r(1)+2*x1)/3)^3+1/(1+0.5*x1)^4);
f2 = 1/(1+0.5*x1)^4;
f = R(2)*f1+f2;
x1 = vpasolve(1==f,x1,[0,Inf]);
x(1) = x1;
f1 = (1/(1+R(1))+1/(1+0.5*(2*r(1)+x1)/3)^2+1/(1+0.5*(r(1)+2*x1)/3)^3+1/(1+0.5*x1)^4);

%solve xi ---- i=2:6
for j = 2:length(x)-1
    for k = 1:dnum(j-1)
      f1 = f1 + 1/(1+0.5*((dnum(j-1)-k)*x(j-1)+k*x(j))/dnum(j-1))^(num(j-1)+k);
    end
    f2 = 1/(1+0.5*x(j))^(num(j));
    f = R(j+1)*f1 + f2;
    x(j) = vpasolve(1==f,x(j),[0,Inf]);
    f1 = subs(f1,x(j));
end

%------------------%----------------%-------------%---------------%-------------%------------%-------%
% 
%%
% *For solving the x7 --- 30yr spot rate, it costs much much longer time than before*

for j = 7
    for k = 1:dnum(j-1)
      f1 = f1 + 1/(1+0.5*((dnum(j-1)-k)*x(j-1)+k*x(j))/dnum(j-1))^(num(j-1)+k);
    end
    f2 = 1/(1+0.5*x(j))^(num(j));
    f = R(j+1)*f1 + f2;
    x(j) = vpasolve(1==f,x(j),[0,Inf]);
end