function y = fun(x)
%函数用于计算粒子适应度值
%x     ·      input           输入粒子 
%y           output          粒子适应度值 
target = [155,258,98,106,123,135,138,140,150,180];%10个航班的target时间，用于计算适应度函数
c = [10,10,30,30,30,30,30,30,30,30];%10个航班的早到/晚到的惩罚，用于计算适应度函数
time_span = [0	3	15	15	15	15	15	15	15	15;
3	0	15	15	15	15	15	15	15	15;
15	15	0	8	8	8	8	8	8	8;
15	15	8	0	8	8	8	8	8	8;
15	15	8	8	0	8	8	8	8	8;
15	15	8	8	8	0	8	8	8	8;
15	15	8	8	8	8	0	8	8	8;
15	15	8	8	8	8	8	0	8	8;
15	15	8	8	8	8	8	8	0	8;
15	15	8	8	8	8	8	8	8	0];%计算时间惩罚
cost = 0.0;
xx = x';
for i=1:10
    %disp(x(i))
    cost = (abs(xx(i)-target(i)))*c(i) + cost;
    for j=1:10
        if (abs(xx(i)-xx(j))-time_span(i,j))<0
            cost = Inf;
        end
    end
end
y=cost;


