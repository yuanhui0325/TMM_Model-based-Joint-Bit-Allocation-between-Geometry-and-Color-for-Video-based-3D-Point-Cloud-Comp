function [target,modelPara] = TMM_convex_method(GQstep,CQstep,Geo_Bits,Col_Bits,MSE,targetRate)
% -----------------------------------------------------------------------------------
% input:
% GQstep: 1 x 3 vector 预编码的几何量化步长
% CQstep: 1 x 3 vector 预编码的属性量化步长
% Geo_Bits: 1 x 3 vector 预编码的几何比特(bits per million points)
% Col_Bits: 1 x 3 vector 预编码的属性比特(bits per million points)
% MSE: 1 x 3 vector 预编码的失真 颜色和几何失真总和 normal
% targetRange: 待搜索目标码率的区间 eg. targetRange = linspace(62000,690000,315)
% -----------------------------------------------------------------------------------
% output:
% target: m x n vector 每一行表示 目标码率 目标QP 目标量化步长
% modelPara: 1 x 4 vector 模型参数 顺序为a b c d
% -----------------------------------------------------------------------------------
target = [];
QP=[22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42];
Qstep=[8 9 10 11.25 12.75 14.25 16 18 20 22.5 25.5 28.5 32 36 40 45 51 57 64 72 80];
Model_GQstep=Qstep;
Model_CQstep=Qstep;
clear g c;
%%%%% %解方程
%Alex的Geo_R=a*power(GQstep,b) 使用(28.5,11.25),(32,36)。
geo_p=[28.5 Geo_Bits(1)
    32 Geo_Bits(2)];
geo_x0=[300000,-1];
option=optimset('MaxIter',100000,'MaxFunEvals',1000000);
%geo_x(1)对应参数a;geo_x(2)对应参数b.
[geo_x,geo_fval,geo_exitflag,geo_output]=fsolve(@(x)seperateratefun(x,geo_p),geo_x0,option);

%Alex的Col_R=c*power(CQstep,d) 使用(28.5,11.25),(32,36)。
col_p=[11.25 Col_Bits(1)
    36 Col_Bits(2)];
col_x0=[10000000,-2];
option=optimset('MaxIter',100000,'MaxFunEvals',1000000);
%col_x(1)对应参数c;col_x(2)对应参数d.
[col_x,col_fval,col_exitflag,col_output]=fsolve(@(x)seperateratefun(x,col_p),col_x0,option);

%Alex的D=d*GQstep+e*CQstep+f使用(28.5,11.25),(32,36),(10,28.5),
syms d e f d1 d2 d3;
d1=MSE(1);
d2=MSE(2);
d3=MSE(3);
equ4=d*28.5+e*11.25+f-d1;
equ5=d*32+e*36+f-d2;
equ6=d*10+e*28.5+f-d3;
[d,e,f]=solve(equ4,equ5,equ6);
d=vpa(d);
e=vpa(e);
f=vpa(f);
%%%%% %解方程
modelPara = [geo_x(1)  geo_x(2) col_x(1) col_x(2) d e f]; %Geo_R=a*power(GQstep,b) ; Col_R=c*power(CQstep,d); D=d*GQstep+e*CQstep+f
%Model data optimal convex operation with Interior-point
for j = 1 : length(targetRate)
    options = optimoptions('fmincon','Algorithm','interior-point', 'MaxIter',1000000, 'MaxFunEvals',1000000);
    A_value=[]; b_value=[];Aeq=[];beq=[];lb=[8,8];ub=[80,80];
    n=[geo_x(1),geo_x(2),col_x(1),col_x(2),targetRate(j)];
    m=[double(d),double(e),double(f)];
    %x(1)表示Geometry Qstep;x(2)表示Color Qstep.
    [x,fval,exitflag] = fmincon(@(x) fmincon_dis(x,m),[8,8],A_value,b_value,Aeq,beq,lb,ub,@(x) fmincon_rate(x,n),options); %x(1)表示GQstep，x(2)表示CQstep
    %找到离散的（GQstep,CQstep)中离最优点（x(1),x(2))最近的点作为最优点
    %初始化最优点的初始值
    min_dis=sqrt(power((80-8),2)+power((80-8),2));
    convex_GQstep=0;
    convex_CQstep=0;
    for GQstep_index=1:21
        for CQstep_index=1:21
            if sqrt(power((Model_GQstep(GQstep_index)-x(1)),2)+power((Model_CQstep(CQstep_index)-x(2)),2))<=min_dis
                min_dis=sqrt(power((Model_GQstep(GQstep_index)-x(1)),2)+power((Model_CQstep(CQstep_index)-x(2)),2));
                convex_GQstep=Model_GQstep(GQstep_index);
                convex_CQstep=Model_CQstep(CQstep_index);
            end
        end
    end
    % 此处将每个target rate得到的几何 属性量化参数信息迭代放进输出 可参考下面的格式
    g_i=find(Qstep==convex_GQstep);
    gQP=QP(g_i);
    c_i=find(Qstep==convex_CQstep);
    cQP=QP(c_i);
    temp= [targetRate(j), gQP, cQP]; %行向量
    target = [target; temp]; %按行插入到输出中
    clear g_i gQP c_i cQP convex_GQstep convex_CQstep;
end
clear Model_GQstep Model_CQstep n m x convex_CQstep convex_GQstep;