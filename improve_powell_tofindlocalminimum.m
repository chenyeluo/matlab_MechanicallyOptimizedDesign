% 定义函数
f = @(x) (x-2)^2;%x^3 - 2*x^2;
% 初始点
x0 = 0;  

% 调用算法
[xmin, fmin, iter] = improved_powell(f, x0);

% 显示结果
disp(['极小值点: x = ', num2str(xmin)]);
disp(['极小值: f(x) = ', num2str(fmin)]);
disp(['迭代次数: ', num2str(iter)]);




%%
% 改进Powell算法
function [xmin, fmin, iter] = improved_powell(f, x0, tol, maxiter)
% 输入:
%   f - 目标函数句柄
%   x0 - 初始点
%   tol - 容差（默认1e-6）
%   maxiter - 最大迭代次数（默认100）
% 输出:
%   xmin - 找到的最小值点
%   fmin - 最小值
%   iter - 实际迭代次数

if nargin < 4
    maxiter = 100;
    if nargin < 3
        tol = 1e-6;
    end
end

x = x0;
fmin = f(x);
iter = 0;
delta = inf;

% 初始方向集（前进和后退方向）
D = [1, -1]; 

while iter < maxiter && delta > tol
    iter = iter + 1;
    fold = fmin;
    
    % 沿两个方向搜索
    for i = 1:2
        d = D(i);
        f1d = @(alpha) f(x + alpha*d);
        alpha = golden_search(f1d, 0, 1, 1e-6);
        x = x + alpha*d;
        fmin = f(x);
    end
    
    delta = abs(fold - fmin);
    
    % 构造新方向
    d_new = x - (x - D(1)*alpha);
    
    % 检查方向替换条件
    f1 = f(x - D(1)*alpha);
    f2 = fmin;
    f3 = f(x + d_new);
    
    if (f1 - 2*f2 + f3)*(f1 - f2 - delta)^2 >= 0.5*delta*(f1 - f3)^2
        D = [d_new, -d_new]; % 替换方向集
        f1d = @(alpha) f(x + alpha*d_new);
        alpha = golden_search(f1d, 0, 1, 1e-6);
        x = x + alpha*d_new;
        fmin = f(x);
    end
end

xmin = x;
end

function alpha = golden_search(f, a, b, tol)
% 黄金分割法一维搜索
gr = (sqrt(5)-1)/2;

c = b - gr*(b-a);
d = a + gr*(b-a);

while abs(c-d) > tol
    if f(c) < f(d)
        b = d;
    else
        a = c;
    end
    
    c = b - gr*(b-a);
    d = a + gr*(b-a);
end

alpha = (a+b)/2;
end