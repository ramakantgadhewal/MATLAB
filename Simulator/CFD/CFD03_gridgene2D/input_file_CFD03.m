include_flags_CFD03;

% 给定翼面上网格点数目，要求为偶数
num = 100;
% 计算翼面网格点坐标(x0,y0)
[x0,y0] = airfoildata;

% 模型的弦长
len  = max(x0)-min(x0);
% 最远距离除以弦长
fac  = 15;
% 代数网格生成指数
mag  = 1;

% 网格数：ksi[i=const,1:2*m+num+1], eta[j=const,1:n]
m = 50;
n = 100;
% Poisson 方程源项系数
% - Format: (fgksi, ksi, fgeta, eta, amplitude, decay factor)
attr_num = 3;
attr = [1,1,      0,1,    100,10
        1,m+1,    1,1,    10000,10
        1,m+num+1,1,1,    10,10];
% 迭代指标：粘合点是(1)/否(0)参与网格迭代
iterflag_glueline = 0;
% 迭代指标：尾迹区是(1)/否(0)参与网格迭代
iterflag_mode = 1;
% 迭代指标：求解过程之中网格结点误差限
qual = 1e-2;