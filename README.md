# UAV-Dynamical-Model-with-Euler-Lagrange-Approach
SUSTech SDM271 System Modeling and Simulation

Build the dynamics model of the UAV in Matlab/simulink, and realize the UAV heavy straight upward motion for 5 seconds by controlling four motors.

```
四旋翼构型为“X”型，螺旋桨序号如下所示：
         3↓   1↑
           \ /
           / \
         2↑   4↓
其中，↑表示螺旋桨逆时针旋转；↓表示螺旋桨顺时针旋转。
```

## 外部封装

![image](https://user-images.githubusercontent.com/117464811/232400041-a02f1691-366f-4566-a015-334dc7a0e8bb.jpeg)

## 总体结构

控制效果模型、动力学模型、运动学模型

![image](https://user-images.githubusercontent.com/117464811/232400144-3a323105-58ab-4e66-9777-2d0e6c4ffe72.jpeg)

### 控制效果模型

![image](https://user-images.githubusercontent.com/117464811/232400240-b8c5aa4b-0ac8-4100-b623-58180f34edf0.jpeg)

```matlab
% 函数描述
% 1.作用：本函数用来计算螺旋桨旋转产生的总拉力和反扭力矩。

function [f, tau_x, tau_y, tau_z] = fcn(w1, w2, w3, w4, c_T, c_M, d)

% 2.函数输入：
%   wi：四个螺旋桨的转速(rad/s)
%   c_T：螺旋桨拉力系数
%   c_M：螺旋桨转矩系数
%   d：机体中心和任一电机的距离(m)

% 3.函数输出：
%   f：螺旋桨拉力（机体轴）
%   tau_x：x轴反扭力矩（机体轴）
%   tau_y：y轴反扭力矩（机体轴）
%   tau_z：z轴反扭力矩（机体轴）

    f = c_T * (w1^2 + w2^2 + w3^2 + w4^2);
    tau_x = d * c_T * (sqrt(2)/2) * (-w1^2 + w2^2 + w3^2 - w4^2);
    tau_y = d * c_T * (sqrt(2)/2) * (w1^2 - w2^2 + w3^2 - w4^2);
    tau_z = c_M * (w1^2 + w2^2 - w3^2 - w4^2);
    
end
```

### 动力学模型

位置动力学模型、姿态动力学模型

![image](https://user-images.githubusercontent.com/117464811/232400341-023d7cfb-9e88-46f8-b208-f4bddd26f33a.jpeg)

#### 位置动力学模型

![image](https://user-images.githubusercontent.com/117464811/232400393-1ff23765-e19f-499c-94db-08cab5e79c7b.jpeg)

```matlab
% 函数描述
% 1.作用：本函数为四旋翼的位置动力学微分方程组，通过四旋翼所受拉力、姿态角计算得到四旋翼飞行器的加速度

function [v_x_dot, v_y_dot, v_z_dot] = fcn(m, g, phi, theta, psi, f)

    % 2.函数输入：
    %   m：四旋翼飞行器质量(kg)
    %   g：重力加速度(m/s^2)
    %   phi：滚转角(rad)
    %   theta：俯仰角(rad)
    %   psi：偏航角(rad)
    %   f：螺旋桨产生的总拉力(N)
    
    %   v_x_dot：地球坐标系下沿x轴的速度的导数
    v_x_dot = -f * (1/m) * (cos(psi) * sin(theta) * cos(phi) + sin(psi) * sin(phi));
    
    %   v_y_dot：地球坐标系下沿y轴的速度的导数
    v_y_dot = -f * (1/m) * (sin(psi) * sin(theta) * cos(phi) - cos(psi) * sin(phi));
    
    %   v_z_dot：地球坐标系下沿z轴的速度的导数
    v_z_dot = g - f * (1/m) * cos(phi) * cos(theta);

end
```

#### 姿态动力学模型

![image](https://user-images.githubusercontent.com/117464811/232400453-02f10743-831f-4939-bb82-9a3da6113299.jpeg)

```matlab
% 函数描述
% 1.作用：本函数为四旋翼的姿态动力学微分方程组，通过四旋翼所受力矩和螺旋桨转速计算得到四旋翼的机体角速度

function [p_dot, q_dot, r_dot] = ...
fcn(tau_x, tau_y, tau_z, w1, w2, w3, w4, I_xx, I_yy, I_zz, J_RP, p, q, r)
    
    % 2.函数输入：
    %   wi：四个螺旋桨的转速(rad/s)
    %   tau_x：x轴反扭力矩（机体轴）(N·m)
    %   tau_y：y轴反扭力矩（机体轴）(N·m)
    %   tau_z：z轴反扭力矩（机体轴）(N·m)
    %   I_xx：四旋翼x轴转动惯量(kg·m^2)
    %   I_yy：四旋翼y轴转动惯量(kg·m^2)
    %   I_zz：四旋翼z轴转动惯量(kg·m^2)
    %   J_RP：整个电机转子和螺旋桨绕转轴的总转动惯量(kg·m^2)
  
    Omega = -w1 + w2 - w3 + w4;

    %   p_dot：四旋翼x轴角加速度（机体轴)
    %   p：四旋翼x轴角速度（机体轴)(rad/s)
    p_dot = (1/I_xx) * (tau_x + q * r * (I_yy - I_zz) - J_RP * q * Omega);

    %   q_dot：四旋翼y轴角加速度（机体轴)
    %   q：四旋翼y轴角速度（机体轴)(rad/s)
    q_dot = (1/I_yy) * (tau_y + p * r * (I_zz - I_xx) + J_RP * p * Omega);

    %   r_dot：四旋翼z轴角加速度（机体轴)
    %   r：四旋翼z轴角速度（机体轴)(rad/s)
    r_dot = (1/I_zz) * (tau_z + p * q * (I_xx - I_yy) );

end
```

### 运动学模型

位置运动学模型、姿态运动学模型

![image](https://user-images.githubusercontent.com/117464811/232400704-9e7bf79b-8248-4d57-889b-4b2071d31630.jpeg)

#### 位置运动学模型

![image](https://user-images.githubusercontent.com/117464811/232400778-64b253e8-e537-4581-97bf-6d043f4749d6.jpeg)

#### 姿态运动学模型

![image](https://user-images.githubusercontent.com/117464811/232400824-65f98fee-01fc-468c-8457-c21329b804d7.jpeg)

```matlab
% 函数描述
% 1.作用：本函数为四旋翼的姿态运动学微分方程组，通过四旋翼的机体角速度计算得到四旋翼的姿态角

function [phi_dot, theta_dot, psi_dot] = fcn(p, q, r, phi, theta)
    
    % 2.函数输入：
    %   p：四旋翼x轴角速度（机体轴)(rad/s)
    %   q：四旋翼y轴角速度（机体轴)(rad/s)
    %   r：四旋翼z轴角速度（机体轴)(rad/s)
    %   phi：滚转角(rad)
    %   theta：俯仰角(rad)

    %   phi_dot：滚转角速度
    phi_dot = p + q * tan(theta) * sin(phi) + r * tan(theta) * cos(phi);

    %   theta_dot：俯仰角速度
    theta_dot = q * cos(phi) - r * sin(phi);

    %   psi_dot：偏航角速度
    psi_dot = (q * sin(phi) + r * cos(phi)) / cos(theta);

end
```

## 运行结果

### 位置变化曲线

![image](https://user-images.githubusercontent.com/117464811/232401047-00438cf8-f45e-4151-8d7b-a7bb1ee523f1.jpeg)

### 平移速度变化曲线

![image](https://user-images.githubusercontent.com/117464811/232401071-6ae52a0f-9f0b-4b51-8419-f85501a95af3.jpeg)

### 姿态角变化曲线

![image](https://user-images.githubusercontent.com/117464811/232401166-02cddcf2-95f0-41cd-a3ea-ebcd4f5d808f.jpeg)

