clc
clear

global R;
R = 6375 * 10^3;
global h_atm;
h_atm = 50 * 10^3;
global h_trop;
h_trop = 18 * 10^3;

global F_pull;
F_pull = 1 * 10^2;
global t_pull;
t_pull = 600;

global x_0;
x_0 = 0;
global y_0;
y_0 = R + 0; %высота над центром Земли
global alpha;
alpha = pi/3;
global len;
len = 10000;

g_0 = g(R);
t = 0;
x = x_0;
y = y_0;
v = 0;
v_x = v * cos(alpha);
v_y = v * sin(alpha);

if alpha == - pi/2 || alpha == pi/2 %нужно т.к. cos(pi/2) не даёт точный 0
    v_x = 0;
end

h = 10^(-1);

T = [0 0];
X = [0 0];
Y = [0 0];
V_x = [0 0];
V_y = [0 0];
W_y = [0,0];

L = [0 0];

T(1) = t;
X(1) = x;
Y(1) = y;
V_x(1) = v_x;
V_y(1) = v_y;

i = 2;

F = figure(1);
A = axes(F);
P = plot(A, X, Y, 'LineWidth',1);
xlabel( 'x' ); ylabel( 'y' );
hold on
%xlim([-8 * 10^6, 8 * 10^6])
%ylim([-8 * 10^6, 8 * 10^6])

U = 0 : pi/200 : 2*pi;
plot(A, R*cos(U), R*sin(U), 'r');

plot(A, (R + h_atm)*cos(U), (R + h_atm)*sin(U), 'B--');
plot(A, (R + h_trop)*cos(U), (R + h_trop)*sin(U), '--');

plot(A, 0, 0, 'o');

t_1 = t;
h_1 = 0.1 * 10^(2);
epsilon = 0.2 * 10^(-2);
t_lim = 1.139000000020220 * 10^5;

while t < t_lim
    while true
                
        if t_1 - t <= epsilon
            pause(0.001);
            P.XData = X;
            P.YData = Y;   
            drawnow
            t_1 = t_1 + h_1;
        end      
        
        t_prev = t;
        x_prev = x;
        y_prev = y;
        v_x_prev = v_x;
        v_y_prev = v_y;
        
        k1 = k_1(t, x, v_x, y, v_y);
        l1 = l_1(t, x, v_x, y, v_y);
        m1 = m_1(t, x, v_x, y, v_y);
        n1 = n_1(t, x, v_x, y, v_y);
        
        k2 = k_2(t, x, v_x, y, v_y, h, k1, l1, m1, n1);
        l2 = l_2(t, x, v_x, y, v_y, h, k1, l1, m1, n1);
        m2 = m_2(t, x, v_x, y, v_y, h, k1, l1, m1, n1);
        n2 = n_2(t, x, v_x, y, v_y, h, k1, l1, m1, n1);
        
        k3 = k_3(t, x, v_x, y, v_y, h, k2, l2, m2, n2);
        l3 = l_3(t, x, v_x, y, v_y, h, k2, l2, m2, n2);
        m3 = m_3(t, x, v_x, y, v_y, h, k2, l2, m2, n2);
        n3 = n_3(t, x, v_x, y, v_y, h, k2, l2, m2, n2);
        
        k4 = k_4(t, x, v_x, y, v_y, h, k3, l3, m3, n3);
        l4 = l_4(t, x, v_x, y, v_y, h, k3, l3, m3, n3);
        m4 = m_4(t, x, v_x, y, v_y, h, k3, l3, m3, n3);
        n4 = n_4(t, x, v_x, y, v_y, h, k3, l3, m3, n3);
        
        x = x + h/6 * (k1 + 2 * k2 + 2 * k3 + k4);
        v_x = v_x + h/6 * (l1 + 2 * l2 + 2 * l3 + l4);
        y = y + h/6 * (m1 + 2 * m2 + 2 * m3 + m4);
        v_y = v_y + h/6 * (n1 + 2 * n2 + 2 * n3 + n4);
        w_y = g_y(t, x, v_x, y, v_y);
        
        t = t + h;

        if abs(t - t_pull) <= epsilon
            text = ['The pull stopped, t = ', num2str(t), ' sec'];
            disp(text)
            plot(A, x, y, 'kx');
        end
        
        if t > t_lim
            break;
        end
        
        T(i) = t;
        X(i) = x;
        Y(i) = y;

        i = i + 1;

        if sqrt(x^2 + y^2) <= R
            text = ['The object landed, t = ', num2str(t), ' sec'];
            disp(text)
            
            P.XData = X;
            P.YData = Y;
            drawnow
            
            plot(A, x, y, 'kx');
            return;
        end
    end
end


function vel_x = f_x(t, x, v_x, y, v_y)
    vel_x = v_x;
end

function vel_y = f_y(t, x, v_x, y, v_y)
    vel_y = v_y;
end

function theta = theta(x, y)
    if x >= 0
        theta = atan(y / x);
    elseif x < 0
        theta = pi + atan(y / x);
    elseif x == 0 && y > 0
        theta = pi/2;
    elseif x == 0 && y < 0
        theta = -pi/2;
    elseif x == 0 && y == 0
        theta = 0;
    end
end

function w_x = g_x(t, x, v_x, y, v_y)
    global R;
    global h_atm;
    global F_pull;
    global t_pull;
    global x_0;
    global y_0;
    global alpha;
    global len;
    
    r = sqrt(x^2 + y^2);
    C_D = 0.47;
    S = pi * 0.01^2;
    
    if t <= t_pull
        F = F_pull;
    else
        F = 0;
    end
    
    if sqrt((x - x_0)^2 + (y - y_0)^2) <= len
        w_x = F * cos(alpha);
    else
        w_x_ord = - g(r) * cos(theta(x, y)) - v_x * m(t) / M(t) + v_x / sqrt(v_x^2 + v_y^2) * F / M(t);
        if r - R <= h_atm
            w_x = w_x_ord - 1/2 * rho(x, y) * C_D * S * v_x * sqrt(v_x^2 + v_y^2) / M(t);
        else
            w_x = w_x_ord;
        end
    end
end

function w_y = g_y(t, x, v_x, y, v_y)
    global R;
    global h_atm;
    global F_pull;
    global t_pull;
    global x_0;
    global y_0;
    global alpha;
    global len;
    
    r = sqrt(x^2 + y^2);
    C_D = 0.47;
    S = pi * 0.01^2;
    
    if t <= t_pull
        F = F_pull;
    else
        F = 0;
    end
    
    if sqrt((x - x_0)^2 + (y - y_0)^2) <= len
        w_y = F * sin(alpha);
    else
        w_y_ord = - g(r) * sin(theta(x, y)) - v_y * m(t) / M(t) + v_y / sqrt(v_x^2 + v_y^2) * F / M(t);
        if r - R <= h_atm
            w_y = w_y_ord - 1/2 * rho(x, y) * C_D * S * v_y * sqrt(v_x^2 + v_y^2) / M(t);
        else
            w_y = w_y_ord;
        end
    end
end

function g = g(r)
    G = 6.67430 * 10^(-11);
    M = 5.9722 * 10^24;
    g = G * M / r^2;
end

function M = M(t)
global t_pull;
    if t <= t_pull
        M = 50 - 50/t_pull * t + 10;
    else
        M = 10;
    end
end

function m = m(t)
    dt = 10^(-12);
    m = (M(t + dt) - M(t - dt)) / (2 * dt);
end

function rho = rho(x, y)
    global R;
    global h_atm;
    global h_trop;
    p_0 = 101325;
    T_0 = 288.15;
    L = 0.0065;
    R_gas = 8.31446;
    mu = 0.0289652;
    H_tr = 6300;
    
    r = sqrt(x^2 + y^2);
    
    if r - R <= h_trop
        rho = p_0 * mu / (R_gas * T_0) * (1 - L/T_0 * (r - R))^(mu * g(r) / (L * R_gas) - 1);
    elseif r - R > h_trop && r - R <= h_atm
        rho = p_0 * mu / (R_gas * T_0) * (1 - L/T_0 * h_trop)^(mu * g(R + h_trop) / (L * R_gas) - 1) * exp(- (r - R - h_trop) / H_tr);
    elseif r - R > h_atm
        rho = 0;
    end
end


function k = k_1(t, x, v_x, y, v_y)
    k = f_x(t, x, v_x, y, v_y);
end

function l = l_1(t, x, v_x, y, v_y)
    l = g_x(t, x, v_x, y, v_y);
end

function m = m_1(t, x, v_x, y, v_y)
    m = f_y(t, x, v_x, y, v_y);
end

function n = n_1(t, x, v_x, y, v_y)
    n = g_y(t, x, v_x, y, v_y);
end

function k = k_2(t, x, v_x, y, v_y, h, k_1, l_1, m_1, n_1)
    k = f_x(t + h/2, x + h/2 * k_1, v_x + h/2 * l_1, y + h/2 * m_1, v_y + h/2 * n_1);
end

function l = l_2(t, x, v_x, y, v_y, h, k_1, l_1, m_1, n_1)
    l = g_x(t + h/2, x + h/2 * k_1, v_x + h/2 * l_1, y + h/2 * m_1, v_y + h/2 * n_1);
end

function m = m_2(t, x, v_x, y, v_y, h, k_1, l_1, m_1, n_1)
    m = f_y(t + h/2, x + h/2 * k_1, v_x + h/2 * l_1, y + h/2 * m_1, v_y + h/2 * n_1);
end

function n = n_2(t, x, v_x, y, v_y, h, k_1, l_1, m_1, n_1)
    n = g_y(t + h/2, x + h/2 * k_1, v_x + h/2 * l_1, y + h/2 * m_1, v_y + h/2 * n_1);
end

function k = k_3(t, x, v_x, y, v_y, h, k_2, l_2, m_2, n_2)
    k = f_x(t + h/2, x + h/2 * k_2, v_x + h/2 * l_2, y + h/2 * m_2, v_y + h/2 * n_2);
end

function l = l_3(t, x, v_x, y, v_y, h, k_2, l_2, m_2, n_2)
    l = g_x(t + h/2, x + h/2 * k_2, v_x + h/2 * l_2, y + h/2 * m_2, v_y + h/2 * n_2);
end

function m = m_3(t, x, v_x, y, v_y, h, k_2, l_2, m_2, n_2)
    m = f_y(t + h/2, x + h/2 * k_2, v_x + h/2 * l_2, y + h/2 * m_2, v_y + h/2 * n_2);
end

function n = n_3(t, x, v_x, y, v_y, h, k_2, l_2, m_2, n_2)
    n = g_y(t + h/2, x + h/2 * k_2, v_x + h/2 * l_2, y + h/2 * m_2, v_y + h/2 * n_2);
end

function k = k_4(t, x, v_x, y, v_y, h, k_3, l_3, m_3, n_3)
    k = f_x(t + h, x + h * k_3, v_x + h * l_3, y + h * m_3, v_y + h * n_3);
end

function l = l_4(t, x, v_x, y, v_y, h, k_3, l_3, m_3, n_3)
    l = g_x(t + h, x + h * k_3, v_x + h * l_3, y + h * m_3, v_y + h * n_3);
end

function m = m_4(t, x, v_x, y, v_y, h, k_3, l_3, m_3, n_3)
    m = f_y(t + h, x + h * k_3, v_x + h * l_3, y + h * m_3, v_y + h * n_3);
end

function n = n_4(t, x, v_x, y, v_y, h, k_3, l_3, m_3, n_3)
    n = g_y(t + h, x + h * k_3, v_x + h * l_3, y + h * m_3, v_y + h * n_3);
end
