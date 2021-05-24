clc
clear

t = 0;
v = 10;
x = 0;
y = 10;
alpha = pi/3;
v_x = v * cos(alpha);
v_y = v * sin(alpha);

if alpha == - pi/2 || alpha == pi/2
    v_x = 0;
elseif alpha == 0
    v_x = v;
end

l = 0;

m = 1;

h = 0.001;

T = [0 0];
X = [0 0];
Y = [0 0];
V_x = [0 0];
V_y = [0 0];

L = [0 0];

T(1) = t;
X(1) = x;
Y(1) = y;
V_x(1) = v_x;
V_y(1) = v_y;

L(1) = l;

i = 2;

F = figure(1);
A = axes(F);
P = plot(A, X, Y, '-', 'LineWidth',1);
xlabel( 'x' ); ylabel( 'y' );
hold on
xlim([-20,20])
% ylim([0,14])

x_land = -20:0.1:20;
plot(A, x_land, land(x_land),'LineWidth',1);

t_1 = t;
h_1 = 0.01;
epsilon = 0.2 * 10^(-2);
t_lim = 30;

while t < t_lim
    while y >= land(x)
                
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
        
        if y < land(x)
            disp(t)
            y = land(x);
            
            P.XData = X;
            P.YData = Y;
            drawnow
            
            break;
        end

        if t > t_lim
            return;
        end
        
        T(i) = t;
        X(i) = x;
        Y(i) = y;
        V_x(i) = v_x;
        V_y(i) = v_y;
        i = i + 1;
    end
    
    if y == land(x)
        if v_x > 0
            beta = atan(v_y / v_x);
            elseif v_x < 0
            beta = pi + atan(v_y / v_x);
        elseif v_x == 0
            beta = pi / 2;
        end
    
        if x - x_prev >= 0
            x_next = x + 10^(-12);
            gamma = atan( (land(x_next) - land(x))/(x_next - x) );
        elseif x - x_prev < 0
            x_next = x + 10^(-12);
            gamma = pi + atan( (land(x_next) - land(x))/(x_next - x) );
        end

        beta = -beta + 2 * gamma;
    
        E = 5;
        v_sq = v_x^2 + v_y^2 - 2 * E / m;
        if v_sq <= 0
            break;
        end
        v_new = sqrt(v_sq);

        v_x = v_new * cos(beta);
        v_y = v_new * sin(beta);
    end
    
end


function vel_x = f_x(t, x, v_x, y, v_y)
    vel_x = v_x;
end

function vel_y = f_y(t, x, v_x, y, v_y)
    vel_y = v_y;
end

function w_x = g_x(t, x, v_x, y, v_y)
    w_x = 0 - 0.09 * sqrt(v_x^2 + v_y^2) * v_x;;
end

function w_y = g_y(t, x, v_x, y, v_y)
    g = 9.81;
    w_y = -g - 0.09 * sqrt(v_x^2 + v_y^2) * v_y;
end

function land = land(x)
    land = 0.05 * (x).^2;
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
