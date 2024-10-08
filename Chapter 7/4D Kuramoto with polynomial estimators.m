n = 4; % dimension
c = 0; % counts how often the indicator in line 29 is positive
d = 0; % counts how often the indicator is positive, but the configuration still returns to  [0 0 ... 0]
splaystate = 0; % counts the number of splaystates for dimension 4 and in general the number of times, where
% the configuration doesn't end up in [0 0 ... 0]
T = 30000;
dt = 0.005;
for Stat = 1:10000
    x = 2*pi*(rand(n-1,1)-0.5);
    x = [x; -sum(x)]; % To ensure that the configuration can converge towards [0 0 ... 0]
    y = x;
    for j = 1:T
         y = KurSchritt(y,dt);
    end
    if max(y)>0.1
        splaystate = splaystate+1;
    end
    a1 = KurSchritt(x,dt);
    a2 = KurSchritt(a1,dt);
    a3 = KurSchritt(a2,dt);
    d1 = a1-x; % Approximation of first derivative
    d2 = a2+x-2*a1; % Approximation of second derivative
    p1 = d1'*x; % First Indicator
    p2 = d2'*x-d1'*d1; % Second Indicator
    p3 = 0; % Third Indicator
    for m = 1:n
        p3 = p3+x(m)*((x(m)*x(m))*(a3(m)-3*a2(m)+3*a1(m)-x(m))-5*x(m)*((a1(m)-x(m))*(a2(m)+x(m)-2*a1(m)))+4*(a1(m)-x(m))*((a1(m)-x(m))*(a1(m)-x(m))));
    end
    if p1>0 || p2>0 || p3>0 % Set here your custom mix of Indicators
        c = c+1;
        if max(y)<0.1
            d = d+1;
        end
        disp("Splaystate probability: "+splaystate/Stat);
        disp("Indicator positive: "+c/Stat);
        disp("Failed: "+d/c);
    end
end

function y = KurSchritt(x,dt)
y = x;
n = length(x);
    if n>2
        for i = 2:n-1
            y(i) = x(i)+dt*(sin(x(i-1)-x(i))+sin(x(i+1)-x(i)));
        end
        y(1) = x(1)+dt*(sin(x(n)-x(1))+sin(x(2)-x(1)));
        y(n) = x(n)+dt*(sin(x(n-1)-x(n))+sin(x(1)-x(n)));
    elseif n == 2
        a = x(1)+2*dt*sin(x(2)-x(1));
        b = x(2)+2*dt*sin(x(1)-x(2));
        y = [a; b];
    end
end
