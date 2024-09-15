% What happens here? We follow different eigenvalues along a straight line of the parameter r and plot their values.
q = [1 1]; % twist
k1 = [0 1];
k2 = [1 1];
k3 = [2 2];
for i=1:1000
    x = i*pi/1000;
    r = [x/2 x]; % the line r follows 
    y1(i) = C1(q,k1,r);
    y2(i) = C1(q,k2,r);
    y3(i) = C1(q,k3,r);
end
x=(pi/1000):(pi/1000):pi;
figure(1)
plot(x,y1,"r-","Linewidth",2);
hold on;
plot(x,y2,"b-","Linewidth",2);
hold on;
plot(x,y3,"b-","Linewidth",2,"Color",[0 0.7 0]);
grid on;
legend('(0,1)','(1,1)','(2,2)');

function c1 = C1(q,k,r)
    c1 = (2*pi)^length(r)*[V(r,q+k)+V(r,q-k)-2*V(r,q)]/4;
end

function v = V(r,k)
% V Calculates the fourier coefficient for the given values
    v = 2;
    for i=1:length(k)
        if k(i) == 0
            v_i = r(i)/pi;
        else
            v_i = sin(k(i)*r(i))/(pi*k(i));
        end
        v = v*v_i;
    end
end
