q = 1;
m = 1000;
A=zeros(m,m);
n = 10;
for i=1:m
    for j=1:m
        for k=1:n                
            r = [pi*i/m pi*j/m];
            Tr = C1(q,k,pi/3)+C1(q,k,r(1))-2*pi*V(r(2),q);
            Det = (C1(q,k,pi/3)-pi*V(r(2),q))*(C1(q,k,r(1))-pi*V(r(2),q))-pi^2*((V(r(2),q-k)+V(r(2),q+k))^2)/4;
            if [Tr>0 || Det<0] && A(i,j)==0
                A(i,j) = 1;
            end
        end
    end
end
A = A';
figure('pos',[100 100 1000 800]);
ax = subplot(1,1,1);
imagesc([0,pi],[0,pi],A);
xlabel("r_{2,2}");
ylabel("r_{1,2}");
ax.YDir = "normal";
title('Largest c1 greater than zero?')
C = colormap();
L = size(C,1);
Gs = round(interp1(linspace(min(A(:)),max(A(:)),L),1:L,A));
H = reshape(C(As,:),[size(As) 3]);
subplot(1,2,2);

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
