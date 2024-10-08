% type colormap(cmap) to get the colors of the thesis

q = [1 1]; % twist
m = 1000; % resolution
Stable = zeros(m,m); % is zero in the stable zone, otherwise 1
n = 10;
C1max = zeros(m,m)-10; % stores the maximal eigenvalue
BiggestFunction = zeros(m,m); % stores the corresponding eigenfunction with the value 10*abs(k1)+abs(k2)
for i = 1:m
    for j = 1:m
        for k1 = 0:n
            for k2 = -n:n
                k = [k1 k2];
                r = [pi*i/m pi*j/m];
                z = C1(q,k,r);
                if k1 ~= 0 || k2 ~= 0
                    if C1max(i,j) < z
                        C1max(i,j) = z;
                        BiggestFunction(i,j) = 10*abs(k1)+abs(k2);
                    end
                end
                if z > 0 && Stable(i,j) == 0
                    Stable(i,j) = 1;
                end
            end
        end
    end
end

% post-processing
for i = 1:m
    for j = 1:m
        if Stable(i,j) ~= 0
            C1max(i,j) = 99;
            BiggestFunction(i,j) = 0;
        end
    end
end
for i = 1:m
    for j = 1:m
        if BiggestFunction(i,j) == 0
            BiggestFunction(i,j) = 0;
        elseif BiggestFunction(i,j) == 1
            BiggestFunction(i,j) = 1;
        elseif BiggestFunction(i,j) == 10
            BiggestFunction(i,j) = 2;
        elseif BiggestFunction(i,j) == 11
            BiggestFunction(i,j) = 3;
        elseif BiggestFunction(i,j) == 12
            BiggestFunction(i,j) = 4;
        elseif BiggestFunction(i,j) == 21
            BiggestFunction(i,j) = 5;
        end
    end
end
Stable = Stable';
B = BiggestFunction';
cmap = zeros(10, 6);
cmap = [1, 0, 0; ...   % Red for 0
  0, 0.8, 0; ...       % lightgreen 1
  0, 1, 0; ...       % darkgreen for 2
  1, 1, 0; ...       % yellow for 3
  0, 0, 0.1; ...       % light blue for 4
  0, 0, 1];       % dark blue for 5
colormap(cmap);
figure('pos',[100 100 1000 800]);
ax = subplot(1,1,1);
imagesc([0,pi],[0,pi],B);
xlabel("r_1");
ylabel("r_2");
ax.YDir = "normal";
title('Largest eigenfunctions?');
C = colormap();
C = cmap();
L = size(C,1);
Gs = round(interp1(linspace(min(B(:)),max(B(:)),L),1:L,B));
H = reshape(C(Bs,:),[size(Bs) 3]);
subplot(1,2,2);
imshow(B, cmap);

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
