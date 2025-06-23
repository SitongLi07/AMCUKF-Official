function [sigma, wm, wc] = OnSigmaPoint(x, P)
alph = 1.0;
beta = 2.0;
ka = 0;
L = size(P, 2);
lmd = (L+ka)*alph*alph - L;
PS = sqrtm((L+lmd)*P);
wm(1) = lmd/(L+lmd);
wc(1) = wm(1) + (1-alph*alph+beta);
sigma(1) = x;
for i=2:(2*L+1)
    wm(i) = 1.0/(2*(L+lmd));
    wc(i) = wm(i);
    if i <= L+1;
        sigma(i) = x + PS(:,i-1);
    else
        sigma(i) = x - PS(:,i-1-L);
    end
end;