function [N,r,dr] = dicretization()
    N = 50*1;
    r_node = linspace(0,1,N+1);
    r = zeros(N,1);
    dr = 1/N;
    for i = 1:N
        r (i) = ( r_node(i)+r_node(i+1) )/2;
    end
end