function [clf,cdf,cmf] = fn_clcdcm(alpha,data)

alpha = alpha*180/pi;
% load data.dat;
alpha1 = data(:,1);
cl = data(:,2);
cd = data(:,3);
cm = data(:,5);
alpha2 = alpha1;
% load data2.dat;
% alpha2 = data2(:,1);
% cd = data2(:,3);

len1 = length(alpha1);
len2 = length(alpha2);

%%
for i = 1:len1-1
    if (alpha >= alpha1(i))&&(alpha <= alpha1(i+1))        
        s = (alpha-alpha1(i))/( alpha1(i+1)-alpha1(i) );
        h1 = 1-s;
        h2 = s;
        clf = h1*cl(i)+h2*cl(i+1);        
    end
end

if (alpha > alpha1(end) )
    clf = cl(end)*cos((alpha-alpha1(i))*pi/180*0);
end
if (alpha < alpha1(1) )
    clf = cl(1)*cos((alpha-alpha1(i))*pi/180*0);
end
%%
for i = 1:len2-1
    if (alpha >= alpha2(i))&&(alpha <= alpha2(i+1))        
        s = (alpha-alpha2(i))/( alpha2(i+1)-alpha2(i) );
        h1 = 1-s;
        h2 = s;
        cdf = h1*cd(i)+h2*cd(i+1);        
    end
end

if (alpha > alpha2(end) )
    cdf = cd(end);
	p = [4.99452440667644e-08,2.74382429769118e-08,-1.54473129502502e-05,0.000105908923449461,0.000449422569440989,0.00721057458308516]';
    n=length(p)-1;    
    sum =0;
    k =n;
    for j=1:n+1
        sum = sum + p(j)*abs(alpha)^k;
        k=k-1;
    end
    cdf = sum;
%     cdf = cd(end);
end
if (alpha < alpha2(1) )
    cdf = cd(1);
	p = [4.99452440667644e-08,2.74382429769118e-08,-1.54473129502502e-05,0.000105908923449461,0.000449422569440989,0.00721057458308516]';
    n=length(p)-1;    
    sum =0;
    k =n;
    for j=1:n+1
        sum = sum + p(j)*abs(alpha)^k;
        k=k-1;
    end
    cdf = sum;
end

%%
for i = 1:len2-1
    if (alpha >= alpha2(i))&&(alpha <= alpha2(i+1))        
        s = (alpha-alpha2(i))/( alpha2(i+1)-alpha2(i) );
        h1 = 1-s;
        h2 = s;
        cmf = h1*cm(i)+h2*cm(i+1);        
    end
end

if (alpha > alpha2(end) )
    cmf = cm(end);
end
if (alpha < alpha2(1) )
    cmf = cm(1);
end
    

end