
function [h,rorig] = Plot_Perfomance_Profiles(A,logplot,colors,lines,markers)

if nargin < 3 || isempty(colors)
    colors  = {'blue' 'orange' 'black' 'red' 'cyan' 'magenta' 'green' 'yellow'};   lines   = {'-' '-.' '--'};
    markers = [ 's' 'o' '^' 'v' 'p' '<' 'x' 'h' '+' 'd' '*' '<' ];
end

%A = A([1; find(~strcmp(A(:,1),''))],:);
np = size(A,1)-1;
ns = size(A,2);
T = zeros(np,ns);


for s = 1:ns
    for p = 1:np
        val = A{p+1,s};
        if ischar(val)
            val = str2double(val);
        end
        if (isempty(val)) || (isnan(val)&&~strcmp(A{p+2,s},'')) ||val == 0
            T(p,s) = NaN;
        else
            T(p,s) = val;
        end
    end
end

% Other colors, lines, and markers are easily possible:

if (nargin < 2); logplot = 0; end
figure(1);
lgnd = A(1,:);
[h,rorig] = DrawProfilePlot(T,ns,np,logplot,'',lgnd,colors,lines,markers);
end





function [hl,rorig] = DrawProfilePlot(T,ns,np,logplot,strtype,lgnd,colors,lines,markers)

strtimeiter = 'time';
%     strtimeiter = '(MAPE)';

% Compute ratios and divide by smallest element in each row.
r = T./repmat(min(T,[],2),1,ns);
rorig = r;
% Replace all NaN's with twice the max_ratio and sort.
max_ratio = max(max(r));
if isnan(max_ratio)
    max_ratio = 1;
end
r(isnan(r)) = 2*max_ratio;
r = sort(r);
numero = 2;%numero = 2;
r = [r;numero*max_ratio*ones(1,ns)];

% Plot stair graphs with markers.
hl = zeros(ns,1);
hold off;
for s = 1:ns
    [xs,ys] = stairs(r(:,s),(1:np+1)/np);
    
    % Only plot one marker at the intercept
    if (xs(1)==1)
        vv = find(xs==1,1,'last');
        xs = xs(vv:end);   ys = ys(vv:end);
    end
    
    sl = mod(s-1,length(lines)) + 1; sc = mod(s-1,length(colors)) + 1; sm = mod(s-1,length(markers)) + 1;
    %option1 = [char(lines(sl)) colors(sc) markers(sm)];
    if (logplot)
        hl(s) = semilogx(xs,ys,'Line',char(lines(sl)),'Marker',markers(sm),'Color',rgb(char(colors(sc))),'MarkerSize',6);
    else
        %hl(s) = plot(xs,ys,'LineStyle',char(lines(sl)),'Marker',markers(sm),'Color',rgb(char(colors(sc))),'MarkerSize',6);
        %varlength = length(xs);
        %mylength=ceil(varlength/3);
        %xs = xs(1:end)   
        %ys = ys(1:end)
        %pause
        hl(s) = plot(xs,ys,'LineStyle',char(lines(sl)),'Color',rgb(char(colors(sc))));
        %hl(s) = plot(xs,ys,'Marker',markers(sm),'Color',rgb(char(colors(sc))),'MarkerSize',6);
    end
    hold on;
end

%title(['Performance Profiles (',num2str(np),' ',strrep(strtype,'WNN_',''),' problems) (',strtimeiter,')']);
title(['Performance Profiles (','100',' ',strrep(strtype,'WNN_',''),' problems) (',strtimeiter,')']);
SetLegend(lgnd);
% Axis properties are set so that failures are not shown, but with the
% max_ratio data points shown. This highlights the "flatline" effect.
if (logplot)
    axis([1 1.1*max_ratio 0 1.05]);
    twop = floor(log2(1.1*max_ratio));
    set(gca,'XTick',2.^(0:twop))
else
    axis([1 1.1*max_ratio 0 1.05]);
end

if ns > 2;
    pluralstr = 's';
else
    pluralstr = '';
end
xlabel(['at most  x \times (the ',strtimeiter,' of the other method',pluralstr,')']);
ylabel('% of problems');
end

function SetLegend(lgnd)
s = ['''',lgnd{1},''''];
for i = 2 : length(lgnd)
    s = [s,',''',lgnd{i},''''];
end
s = ['legend(',s,')'];
eval(s);
end