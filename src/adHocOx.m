bigexp1 = [];
for i=1:length(AGORAModelWithExprArr1)
    if ~isempty(AGORAModelWithExprArr1{i})
	Expr1 = AGORAModelWithExprArr1{i};
        for j=1:length(Expr1.expression)
	    bigexp1(end+1) = Expr1.expression(j);
        end
    end
end
	
bigexp2 = [];
for i=1:length(AGORAModelWithExprArr2)
    if ~isempty(AGORAModelWithExprArr2{i})
	Expr2 = AGORAModelWithExprArr2{i};
        for j=1:length(Expr2.expression)
	    bigexp2(end+1) = Expr2.expression(j);
        end
    end
end

bigexp1 = bigexp1(bigexp1~=0);
bigexp2 = bigexp2(bigexp2~=0);
step1 = floor(length(bigexp1)/60)-1;
bigexp3 = [];
for i=1:step1
    bigexp3(i) = sum(bigexp1(i*60:(i+1)*60));
end
step2 = floor(length(bigexp2)/60)-1;
bigexp4 = [];
for i=1:step2
    bigexp4(i) = sum(bigexp2(i*60:(i+1)*60));
end

