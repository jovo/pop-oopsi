function z = tester(x)

for i=1:numel(x)
    y2=test2(x(i));
    y1=test1(x,i);
end
z=[y1 y2];

    function y=test1(x,i)
        y=exp(x(i));
    end

    function y=test2(x)
        y=exp(x);
    end

end