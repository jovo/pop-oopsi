function z = tester1(x)

for i=1:numel(x)
    z=test1(x,i);
end

    function y=test1(x,i)
        y=exp(x(i));
    end
end