function z = tester2(x)

for i=1:numel(x)
    z=test2(x(i));
end

    function y=test2(x)
        y=exp(x);
    end

end