function pde = KSdata2d

pde = struct('u0',@u0,'v0',@v0);

    function s = u0(p)
        x = p(:, 1);y = p(:, 2);
        s = 1000*exp( -100*( (x - 0.5).^2 + (y - 0.5).^2 ) );
        % s = 4*exp(-(x.^2+y.^2));
        % s = 100*exp( -( x.^2 + y.^2 ) );
    end


    function s = v0(p)
        x = p(:, 1);y=p(:, 2);
        s = 500*exp( -50*( (x - 0.5).^2 + (y - 0.5).^2) );
        % s = exp(-(x.^2+y.^2)/2);
        % s = 0*exp( -( x.^2 + y.^2 )/2 );
    end
end