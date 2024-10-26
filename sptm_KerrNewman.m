classdef sptm_KerrNewman
    properties
        M = 1;
        a = 0.999
        Q = 0;
    end

    methods
        function guv = g_uv(obj,x)
            t = x(1); r = x(2); th = x(3); ph = x(4);
            
            Del = r^2 + obj.a^2 - 2*obj.M*r + obj.Q^2;
            Sig = r^2 + ((obj.a^2) * (cos(th)^2));
            
            gtt = -( 1 - (2*obj.M*r - obj.Q^2)/Sig );
            gtph = -(obj.a*(2*obj.M*r - obj.Q^2)*(sin(th)^2))/Sig;
            gphph = (r^2 + obj.a^2 + ((2*obj.M*r - obj.Q^2)*(obj.a^2)*(sin(th)^2)/Sig) )*sin(th)^2;
            
            guv =   [   [gtt,   0,          0,      gtph    ];...
                        [0,     Sig/Del,    0,      0       ];...
                        [0,     0,          Sig,    0       ];...
                        [gtph,  0,          0,      gphph   ]   ];
        end

        function comp = g_uv_comp(obj,x,i,j)
            guv = obj.g_uv(x);
            comp = guv(i,j);
        end

    end
end
        