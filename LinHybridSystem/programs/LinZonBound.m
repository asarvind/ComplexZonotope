function [ retval ] = LinZonBound( T,acz )
retval = abs(T*[acz.V, acz.W]*[acz.s; (acz.u-acz.l)/2])+T*(acz.c+(acz.W*(acz.u+acz.l)/2));

end

