function [ retval1,retval2 ] = init_maxvect( v )
bool =  (v>-inf);
retval1 = bool;
retval2 = approx(bool,v);
end

