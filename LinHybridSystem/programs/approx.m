function [ retval ] = approx( bool,v )
    if isrow(v)
        v = v';
    end
  retval = v;
  for i = 1:numel(v)
    if bool(i) == 0
        retval(i,1) = 0;
    end
  end
end

