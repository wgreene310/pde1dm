function varargout=callVarargFunc(func, args)
numExpectedArgs=nargin(func);
if numExpectedArgs>length(args)
  fn = func2str(func);
  error(['Function \"' fn '\" is expecting %d arguments but is being called with' ...
    ' only %d arguments.'], numExpectedArgs, length(args));
end
[varargout{1:nargout}]=func(args{1:numExpectedArgs});
end