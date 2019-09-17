function z=zerosLike(varargin)
persistent isOctave
if isempty(isOctave)
  isOctave = exist('OCTAVE_VERSION', 'builtin');
end
x=varargin{end};
%z = builtin('repmat', x(1), m, n);
if(isOctave)
  z = x(1)*zeros(varargin{1:end-1});
else
  z = zeros(varargin{1:end-1}, 'like', x);
end
end