function mustBePositiveOrNan(x)
% mustBePositiveOrNan throw error if not positive or nan
% intended for use as an argument validation function
if ~(all(x>0 | isnan(x)))
  error("mustBePositiveOrNan:mustBePositiveOrNan", "Input must be positive or nan.");
end
end