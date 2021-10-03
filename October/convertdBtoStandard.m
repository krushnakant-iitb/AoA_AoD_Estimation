function val= convertdBtoStandard(valtoconvert)
% Converts values to linear from dB
    val= 10.^(valtoconvert/10);
end