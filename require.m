function require(bool_condition,message)
% require(bool_condition: Boolean, message: String): Null/Error
% Checks that a given condition is true, else gives an error with message

if ~bool_condition
    error(message)
end

end
    