function print_to_file(dst_path::AbstractString, fnc::Function, args...)
    open(dst_path, "w") do input
        redirect_stdout(input) do
            fnc(args...)
        end
    end
end