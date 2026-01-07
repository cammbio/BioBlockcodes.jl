function println_to_file(dst_path::AbstractString, func::Function)
    open(dst_path, "w") do input
        redirect_stdout(input) do
            func()
        end
    end
end