function println_to_file(path::AbstractString, func::Function)
    open(path, "w") do io
        redirect_stdout(io) do
            func()
        end
    end
end