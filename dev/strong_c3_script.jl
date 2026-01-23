using GCATCodes
using Base.Threads
using Logging

# activate logging
global_logger(ConsoleLogger(Logging.Debug))

# cancel flag
const stop_flag = Atomic{Bool}(false)
# track status per combination_size
const job_states = Dict{Int, Symbol}()  # :queued, :running, :finished, :failed
# track thread assignment per combination_size
const job_thread_ids = Dict{Int, Int}()
const job_thread_lock = SpinLock()


function run_jobs_ordered(;
    min_combination_size = 1,
    max_combination_size = 20,
    nworkers = Threads.nthreads(),
)
    stop_flag[] = false

    # reset job states and thread ids
    empty!(job_states)
    empty!(job_thread_ids)

    sizes = collect(min_combination_size:max_combination_size)
    order_channel = Channel{Int}(length(sizes))
    for combination_size in sizes
        put!(order_channel, combination_size)
        job_states[combination_size] = :queued
    end
    close(order_channel)

    jobs = Dict{Int, Task}()
    for worker in 1:min(nworkers, length(sizes))
        jobs[worker] = @spawn begin
            for combination_size in order_channel
                lock(job_thread_lock)
                job_thread_ids[combination_size] = Threads.threadid()
                job_states[combination_size] = :running
                unlock(job_thread_lock)
                try
                    process_strong_c3_combinations_by_combination_size(
                        GCATCodes.ALL_CODONS,
                        combination_size,
                        "files/results/test$(combination_size).txt",
                        "files/checkpoints/test$(combination_size)_cp.txt",
                        stop_flag;
                        show_debug = true,
                    )
                    lock(job_thread_lock)
                    job_states[combination_size] = :finished
                    unlock(job_thread_lock)
                catch err
                    lock(job_thread_lock)
                    job_states[combination_size] = :failed
                    unlock(job_thread_lock)
                    @warn "Combination $(combination_size) failed: $(sprint(showerror, err))"
                end
            end
        end
    end

    println("Jobs scheduled (ordered channel) from $min_combination_size to $max_combination_size.")
    return jobs
end


# schedule jobs
# function run_jobs(; min_combination_size::Int = 1, max_combination_size::Int = 20)
#     # always start with a cleared cancel flag
#     stop_flag[] = false

#     empty!(job_states)
#     empty!(job_thread_ids)
#     jobs = Dict{Int, Task}()
#     for combination_size in min_combination_size:max_combination_size
#         jobs[combination_size] = @spawn begin
#             thread_id = Threads.threadid()
#             lock(job_thread_lock)
#             job_thread_ids[combination_size] = thread_id
#             job_states[combination_size] = :running
#             unlock(job_thread_lock)
#             try
#                 process_strong_c3_combinations_by_combination_size(
#                     GCATCodes.ALL_CODONS,
#                     combination_size,
#                     "files/results/test$(combination_size).txt",
#                     "files/checkpoints/test$(combination_size)_cp.txt",
#                     stop_flag;
#                     show_debug = true,
#                 )
#                 lock(job_thread_lock)
#                 job_states[combination_size] = :finished
#                 unlock(job_thread_lock)
#             catch err
#                 lock(job_thread_lock)
#                 job_states[combination_size] = :failed
#                 unlock(job_thread_lock)
#                 @warn "Combination $(combination_size) failed: $(sprint(showerror, err))"
#             end
#         end
#     end

#     println("Jobs scheduled from combination_size $min_combination_size to $max_combination_size.")
#     return jobs
# end


# cancel all tasks
function cancel_tasks(jobs::Dict{Int, Task})
    isempty(jobs) && throw(ArgumentError("No jobs provided."))

    stop_flag[] = true

    for (combination_size, task) in jobs
        try
            fetch(task)
        catch err
            err isa InterruptException && @info "Task $combination_size cancelled."
            @warn "Task $combination_size encountered an error: $(sprint(showerror, err))"
        end
    end

    stop_flag[] = false
end


# collect error messages from jobs dict
function get_error_messages(jobs::Dict{Int, Task})
    isempty(jobs) && throw(ArgumentError("No jobs provided."))

    error_messages = String[]
    for (combination_size, task) in jobs
        try
            fetch(task)
        catch err
            push!(error_messages, "combination_size=$(combination_size): $(sprint(showerror, err))")
        end
    end
    return error_messages
end


# get status for a single task bound to combination_size
function get_status(jobs::Dict{Int, Task}, combination_size::Int)
    isempty(jobs) && throw(ArgumentError("No jobs provided."))
    haskey(job_states, combination_size) ||
        haskey(jobs, combination_size) ||
        return "Unavailable combination_size $(combination_size)"

    # get task (if present) and thread id snapshot
    task = haskey(jobs, combination_size) ? jobs[combination_size] : nothing
    lock(job_thread_lock)
    thread_id = get(job_thread_ids, combination_size, missing)
    state_symbol = get(job_states, combination_size, nothing)
    unlock(job_thread_lock)

    # derive status label and flags
    status_label = "Running"
    started_flag = false
    failed_flag = false
    done_flag = false

    if state_symbol !== nothing
        if state_symbol === :failed
            status_label = "Failed"
            started_flag = true
            failed_flag = true
            done_flag = true
        elseif state_symbol === :finished
            status_label = "Finished"
            started_flag = true
            done_flag = true
        elseif state_symbol === :queued
            status_label = "Waiting"
            started_flag = false
        else
            status_label = "Running"
            started_flag = true
        end
    end

    if task !== nothing
        started_flag = istaskstarted(task)
        failed_flag = istaskfailed(task)
        done_flag = istaskdone(task)
        status_label = if done_flag && failed_flag
            "Failed"
        elseif done_flag
            "Finished"
        elseif task.state === :waiting
            "Waiting"
        else
            "Running"
        end
    end

    # collect status info
    lines = String[]
    push!(lines, "combination_size=$(combination_size)")
    push!(lines, "status=$(status_label)")
    push!(lines, "thread_id=$(thread_id === missing ? "-" : string(thread_id))")
    push!(lines, "started=$(started_flag)")
    push!(lines, "failed=$(failed_flag)")
    push!(lines, "done=$(done_flag)")

    return join(lines, "\n")
end


# get all statuses for all tasks bound to combination_size
function get_status_all(jobs::Dict{Int, Task})
    isempty(jobs) && throw(ArgumentError("No jobs provided."))

    formatted = IOBuffer()
    combos = !isempty(job_states) ? sort(collect(keys(job_states))) : sort(collect(keys(jobs)))
    for combination_size in combos
        println(formatted, get_status(jobs, combination_size))
        println(formatted, "--------------------")
    end
    return String(take!(formatted))
end


println("Script loaded. You can now call run_jobs() to start processing strong C3 combinations.")
