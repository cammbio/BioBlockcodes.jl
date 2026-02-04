using GCATCodes
using BioSequences: LongDNA
using Base.Threads: Atomic, nthreads, @spawn
using BenchmarkTools

function benchmark_channel_size(;
    channel_size = 1024,
    items = 2000,
    workers = Threads.nthreads(),
    work_sleep = 5e-5,
)
    ch = Channel{Int}(channel_size)
    prod_blocks = Threads.Atomic{Int}(0)

    prod_time = 0.0
    total_time = @elapsed begin
        @sync begin
            prod_task = @async begin
                t0 = time()
                for i in 1:items
                    t_put = time()
                    put!(ch, i)
                    if time() - t_put > 0
                        prod_blocks[] += 1
                    end
                end
                close(ch)
                return time() - t0
            end

            for _ in 1:workers
                @spawn begin
                    for _ in ch
                        sleep(work_sleep)
                    end
                end
            end

            prod_time = fetch(prod_task)
        end
    end

    return (
        channel_size = channel_size,
        items = items,
        workers = workers,
        work_sleep = work_sleep,
        producer_time = prod_time,
        total_time = total_time,
        prod_block_events = prod_blocks[],
    )
end


function sweep_channel_sizes(
    sizes;
    items::Int = 2000,
    workers::Int = Threads.nthreads(),
    work_sleep::Float64 = 5e-5,
)
    results = [
        benchmark_channel_size(channel_size = s, items = items, workers = workers, work_sleep = work_sleep) for s in sizes
    ]
    println("size\tprod_s\ttotal_s\tblocks")
    for r in results
        println(
            "$(r.channel_size)\t$(round(r.producer_time, digits = 4))\t$(round(r.total_time, digits = 4))\t$(get(r, :prod_block_events, 0))",
        )
    end
    return results
end

# sweep_channel_sizes([0, 64, 256, 1024, 4096]; items = 5000, workers = Threads.nthreads(), work_sleep = 5e-4)

println("Done loading script.")