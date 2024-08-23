# What is Multi-Threading

Multi-threading is a way to execute a single task across multiple CPU's. Each CPU typically has a thread, which is an instrution for the task that for a specific CPU. 

Each CPU usually has a single thread. To check how many threads are availible to you, open the Julia REPL and enter,

```julia
julia> versioninfo()
Julia Version 1.9.3
Commit bed2cd540a (2023-08-24 14:43 UTC)
Build Info:

    Note: This is an unofficial build, please report bugs to the project
    responsible for this build and not to the Julia project unless you can
    reproduce the issue using official builds available at https://julialang.org/downloads

Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: 32 × AMD EPYC 7F32 8-Core Processor
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-14.0.6 (ORCJIT, znver2)
  Threads: 1 on 32 virtual cores
Environment:
  LD_LIBRARY_PATH = /share/spack/linux-rhel7-x86_64/gcc-10.3.0/pcre2-10.42-ulnwntf53aqpkc53zfsx5ioexq7cab2d/lib
  LD_LIBRARY_PATH_modshare = /share/spack/linux-rhel7-x86_64/gcc-10.3.0/pcre2-10.42-ulnwntf53aqpkc53zfsx5ioexq7cab2d/lib:1)
```
Here, the line `Threads: 1 on 32 virtual cores` tells me I'm currently using 1 of 32 possible threads availble to me. This is simimilar to me having 32 laptops on a desk and I'm just sitting at one with Julia open.

Another way to verify this is by typing

```julia
julia> Threads.nthreads()
1
```

By default, Julia starts up with a single thread of execution. If you want to us more threads, which would be analagous to more people sitting down at each one of your CPU's, you need to get out of the current Julia REPL session back to your shell and use

```bash
export JULIA_NUM_THREADS=32
```

(The above command works on bourne shells on Linux and OSX. Note that if you're using a C shell on these platforms, you should use the keyword `set` instead of `export`. If you're on Windows, start up the command line in the location of `julia.exe` and use `set` instead of `export`.)

Now let's start a new Julia REPL session and see if we're using all 32 threads.

```julia
julia> Threads.nthreads()
32
```

Great, again, the analogy is now we have 32 laptops (CPU's) with a person (thread) at each one.

Of all the people sitting at each laptop, we have to have a leader. To find who the leader, called a "master", we type

```julia
julia> Threads.threadid()
1
```
So the master thread has the ID `1`.

## Single Thread Example

So we have 32 laptops with a person sitting at each one. Let's say we just want on person, the master, to execute a task 32 times. We'll define that task in the loop below.

```julia
julia> 
for task_number in 1:Threads.nthreads()
           println("Task Number: ", task_number, "\t Done by Thread ID: ", Threads.threadid())
 end
Task Number: 1   Done by Thread ID: 1
Task Number: 2   Done by Thread ID: 1
Task Number: 3   Done by Thread ID: 1
Task Number: 4   Done by Thread ID: 1
Task Number: 5   Done by Thread ID: 1
Task Number: 6   Done by Thread ID: 1
Task Number: 7   Done by Thread ID: 1
Task Number: 8   Done by Thread ID: 1
Task Number: 9   Done by Thread ID: 1
Task Number: 10  Done by Thread ID: 1
Task Number: 11  Done by Thread ID: 1
Task Number: 12  Done by Thread ID: 1
Task Number: 13  Done by Thread ID: 1
Task Number: 14  Done by Thread ID: 1
Task Number: 15  Done by Thread ID: 1
Task Number: 16  Done by Thread ID: 1
Task Number: 17  Done by Thread ID: 1
Task Number: 18  Done by Thread ID: 1
Task Number: 19  Done by Thread ID: 1
Task Number: 20  Done by Thread ID: 1
Task Number: 21  Done by Thread ID: 1
Task Number: 22  Done by Thread ID: 1
Task Number: 23  Done by Thread ID: 1
Task Number: 24  Done by Thread ID: 1
Task Number: 25  Done by Thread ID: 1
Task Number: 26  Done by Thread ID: 1
Task Number: 27  Done by Thread ID: 1
Task Number: 28  Done by Thread ID: 1
Task Number: 29  Done by Thread ID: 1
Task Number: 30  Done by Thread ID: 1
Task Number: 31  Done by Thread ID: 1
Task Number: 32  Done by Thread ID: 1
```

As you can see, the master `Thread ID: 1` when started with `Task Number: 1` and printed the line and moved on to the next. Nothing special. But we have 31 other threads just sitting around doing nothing while `Thread ID: 1` did all the work.

Let's change it so that each of the 32 `task_number`'s are executed by each thread at the same time. Changing the script is easy. We simply use the macro `@threads` to brodcast the task across all threads availible.

```julia
Threads.@threads for task_number in 1:Threads.nthreads()
           println("Task Number: ", task_number, "\t Done by Thread ID: ", Threads.threadid())
       end
```
What we get is,

```julia
Task Number: 1   Done by Thread ID: 2
Task Number: 32  Done by Thread ID: 26
Task Number: 30  Done by Thread ID: 26
Task Number: 31  Done by Thread ID: 26
Task Number: 22  Done by Thread ID: 11
Task Number: 3   Done by Thread ID: 3
Task Number: 28  Done by Thread ID: 32
Task Number: 14  Done by Thread ID: 18
Task Number: 27  Done by Thread ID: 31
Task Number: 13  Done by Thread ID: 15
Task Number: 12  Done by Thread ID: 17
Task Number: 2   Done by Thread ID: 1
Task Number: 17  Done by Thread ID: 21
Task Number: 4   Done by Thread ID: 6
Task Number: 24  Done by Thread ID: 28
Task Number: 26  Done by Thread ID: 30
Task Number: 19  Done by Thread ID: 13
Task Number: 7   Done by Thread ID: 5
Task Number: 15  Done by Thread ID: 19
Task Number: 11  Done by Thread ID: 16
Task Number: 25  Done by Thread ID: 29
Task Number: 23  Done by Thread ID: 10
Task Number: 6   Done by Thread ID: 9
Task Number: 5   Done by Thread ID: 4
Task Number: 10  Done by Thread ID: 14
Task Number: 21  Done by Thread ID: 12
Task Number: 20  Done by Thread ID: 24
Task Number: 9   Done by Thread ID: 8
Task Number: 18  Done by Thread ID: 22
Task Number: 29  Done by Thread ID: 27
Task Number: 8   Done by Thread ID: 7
Task Number: 16  Done by Thread ID: 20
```
Each loop for the for-loop was executed by a different thread - each thread only did a single task. It was executed in *parallel*.

## Data Races

Because we have multiple threads executing tasks at the same time, it's possible that one thread is reading values stored in a memory location, while another thread is writing to it. This is called a *data race*.

Here's an example. Let's create a vector,

```julia
julia> n = 1000
1000

julia> vector  = collect(1:n)
5-element Vector{Int64}:
 1
 2
⋮
 999
 1000
 ```
And we know the sum of the values in each index of that vector to be

```julia
julia> sum(vector)
500500
```
Let's do that same caculation using loops similar to before. With just one thread, we define a function that starts with `0` as the and takes the value at each index and adds it to `sum_out`

```julia
function oneThreadSum(vector)
    sum_out = 0
    for index in eachindex(vector)
        sum_out += vector[index]
    end
    return sum_out
end
```
This returns,

```julia
julia> oneThreadSum(vector)
500500
```

Great. But now let's do it with multiple threads. Just as we did before, we broadcast the execution of the loop with `Threads.@threads`

```julia
function multThreadSumWrong(vector)
    sum_out = 0
    Threads.@threads for index in eachindex(vector)
        sum_out += vector[index]
    end
    return sum_out
end
```
which returns,

```julia
julia> multThreadSumWrong(vector)
393529
```
That's not right. If we keep doing it,

```julia
julia> multThreadSumWrong(vector)
393529

julia> multThreadSumWrong(vector)
439730

julia> multThreadSumWrong(vector)
500500

julia> multThreadSumWrong(vector)
402128
```
It get gives a different answer everytime. Even though it got it right once, it's still inconsistant. The more loops it has to do, the more inconsistant it becomes. This is a data race.

The reason why it happens is because all threads have access to the shared memory `sum_out`. Each thread is reading `sum_out` at the same time another time another one is writing to it. So it is constantly getting overwritten.

## How to Prevent Data Races

Instead of doing a single shared memory location to accumulate the values, we set up a separate memory location `sum_out` for each thread `Threads.nthreads()` to store is values with

```julia
sum_out = zeros(Int, Threads.nthreads())
```
and when all threads are done with their tasks we add up the values in each respective memory location "owned" by each thread with `sum(sum_out)`. The function is as follows,

```julia
function multThreadSumCorrect(vector)
    sum_out = zeros(Int, Threads.nthreads())
    Threads.@threads for index in eachindex(vector)
        sum_out[Threads.threadid()] += vector[index]
    end
    return sum(sum_out)
end
```
The result is

```julia
julia> multThreadSumCorrect(vector)
500500

julia> multThreadSumCorrect(vector)
500500

julia> multThreadSumCorrect(vector)
500500
```

## Going Forward

We can use this as a model for other jobs, such as compiling data or running simulations. 