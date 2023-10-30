"""
Example for scriptflow 

Adapt hpc parameters to your specifications. Then run via the command:
> scriptflow run mysim
"""

import scriptflow as sf
import os

# set executor to slurm or hpc 

# sf.init({ # Runner for Slurm
#     "executors":{
#         "slurm":{
#             "maxsize": 3,
#             "account": 'pi-chansen1',
#             "user": 'wiemann',
#             "partition": 'standard',
#             "modules": 'julia/1.8',
#             "walltime": '00:01:00'
#         } 
#     },
#     'debug': True
# })

sf.init({ # Runner for PBS
    "executors":{
        "hpc":{
            "account": "abram",
            "maxsize": 20,
            "modules": 'julia/1.8.3',
            "walltime": '00:15:00'
        } 
    },
    'debug': True
})

# sf.init({ # local runner
#     "executors":{
#         "local": {
#             "account" : "abram",
#             "maxsize" : 5
#         } 
#     },
#     'debug':True
# })

# create temp-directory to store results in
temp_dir = 'temp'
if not os.path.exists(temp_dir):
    os.mkdir(temp_dir)

# define a flow called Rit
async def flow_mysim():

    # Initialize output
    task_init = [
        sf.Task(
        cmd = f"julia init_output.jl",
        # outputs = f"{temp_dir}/res_{i}.csv",
        name = f"initialize_output")
        # for i in range(5)
        # for i in [2,3,4,5]
    ]

    await sf.bag(*task_init)

    # Calibrate each version of model
    tasks = [
        sf.Task(
        cmd = f"julia calibrate_model.jl {i}",
        # outputs = f"{temp_dir}/res_{i}.csv",
        name = f"sim-cal-{i}").set_cpu(2).set_memory(4)
        # for i in range(5)
        # for i in range(2, 113, 1)
        for i in range(2, 11, 1)
    ]

    await sf.bag(*tasks)

    # Decompose each change in model
    tasks = [
        sf.Task(
        cmd = f"julia decomposition.jl {i}",
        # outputs = f"{temp_dir}/res_{i}.csv",
        name = f"sim-decomp-{i}").set_cpu(2).set_memory(4)
        # for i in range(5)
        # for i in range(2, 112, 2)
        for i in range(2, 10, 2)
    ]

    await sf.bag(*tasks)

    # Aggregates the simulation results and stores as .csv
    # t_agg = sf.Task(
    #     cmd = f"julia agg_results.jl {temp_dir}",
    #     outputs = "results.csv",
    #     name = "agg-results")
    
    # await t_agg

    # Combine output
    task_combine = [
        sf.Task(
        cmd = f"julia combine_output.jl",
        # outputs = f"{temp_dir}/res_{i}.csv",
        name = f"combine_output")
        # for i in range(5)
        # for i in [2,3,4,5]
    ]

    await sf.bag(*task_combine)


    # Can add cleanup-tasks here (e.g., to remove temp_dir)
    # for f in os.listdir(temp_dir):
    #     os.remove(os.path.join(temp_dir, f))