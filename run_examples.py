"""
Run the target in the build folder to execute the generated code examples.

1. Navigate to the build folder and build the project.
2. Run all the executables with an optional argument for simulation time steps (default is 100):
   - python run_examples.py [t_sim]

This script will run all compiled executables for examples: toyExample, robotArm, helicopter, ...
"""

import subprocess
import sys
import os
import re
from pathlib import Path

def get_binary_size(binary_path):
    """Get the size of the binary file in bytes."""
    if not Path(binary_path).exists():
        print("Binary not found: {}".format(binary_path))
        return None
    return Path(binary_path).stat().st_size

def parse_results(output):
    """Parse total runtime and average runtime from executable output."""
    total_match = re.search(r'Total runtime over \d+ steps: ([\d.]+) ms', output)
    average_match = re.search(r'Average runtime per step: ([\d.]+) ms', output)
    cost_match = re.search(r'Cumulative cost: ([\d.eE+-]+)', output)

    total = float(total_match.group(1)) if total_match else None
    average = float(average_match.group(1)) if average_match else None
    cost = float(cost_match.group(1)) if cost_match else None

    return total, average, cost

def run_executable(exe_path, t_sim):
    """Run the executable and capture output and runtime."""
    try:
        result = subprocess.run([exe_path, str(t_sim)], capture_output=True, text=True, timeout=1200)
        return result.stdout + result.stderr
    except subprocess.TimeoutExpired:
        print(f"Timeout for {exe_path}")
        return None
    except FileNotFoundError:
        print(f"Executable not found: {exe_path}")
        return None

def main():
    # Default simulation steps
    t_sim = 1000
    if len(sys.argv) > 1:
        try:
            t_sim = int(sys.argv[1])
        except ValueError:
            print("Invalid simulation steps argument")
            sys.exit(1)

    # Assume build directory is at generated_code/build
    build_dir = Path(__file__).parent / "build"
    if not build_dir.exists():
        print(f"Build directory {build_dir} does not exist. Please build the project first.")
        sys.exit(1)

    examples = ["toyExample", "robotArm", "helicopter", "springMass", "spacecraft", "quadcopter"]
    solvers = ["ALADIN", "osqp", "qpOASES", "hpipm"]

    # Dictionary to map solver names to executable prefixes
    solver_prefixes = {
        "ALADIN": "ALADIN_",
        "osqp": "osqp_",
        "qpOASES": "qpOASES_",
        "hpipm": "hpipm_",
    }


    res_dict = {}
    for example in examples:
        for solver in solvers:
            print("Running example '{}' with solver '{}'...".format(example, solver))
            prefix = solver_prefixes[solver]
            exe_name = f"{prefix}{example}"
            exe_path = build_dir / solver / exe_name

            # Get binary size
            size = get_binary_size(exe_path)
            total_runtime = None
            avg_runtime = None
            cost = -1

            if size is not None:
                # Run the executable
                output = run_executable(exe_path, t_sim)
                if output:
                    total_runtime, avg_runtime, cost = parse_results(output)
            else:
                print(f"Binary not found for {exe_name}, skipping...")

            res_dict[(solver, example)] = {
                "size": size,
                "total_runtime": total_runtime,
                "avg_runtime": avg_runtime,
                "cost": cost
            }

    print("solver,case,binary_size_kb,total_runtime_ms,avg_runtime_per_step_ms,J,normed_J")
    for example in examples:
        baseline_solver = ["qpOASES", "ALADIN"][0]
        baseline_cost = res_dict.get((baseline_solver, example), {}).get("cost", -1)
        for solver in solvers:
            result = res_dict.get((solver, example), {})
            size = result.get("size", "N/A")
            total_runtime = result.get("total_runtime", "N/A")
            avg_runtime = result.get("avg_runtime", "N/A")
            cost = result.get("cost", -1)
            normized_cost = cost / baseline_cost if baseline_cost > 0 and cost != -1 else "N/A"
            # Print CSV line
            print(f"{solver},{example},{int(size / 1000.0) if size is not None else 'N/A'},{total_runtime},{avg_runtime},{cost},{normized_cost}")
        print("")

if __name__ == "__main__":
    main()