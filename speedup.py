import subprocess
import json
from collections import defaultdict

N_RUNS = 3 

serial_program = "serial/foxes-rabbits"
parallel_program = "omp/foxes-rabbits-omp"

params = ["300 6000 900 8000 200000 7 100000 12 20 9999",
        "4000 900 2000 100000 1000000 10 400000 30 30 12345",
        "20000 1000 800 100000 80000 10 1000 30 8 500",
        "100000 200 500 500 1000 3 600 6 10 1234"]   

results = defaultdict(lambda: defaultdict(lambda : defaultdict(list)))

print(results)
for par in params:
    for version,name in ((serial_program,"serial"), (parallel_program, "omp")):
        for _ in range(N_RUNS):
            p = subprocess.run([version] + par.split(),capture_output=True)
            time = float(p.stderr.decode().split("s")[0]) 
            out = p.stdout.decode()
            results[par][name]["time"].append(time)
            results[par][name]["time"].append(time)
            results[par][name]["stdout"].append(out.strip())
            with open("speedup.json", "w") as f:
                json.dump(results, f, indent=2)
