import math
import os
import pickle

try:
    import libconf
except ImportError:
    libconf = None


def get_varied_keys(cases):
    first_case = cases[0]
    varied_keys = set()
    for case in cases[1:]:
        for key, val in case.items():
            if first_case.get(key, None) != val and '/' not in str(val) and ':' not in str(val) and key != 'tasks_per_node':
              varied_keys.add(key)
    return list(sorted(varied_keys))


def gen_folders(campaign_name, cases, generators):
    varied_keys = get_varied_keys(cases)

    def get_case_dir(cfg):
        out = []
        for key in varied_keys:
            out.append(key + str(cfg[key]))

        if len(out) == 0:
            return 'case'
        else:
            return os.path.join(campaign_name, "_".join(out))

    # generate folders
    for case in cases:
        job_dir = get_case_dir(case)
        if os.path.exists(job_dir):
            print("Skipping gen_folders for case '{}' (folder already exists).".format(job_dir))
            continue
        else:
            print("Creating folder for case '{}'".format(job_dir))

        os.makedirs(job_dir)

        # compute nnodes automatically from tasks per node + number of tasks (used in job script)
        if 'ntasks' in case:
            case['nnodes'] = int(math.ceil(float(case['ntasks']) / float(case['tasks_per_node'])))

        # save case data as a Python pickle file
        with open(os.path.join(job_dir, 'case.pkl'), 'wb') as f:
            pickle.dump(case, f)

        for generator in generators:
            generator(case, job_dir)
