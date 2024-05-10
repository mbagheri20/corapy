import yaml
from yaml import Loader


def configs_reader(file_name):
    with open(file_name, 'r') as inputfile:
        data = yaml.load(inputfile, Loader=Loader)

    MPAPI_KEY = data['MPAPI_KEY']
    SLURM_CONFIGURATION = data['SLURM_CONFIGURATION']

    return [MPAPI_KEY, SLURM_CONFIGURATION[0]]
