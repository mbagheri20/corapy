import yaml
from yaml import Loader



def yaml_reader(file_name):
    """read yaml files contains raman tensors that downloaded from ramandb.oulu.fi

    Parameters:

    file_name: String
     name of the yaml file

    Returns:

    freqs: array
     list of frequncies
    tensors: array
     list of raman tensors
    num_activity: Integer
     number of activities
    """
    with open(file_name, 'r') as inputfile:
        data = yaml.load(inputfile, Loader=Loader)

    num_activity = len(data['raman_activities'])

    freqs = []
    tensors = []
    for i in range(num_activity):
        freqs.append(data['raman_activities'][i]['frequency(1/cm)'])
        tensors.append(data['raman_activities'][i]['raman_tensor'])

    return freqs, tensors, num_activity
