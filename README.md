# CoRaPy (Computational Raman Python tools)

Python tools for computational raman studies (Developing in progress...)

Right now, only tools for plotting Raman spectra available which works with tensor files downloaded from the [Computational Raman Database](https://ramandb.oulu.fi/)

## Installation

### system requirements
* Python 3.X
* Numpy
* Scipy
* Phonopy
* ASE
* matplotlib


### Developer Installation

We highly recommend to use a python virtual enviroment for instaling the requirements and code package. to do this you can use this command:

```bash
python3 -m venv corapy-env
```
then activate virtual enviroment with following command:

```bash
source corapy-env/bin/activate
```
Get the source code from github:

```bash
git clone https://github.com/mbagheri20/corapy.git
cd corapy
```

You can use following command for installing system requirements packages:
```bash
pip install -r requirements.txt
```
Go to the directory with source files and run this command:

```bash
pip install -e .  

```

## How to use

Go to material page you would like to plot spectra on [Computational Raman Database](https://ramandb.oulu.fi/) and download the yaml file from Raman spectra tab "Download tensors for python"
Then you can use `corapy-plot -h` to see how to use.


## examples

There is an example in the example directory that can be used like this:

```bash
corapy-plot -i example/mp-1434.yaml
```
