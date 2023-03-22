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
```

You can use following command for installing system requirements packages:
```bash
cd corapy
pip install -r requirements.txt
```
To install the code go to the directory with source files and run this command:

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

## How to cite

If you have used Corapy, please cite the following article:

- "High-throughput computation of Raman spectra from first principles",

  Mohammad Bagheri & Hannu-Pekka Komsa, Sci. Data, **10**, 80 (2023)

  https://doi.org/10.1038/s41597-023-01988-5 (Open access)

  ```
  @article{corapy,
     doi = {10.1038/s41597-023-01988-5},
     url = {https://doi.org/10.1038/s41597-023-01988-5},
     journal = {Scientific Data},
     year = {2023},
     title = {High-throughput computation of Raman spectra from first principles},
     author = {Mohammad Bagheri and Hannu-Pekka Komsa},
     pages = {80},
     issue = {14},
     volume = {10}
  }
  ```
