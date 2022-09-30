# setup.py
import setuptools


setuptools.setup(
name="corapy",
version="0.1",
author="Mohammad Bagheri",
author_email="mohammad.bagheri@oulu.fi",
description="Computational Raman Python tools",
packages=["corapy","corapy/spectrum"],
install_requires=["numpy","phonopy"],
package_data={"": ["README.md","LICENSE","doc","example"]},
url="http://ramandb.oulu.fi/",
classifiers=[
"Development Status :: ",
"Intended Audience :: Science/Research",
"License :: OSI Approved :: Creative Commons Attribution-ShareAlike 4.0 International License",
"Operating System :: OS Independent",
"Programming Language :: Python :: 3.5",
"Programming Language :: Python :: 3.6",
"Programming Language :: Python :: 3.7",
"Programming Language :: Python :: 3.8",
"Programming Language :: Python :: 3.9",
"Topic :: Scientific/Engineering :: Materials",
"Topic :: Scientific/Engineering :: Physics"]
)
