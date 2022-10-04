import os
from distutils.core import Extension, setup

setup(
    name="corapy",
    version="0.1.0",
    description="Computational Raman Python tools",
    author="Mohammad Bagheri",
    author_email="Mohammad.Bagheri@oulu.fi",
    url="https://ramandb.oulu.fi",
    python_requires=">=3.7",
    install_requires=["numpy>=1.19.0","matplotlib>=2.2.2", "ase"],
    packages=[
        "corapy",
        "corapy.spectrum",
    ],
    scripts=[
        "scripts/corapy-plot",
    ],
    classifiers=[
          "Development Status :: 1 - Beta",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: BSD 3-Clause License",
          "Operating System :: OS Independent",
          "Programming Language :: Python :: 3.5",
          "Programming Language :: Python :: 3.6",
          "Programming Language :: Python :: 3.7",
          "Programming Language :: Python :: 3.8",
          "Programming Language :: Python :: 3.9",
          "Topic :: Scientific/Engineering :: Materials",
          "Topic :: Scientific/Engineering :: Physics",
          ],
)
