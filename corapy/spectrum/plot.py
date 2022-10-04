#! /usr/bin/env python

from argparse import ArgumentParser
from corapy.spectrum.reader import yaml_reader
from corapy.spectrum.intensity import calculate_intensity
from corapy.spectrum.plotter import plot_spectrum



def run():
    parser = ArgumentParser(description="Example: corapy-plot -i tensors.yaml -w 4.0")
    parser.set_defaults(
        width=4.0,
        function='lorentzian',
        save=None,
        spectype='intensity',
    )
    parser.add_argument("-i", "--input_file", help="tensors files downloaded from ramandb.oulu.fi")
    parser.add_argument("-st", "--spectype", help="type of spectra activity or intensity")
    parser.add_argument("-w", "--width", help="width of the Gaussian or Lorentzian function,default= 4 cm^-1 ", type=float)
    parser.add_argument("-f", "--function", help="name of broadening function(lorentzian or gaussian),default= lorentzian")
    parser.add_argument("-s", "--save", help="name file for saving plot with extension (e.g. plot.png)")
    args = parser.parse_args()
    
    freqs, tensors, num_activity = yaml_reader(file_name=args.input_file)
    activity, intensity = calculate_intensity(freqs,tensors, num_activity)
    if args.spectype == activity:
        plot_spectrum(freqs, activity, args.width, args.function, args.save)
    else:
        plot_spectrum(freqs, intensity, args.width, args.function, args.save)



if __name__ == "__main__":
    run()
