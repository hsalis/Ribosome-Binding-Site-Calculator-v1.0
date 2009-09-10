#!/usr/bin/env python

#This file is part of the Ribosome Binding Site Calculator.

#The Ribosome Binding Site Calculator is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#The Ribosome Binding Site Calculator is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with Ribosome Binding Site Calculator.  If not, see <http://www.gnu.org/licenses/>.

#Copyright 2008-2009 by Howard Salis.

from RBS_MC_Design import Monte_Carlo_Design
import sys, time

def all(x):
    """The all() function from Python 2.5, for backwards compatibility with previous versions of Python."""

    output = True
    for bool in x:
        if not bool:
            output = False
            break
    return output

def run():

    #Read command line arguments
    input_args = []
    for arg in sys.argv:
        input_args.append(arg)

    if len(input_args) == 4:
        pre_seq = input_args[1]
        post_seq = input_args[2]
        TIR = float(input_args[3])
    else:
        output = "Usage: " + input_args[0] + " [Pre-Sequence] [Post-Sequence] [Translation Initiation Rate]" + "\n"
        return (False,output)

    #Check Sequence and TIR Inputs
    valid_chars = ('A','C','T','G','U','a','c','t','g','u')
    valid_start_codons = ('ATG','AUG','GTG','GUG')

    if TIR <= 0: return (False,"invalid TIR")
    if len(pre_seq) <= 0 or not all([s in valid_chars for s in pre_seq[:]]): return (False,"invalid pre_seq")
    if len(post_seq) <= 0 or not all([s in valid_chars for s in post_seq[:]]) or post_seq[0:3].upper() not in valid_start_codons: return (False,"invalid post_seq")

    total_iterations=0
    MaxIter = 10000
    (TIR_out, RBS_out, estimator, iterations) = Monte_Carlo_Design(pre_seq,post_seq, None, TIR, MaxIter, verbose = False)
    total_iterations+=iterations

    #If RBS is not designed to specs within 10000 iterations, then start from a different initial condition
    while iterations == MaxIter:
        (TIR_out, RBS_out, estimator, iterations) = Monte_Carlo_Design(pre_seq,post_seq, None, TIR, MaxIter, verbose = False)
        total_iterations+=iterations

    output_string = "Program Executed" + "\n" + RBS_out + "\n" + str(TIR_out) + "\n" + str(total_iterations) + "\n"
    return (True,output_string)

def ReadOutput(output):
    import stringIO

    handle = stringIO.stringIO(output)

    for line in handle.readlines():
        if line == "Program Executed":
            RBS_out = handle.readlines()
            TIR_out = float(handle.readlines())
            total_iterations = int(handle.readlines())
            break
        elif line == "Program NOT Executed":
            RBS_out = ""
            TIR_out = -1.0
            total_iterations = 0
            break

    return (RBS_out, TIR_out, total_iterations)

if __name__ == "__main__":

    #Runs RBS Design optimization
    #Prints to standard output -- captured by globusrun interactive and transferred to the submitting computer.
    start = time.clock()

    (success,output) = run()

    end = time.clock()
    time_elapsed = end - start

    if success:
        print output
    else:
        print "Error"
        print output
