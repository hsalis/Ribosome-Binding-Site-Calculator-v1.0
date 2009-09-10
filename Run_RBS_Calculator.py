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

#Copyright 2008-2009 Howard Salis

from RBS_Calculator import RBS_Calculator
import sys, math

if __name__ == "__main__":
    #Read command line arguments
    input_args = []
    for arg in sys.argv:
        input_args.append(arg)

    start = 0

    if len(input_args) == 3:
        seq = input_args[1]
        start = input_args[2]
    elif len(input_args) == 2:
        seq = input_args[1]
    else:
        output = "Usage: " + input_args[0] + " [RNA/DNA Sequence] (start position)" + "\n"
        print output

    if start == 0:
        start_range = [0, len(seq)]
    else:
        start_range = [int(start), int(start)+1]

    name = "no name"

    #Create instance of RBS Calculator
    calcObj = RBS_Calculator(seq, start_range, name)
    calcObj.calc_dG()

    dG_total_list = calcObj.dG_total_list[:]
    start_pos_list = calcObj.start_pos_list[:]
    kinetic_score_list = calcObj.kinetic_score_list[:]

    expr_list = []
    for dG in dG_total_list:
        expr_list.append(calcObj.K * math.exp(-dG/calcObj.RT_eff))

    print len(expr_list)
    for (expr,start_pos,ks) in zip(expr_list,start_pos_list,kinetic_score_list):
        print start_pos, expr, ks