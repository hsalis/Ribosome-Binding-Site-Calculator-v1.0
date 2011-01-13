#Python wrapper for the Vienna RNA Package by Andreas R. Gruber, Ronny Lorenz, Stephan H. Bernhart, Richard Neuböck, and Ivo L. Hofacker (NAR, 2008).

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

#This Python wrapper is written by Howard Salis. Copyright 2008-2009 is owned by the University of California Regents. All rights reserved. :)
#Use at your own risk.

import os.path
import os, popen2, time

current_dir = os.path.dirname(os.path.abspath(__file__)) + "/tmp"
if not os.path.exists(current_dir): os.mkdir(current_dir)

debug=0

#Class that encapsulates all of the functions from NuPACK 2.0
class ViennaRNA(dict):

    debug_mode = 0
    RT = 0.61597 #gas constant times 310 Kelvin (in units of kcal/mol)

    def __init__(self,Sequence_List,material = "rna37"):

        self.ran = 0

        import re
        import random
        import string

        exp = re.compile('[ATGCU]',re.IGNORECASE)

        for seq in Sequence_List:
            if exp.match(seq) == None:
                error_string = "Invalid letters found in inputted sequences. Only ATGCU allowed. \n Sequence is \"" + str(seq) + "\"."
                raise ValueError(error_string)

        self["sequences"] = Sequence_List
        self["material"] = material

        random.seed(time.time())
        long_id = "".join([random.choice(string.letters + string.digits) for x in range(10)])
        self.prefix = current_dir + "/temp_" + long_id

    def mfe(self, strands,Temp = 37.0, dangles = "all",outputPS = False):

        self["mfe_composition"] = strands

        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")

        seq_string = "&".join(self["sequences"])
        input_string = seq_string + "\n"

        handle = open(self.prefix,"w")
        handle.write(input_string)
        handle.close()

        #Set arguments
        material = self["material"]
        if dangles is "none":
            dangles = " -d0 "
        elif dangles is "some":
            dangles = " -d1 "
        elif dangles is "all":
            dangles = " -d2 "

        if outputPS:
            outputPS_str = " "
        else:
            outputPS_str = " -noPS "

        #Call ViennaRNA C programs
        cmd = "RNAcofold"
        args = outputPS_str + dangles + " < " + self.prefix

        output = popen2.Popen3(cmd + args)
        #output.tochild.write(input_string)

        while output.poll() < 0:
            try:
                output.wait()
                time.sleep(0.001)
            except:
                break

        if debug == 1: print output.fromchild.read()

        #Skip the unnecessary output lines
        line = output.fromchild.readline()

        line = output.fromchild.readline()
        words = line.split(" ")
        bracket_string = words[0]
        (strands,bp_x, bp_y) = self.convert_bracket_to_numbered_pairs(bracket_string)

        energy = float(words[len(words)-1].replace(")","").replace("(","").replace("\n",""))

        self._cleanup()
        self["program"] = "mfe"
        self["mfe_basepairing_x"] = [bp_x]
        self["mfe_basepairing_y"] = [bp_y]
        self["mfe_energy"] = [energy]
        self["totalnt"]=strands


        #print "Minimum free energy secondary structure has been calculated."

    def subopt(self, strands,energy_gap,Temp = 37.0, dangles = "all", outputPS = False):

        self["subopt_composition"] = strands

        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")

        seq_string = "&".join(self["sequences"])
        input_string = seq_string + "\n"

        handle = open(self.prefix,"w")
        handle.write(input_string)
        handle.close()

        #Set arguments
        material = self["material"]
        if dangles is "none":
            dangles = " -d0 "
        elif dangles is "some":
            dangles = " -d1 "
        elif dangles is "all":
            dangles = " -d2 "

        if outputPS:
            outputPS_str = " "
        else:
            outputPS_str = " -noPS "

        #Call ViennaRNA C programs
        cmd = "RNAsubopt"
        args = " -e " + str(energy_gap) + outputPS_str + dangles + " < " + self.prefix

        output = popen2.Popen3(cmd + args)

        while output.poll() < 0:
            try:
                output.wait()
                time.sleep(0.001)
            except:
                break

        #print output.fromchild.read()
        if debug == 1: print output.fromchild.read()

        #Skip unnecessary line
        line = output.fromchild.readline()

        self["subopt_basepairing_x"] = []
        self["subopt_basepairing_y"] = []
        self["subopt_energy"] = []
        self["totalnt"]=[]
        counter=0

        while len(line)>0:
            line = output.fromchild.readline()
            if len(line) > 0:
                counter+=1
                words = line.split(" ")
                bracket_string = words[0]
                energy = float(words[len(words)-1].replace("\n",""))

                (strands,bp_x, bp_y) = self.convert_bracket_to_numbered_pairs(bracket_string)

                self["subopt_energy"].append(energy)
                self["subopt_basepairing_x"].append(bp_x)
                self["subopt_basepairing_y"].append(bp_y)

        self["subopt_NumStructs"] = counter

        self._cleanup()
        self["program"] = "subopt"

        #print "Minimum free energy and suboptimal secondary structures have been calculated."

    def energy(self, strands, base_pairing_x, base_pairing_y, Temp = 37.0, dangles = "all"):

        self["energy_composition"] = strands

        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")

        seq_string = "&".join(self["sequences"])
        strands = [len(seq) for seq in self["sequences"]]
        bracket_string = self.convert_numbered_pairs_to_bracket(strands,base_pairing_x,base_pairing_y)
        input_string = seq_string + "\n" + bracket_string + "\n"

        handle = open(self.prefix,"w")
        handle.write(input_string)
        handle.close()

        #Set arguments
        material = self["material"]
        if dangles is "none":
            dangles = " -d0 "
        elif dangles is "some":
            dangles = " -d1 "
        elif dangles is "all":
            dangles = " -d2 "

        #Call ViennaRNA C programs
        cmd = "RNAeval"
        args = dangles + " < " + self.prefix

        output = popen2.Popen3(cmd + args)

        while output.poll() < 0:
            try:
                output.wait()
                time.sleep(0.001)
            except:
                break

        #if debug == 1: print output.fromchild.read()

        self["energy_energy"] = []

        #Skip the unnecessary output lines
        line = output.fromchild.readline()

        line = output.fromchild.readline()
        words = line.split(" ")
        energy = float(words[len(words)-1].replace("(","").replace(")","").replace("\n",""))

        self["program"] = "energy"
        self["energy_energy"].append(energy)
        self["energy_basepairing_x"] = [base_pairing_x]
        self["energy_basepairing_y"] = [base_pairing_y]
        self._cleanup()

        return energy

    def convert_bracket_to_numbered_pairs(self,bracket_string):

        all_seq_len = len(bracket_string)
        bp_x = []
        bp_y = []
        strands = []

        for y in range(bracket_string.count(")")):
            bp_y.append([])

        last_nt_x_list = []
        counter=0
        num_strands=0
        for (pos,letter) in enumerate(bracket_string[:]):
            if letter is ".":
                counter += 1

            elif letter is "(":
                bp_x.append(pos-num_strands)
                last_nt_x_list.append(pos-num_strands)
                counter += 1

            elif letter is ")":
                nt_x = last_nt_x_list.pop()
                nt_x_pos = bp_x.index(nt_x)
                bp_y[nt_x_pos] = pos-num_strands
                counter += 1

            elif letter is "&":
                strands.append(counter)
                counter=0
                num_strands+=1

            else:
                print "Error! Invalid character in bracket notation."

        if len(last_nt_x_list) > 0:
            print "Error! Leftover unpaired nucleotides when converting from bracket notation to numbered base pairs."

        strands.append(counter)
        bp_x = [pos+1 for pos in bp_x[:]] #Shift so that 1st position is 1
        bp_y = [pos+1 for pos in bp_y[:]] #Shift so that 1st position is 1

        return (strands,bp_x, bp_y)

    def convert_numbered_pairs_to_bracket(self,strands,bp_x,bp_y):

        bp_x = [pos-1 for pos in bp_x[:]] #Shift so that 1st position is 0
        bp_y = [pos-1 for pos in bp_y[:]] #Shift so that 1st position is 0

        bracket_notation = []
        counter=0
        for (strand_number,seq_len) in enumerate(strands):
            if strand_number > 0: bracket_notation.append("&")
            for pos in range(counter,seq_len+counter):
                if pos in bp_x:
                    bracket_notation.append("(")
                elif pos in bp_y:
                    bracket_notation.append(")")
                else:
                    bracket_notation.append(".")
            counter+=seq_len

        return "".join(bracket_notation)

    def _cleanup(self):

        if os.path.exists(self.prefix): os.remove(self.prefix)
        return

if __name__ == "__main__":

    sequences = ["AGGGGGGATCTCCCCCCAAAAAATAAGAGGTACACATGACTAAAACTTTCAAAGGCTCAGTATTCCCACT"] #,"acctcctta"]
    test = ViennaRNA(sequences,material = "rna37")

    test.mfe([1],Temp = 37.0, dangles = "all")

    bp_x = test["mfe_basepairing_x"][0]
    bp_y = test["mfe_basepairing_y"][0]
    strands = test["totalnt"]
    bracket_string = test.convert_numbered_pairs_to_bracket(strands,bp_x,bp_y)
    print bracket_string

    (strands,bp_x, bp_y) = test.convert_bracket_to_numbered_pairs(bracket_string)

    print "Strands = ", strands
    print "bp_x = ", bp_x
    print "bp_y = ", bp_y

    print test.energy(strands, bp_x, bp_y, dangles = "all")
    test.subopt(strands,3.5,dangles = "all")
    print test

#    print bracket_string
#    print test.convert_numbered_pairs_to_bracket(strands,bp_x,bp_y)
