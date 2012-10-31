#Python wrapper for NUPACK 2.0 by Dirks, Bois, Schaeffer, Winfree, and Pierce (SIAM Review)

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
import os, popen2, time, random, string

tempdir = "/tmp" + "".join([random.choice(string.digits) for x in range(6)])

current_dir = os.path.dirname(os.path.abspath(__file__)) + tempdir
if not os.path.exists(current_dir): os.mkdir(current_dir)

debug=0

#Class that encapsulates all of the functions from NuPACK 2.0
class NuPACK(dict):

    debug_mode = 0
    RT = 0.61597 #gas constant times 310 Kelvin (in units of kcal/mol)

    def __init__(self,Sequence_List,material):

        self.ran = 0

        import re
        import string

        exp = re.compile('[ATGCU]',re.IGNORECASE)

        for seq in Sequence_List:
            if exp.match(seq) == None:
                error_string = "Invalid letters found in inputted sequences. Only ATGCU allowed. \n Sequence is \"" + str(seq) + "\"."
                raise ValueError(error_string)

        if not material == 'rna' and not material == 'dna' and not material == "rna1999": raise ValueError("The energy model must be specified as either ""dna"", ""rna"", or ""rna1999"" .")

        self["sequences"] = Sequence_List
        self["material"] = material

        random.seed(time.time())
        long_id = "".join([random.choice(string.letters + string.digits) for x in range(10)])
        self.prefix = current_dir + "/nu_temp_" + long_id

    def complexes(self,MaxStrands, Temp = 37.0, ordered = "", pairs = "", mfe = "", degenerate = "", dangles = "some", timeonly = "", quiet="", AdditionalComplexes = []):
        """A wrapper for the complexes command, which calculates the equilibrium probability of the formation of a multi-strand
        RNA or DNA complex with a user-defined maximum number of strands.  Additional complexes may also be included by the user."""

        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")
        if int(MaxStrands) <= 0: raise ValueError("The maximum number of strands must be greater than zero.")

        #Write input files
        self._write_input_complexes(MaxStrands, AdditionalComplexes)

        #Set arguments
        material = self["material"]
        if ordered: ordered = " -ordered "
        if pairs: pairs = " -pairs "
        if mfe: mfe = " -mfe "
        if degenerate: degenerate = " -degenerate "
        if timeonly: timeonly = " -timeonly "
        if quiet: quiet = " -quiet "
        dangles = "-dangles " + dangles + " "


        #Call NuPACK C programs
        cmd = "complexes"
        args = " -T " + str(Temp) + " -material " + material + " " + ordered + pairs + mfe + degenerate \
        + dangles + timeonly + quiet + " "

        output = popen2.Popen3(cmd + args + self.prefix)
        while output.poll() < 0:
            try:
                output.wait()
                time.sleep(0.001)
            except:
                break

        if debug == 1: print output.fromchild.read()

        #Read output files
        self._read_output_cx()
        self._cleanup("cx")


        if ordered:
            self._read_output_ocx()
            self._read_output_ocx_mfe()

            self._cleanup("ocx")
            self._cleanup("ocx-mfe")
            self._cleanup("ocx-key")
        self._cleanup("in")

        #print "Complex energies and secondary structures calculated."
        self.ran = 1
        self["program"] = "complexes"

    def mfe(self, strands,Temp = 37.0, multi = " -multi", pseudo = "", degenerate = "", dangles = "some"):

        self["mfe_composition"] = strands

        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")

        if (multi == 1 and pseudo == 1): raise ValueError("The pseudoknot algorithm does not work with the -multi option.")

        #Write input files
        self._write_input_mfe(strands)

        #Set arguments
        material = self["material"]
        if multi == "": multi = ""
        if pseudo: pseudo = " -pseudo"
        if degenerate: degenerate = " -degenerate "
        dangles = " -dangles " + dangles + " "

        #Call NuPACK C programs
        cmd = "mfe"
        args = " -T " + str(Temp) + multi + pseudo + " -material " + material + degenerate + dangles + " "

        output = popen2.Popen3(cmd + args + self.prefix)
        while output.poll() < 0:
            try:
                output.wait()
                time.sleep(0.001)
            except:
                break

        if debug == 1: print output.fromchild.read()

        self._read_output_mfe()
        self._cleanup("mfe")
        self._cleanup("in")
        self["program"] = "mfe"
        #print "Minimum free energy secondary structure has been calculated."


    def subopt(self, strands,energy_gap,Temp = 37.0, multi = " -multi", pseudo = "", degenerate = "", dangles = "some"):

        self["subopt_composition"] = strands

        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")

        if (multi == 1 and pseudo == 1): raise ValueError("The pseudoknot algorithm does not work with the -multi option.")

        #Write input files
        self._write_input_subopt(strands,energy_gap)

        #Set arguments
        material = self["material"]
        if multi == "": multi = ""
        if pseudo: pseudo = " -pseudo"
        if degenerate: degenerate = " -degenerate "
        dangles = " -dangles " + dangles + " "

        #Call NuPACK C programs
        cmd = "subopt"
        args = " -T " + str(Temp) + multi + pseudo + " -material " + material + degenerate + dangles + " "

        output = popen2.Popen3(cmd + args + self.prefix)
        while output.poll() < 0:
            try:
                output.wait()
                time.sleep(0.001)
            except:
                break

        if debug == 1: print output.fromchild.read()

        self._read_output_subopt()
        self._cleanup("subopt")
        self._cleanup("in")
        self["program"] = "subopt"

        #print "Minimum free energy and suboptimal secondary structures have been calculated."

    def energy(self, strands, base_pairing_x, base_pairing_y, Temp = 37.0, multi = " -multi", pseudo = "", degenerate = "", dangles = "some"):

        self["energy_composition"] = strands

        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")

        if (multi == 1 and pseudo == 1): raise ValueError("The pseudoknot algorithm does not work with the -multi option.")

        #Write input files
        self._write_input_energy(strands,base_pairing_x,base_pairing_y)

        #Set arguments
        material = self["material"]
        if multi == "": multi = ""
        if pseudo: pseudo = " -pseudo"
        if degenerate: degenerate = " -degenerate "
        dangles = " -dangles " + dangles + " "

        #Call NuPACK C programs
        cmd = "energy"
        args = " -T " + str(Temp) + multi + pseudo + " -material " + material + degenerate + dangles + " "

        output = popen2.Popen3(cmd + args + self.prefix)
        while output.poll() < 0:
            try:
                output.wait()
                time.sleep(0.001)
            except:
                break

        #if debug == 1: print output.fromchild.read()

        self["energy_energy"] = []

        #Skip the comments of the text file
        line = output.fromchild.readline()
        while line[0]=="%":
            line = output.fromchild.readline()


        energy = float(line)
        self["program"] = "energy"
        self["energy_energy"].append(energy)
        self["energy_basepairing_x"] = [base_pairing_x]
        self["energy_basepairing_y"] = [base_pairing_y]
        self._cleanup("in")

        return energy

    def pfunc(self, strands, Temp = 37.0, multi = " -multi", pseudo = "", degenerate = "", dangles = "some"):

        self["pfunc_composition"] = strands

        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")

        if (multi == 1 and pseudo == 1): raise ValueError("The pseudoknot algorithm does not work with the -multi option.")

        #Write input files
        #Input for pfunc is the same as mfe
        self._write_input_mfe(strands)

        #Set arguments
        material = self["material"]
        if multi == "": multi = ""
        if pseudo: pseudo = " -pseudo"
        if degenerate: degenerate = " -degenerate "
        dangles = " -dangles " + dangles + " "

        #Call NuPACK C programs
        cmd = "pfunc"
        args = " -T " + str(Temp) + multi + pseudo + " -material " + material + degenerate + dangles + " "

        output = popen2.Popen3(cmd + args + self.prefix)
        while output.poll() < 0:
            try:
                output.wait()
                time.sleep(0.001)
            except:
                break

        #if debug == 1: print output.fromchild.read()

        #Skip the comments of the text file
        line = output.fromchild.readline()
        words = line.split(" ")
        while line[0]=="%" or words[0] == "Attempting":
            line = output.fromchild.readline()
            words = line.split(" ")

        energy = float(line)

        line = output.fromchild.readline()
        partition_function = float(line)

        self["program"] = "pfunc"
        self["pfunc_energy"] = energy
        self["pfunc_partition_function"] = partition_function
        self._cleanup("in")

        return partition_function

    def count(self, strands, Temp = 37.0, multi = " -multi", pseudo = "", degenerate = "", dangles = "some"):

        self["count_composition"] = strands

        if (multi == 1 and pseudo == 1): raise ValueError("The pseudoknot algorithm does not work with the -multi option.")

        #Write input files
        #Input for count is the same as mfe
        self._write_input_mfe(strands)

        #Set arguments
        material = self["material"]
        if multi == "": multi = ""
        if pseudo: pseudo = " -pseudo"
        if degenerate: degenerate = " -degenerate "
        dangles = " -dangles " + dangles + " "

        #Call NuPACK C programs
        cmd = "count"
        args = " -T " + str(Temp) + multi + pseudo + " -material " + material + degenerate + dangles + " "

        output = popen2.Popen3(cmd + args + self.prefix)
        while output.poll() < 0:
            try:
                output.wait()
                time.sleep(0.001)
            except:
                break

        #if debug == 1: print output.fromchild.read()

        #Skip the comments of the text file
        line = output.fromchild.readline()
        words = line.split(" ")
        while line[0]=="%" or words[0] == "Attempting":
            line = output.fromchild.readline()
            words = line.split(" ")

        number = float(line)

        self["program"] = "count"
        self["count_number"] = number
        self._cleanup("in")

        return number

    def _write_input_energy(self,strands,base_pairing_x,base_pairing_y):
        #Creates the input file for energy NUPACK functions
        #strands is a list containing the number of each strand in the complex (assumes -multi flag is used)
        #base_pairing_x and base_pairing_y is a list of base pairings of the strands s.t. #x < #y are base paired

        NumStrands = len(self["sequences"])
        input_str = str(NumStrands) + "\n"
        for seq in self["sequences"]:
            input_str = input_str + seq + "\n"

        NumEachStrands = ""
        for num in strands:
            NumEachStrands = NumEachStrands + str(num) + " "

        input_str = input_str + NumEachStrands + "\n"
        for pos in range(len(base_pairing_x)):
            input_str = input_str + str(base_pairing_x[pos]) + "\t" + str(base_pairing_y[pos]) + "\n"

        handle = open(self.prefix + ".in", "w")
        handle.writelines(input_str)
        handle.close()

    def _write_input_subopt(self,strands,energy_gap):
        #Creates the input file for mfe and subopt NUPACK functions
        #strands is a list containing the number of each strand in the complex (assumes -multi flag is used)

        NumStrands = len(self["sequences"])
        input_str = str(NumStrands) + "\n"
        for seq in self["sequences"]:
            input_str = input_str + seq + "\n"

        NumEachStrands = ""
        for num in strands:
            NumEachStrands = NumEachStrands + str(num) + " "

        input_str = input_str + NumEachStrands + "\n"
        input_str = input_str + str(energy_gap) + "\n"

        handle = open(self.prefix + ".in", "w")
        handle.writelines(input_str)
        handle.close()

    def _write_input_mfe(self,strands):
        #Creates the input file for mfe and subopt NUPACK functions
        #strands is a list containing the number of each strand in the complex (assumes -multi flag is used)

        NumStrands = len(self["sequences"])
        input_str = str(NumStrands) + "\n"
        for seq in self["sequences"]:
            input_str = input_str + seq + "\n"

        NumEachStrands = ""
        for num in strands:
            NumEachStrands = NumEachStrands + str(num) + " "

        input_str = input_str + NumEachStrands + "\n"

        handle = open(self.prefix + ".in", "w")
        handle.writelines(input_str)
        handle.close()

    def _write_input_complexes(self, MaxStrands, AdditionalComplexes = [] ):

        #First, create the input string for file.in to send into NUPACK
        NumStrands = len(self["sequences"])
        input_str = str(NumStrands) + "\n"
        for seq in self["sequences"]:
            input_str = input_str + seq + "\n"
        input_str = input_str + str(MaxStrands) + "\n"

        handle = open(self.prefix + ".in", "w")
        handle.writelines(input_str)
        handle.close()

        if len(AdditionalComplexes) > 0:
            #The user may also specify additional complexes composed of more than MaxStrands strands. Create the input string detailing this.
            counter=0
            counts = [[]]
            added=[]
            for (complexes,i) in zip(AdditionalComplexes,range(len(AdditionalComplexes))):

                if len(complexes) <= MaxStrands: #Remove complexes if they have less than MaxStrands strands.
                    AdditionalComplexes.pop(i)
                else:
                    counts.append([])
                    added.append(0)
                    for j in range(NumStrands): #Count the number of each unique strand in each complex and save it to counts
                        counts[counter].append(complexes.count(j+1))
                    counter += 1

            list_str = ""
            for i in range(len(counts)-1):
                if added[i] == 0:
                    list_str = list_str + "C " + " ".join([str(count) for count in counts[i]]) + "\n"
                    list_str = list_str + " ".join([str(strand) for strand in AdditionalComplexes[i]]) + "\n"
                    added[i] = 1
                    for j in range(i+1,len(counts)-1):
                        if counts[i] == counts[j] and added[j] == 0:
                            list_str = list_str + " ".join([str(strand) for strand in AdditionalComplexes[j]]) + "\n"
                            added[j]=1

            handle = open(self.prefix + ".list", "w")
            handle.writelines(list_str)
            handle.close()

    def _read_output_cx(self):

        #Read the prefix.cx output text file generated by NuPACK and write its data to instanced attributes
    #Output: energies of unordered complexes in key "unordered_energies"
    #Output: strand composition of unordered complexes in key "unordered_complexes"

        handle = open(self.prefix+".cx", "rU")

        line = handle.readline()

        #Read some useful data from the comments of the text file
        while line[0]=="%":

            words=line.split()

            if len(words) > 7 and words[1] == "Number" and words[2] == "of" and words[3] == "complexes" and words[4] == "from" and words[5] == "enumeration:":
                self["numcomplexes"] =int(words[6])

            elif len(words) > 8 and words[1] == "Total" and words[2] == "number" and words[3] == "of" and words[4] =="permutations" and words[5] == "to" and words[6] == "calculate:":
                self["num_permutations"] = int(words[7])

            line = handle.readline()

        self["unordered_energies"] = []
        self["unordered_complexes"] = []
        self["unordered_composition"] = []

        while line:
            words=line.split()

            if not words[0] == "%":

                complex = words[0]
                strand_compos = [int(f) for f in words[1:len(words)-1]]
                energy = float(words[len(words)-1])

                self["unordered_complexes"].append(complex)
                self["unordered_energies"].append(energy)
                self["unordered_composition"].append(strand_compos)

            line = handle.readline()
        handle.close()

    def _read_output_ocx(self):

    #Read the prefix.ocx output text file generated by NuPACK and write its data to instanced attributes
    #Output: energies of ordered complexes in key "ordered_energies"
    #Output: number of permutations and strand composition of ordered complexes in key "ordered_complexes"

        handle = open(self.prefix+".ocx", "rU")

        line = handle.readline()

        #Read some useful data from the comments of the text file
        while line[0]=="%":

            words=line.split()

            if len(words) > 7 and words[1] == "Number" and words[2] == "of" and words[3] == "complexes" and words[4] == "from" and words[5] == "enumeration:":
                self["numcomplexes"] =int(words[6])

            elif len(words) > 8 and words[1] == "Total" and words[2] == "number" and words[3] == "of" and words[4] =="permutations" and words[5] == "to" and words[6] == "calculate:":
                self["num_permutations"] = int(words[7])

            line = handle.readline()

        self["ordered_complexes"] = []
        self["ordered_energies"] = []
        self["ordered_permutations"] = []
        self["ordered_composition"] = []


        while line:
            words=line.split()

            if not words[0] == "%":

                complex = words[0]
                permutations = words[1]
                strand_compos = [int(f) for f in words[2:len(words)-1]]
                energy = float(words[len(words)-1])

                self["ordered_complexes"].append(complex)
                self["ordered_permutations"].append(permutations)
                self["ordered_energies"].append(energy)
                self["ordered_composition"].append(strand_compos)

            line = handle.readline()
        handle.close()

    def _read_output_ocx_mfe(self):
    #Read the prefix.ocx output text file generated by NuPACK and write its data to instanced attributes
    #Output: energy of mfe of each complex in key "ordered_energy"


        #Make sure that the ocx file has already been read.
        if not (self.has_key("ordered_complexes") and self.has_key("ordered_permutations") and self.has_key("ordered_energies") and self.has_key("ordered_composition")):
            self._read_output_ocx(self,prefix)

        handle = open(self.prefix+".ocx-mfe", "rU")

        #Skip the comments of the text file
        line = handle.readline()
        while line[0]=="%":
            line = handle.readline()

        self["ordered_basepairing_x"] = []
        self["ordered_basepairing_y"] = []
        self["ordered_energy"] = []
        self["ordered_totalnt"]=[]

        while line:
            words=line.split()

            if not line == "\n" and not words[0] == "%" and not words[0] == "":

                #Read the line containing the number of total nucleotides in the complex
                totalnt = words[0]

                self["ordered_totalnt"].append(totalnt)

                #Read the line containing the mfe
                words = handle.readline().split()
                mfe = float(words[0])

                self["ordered_energy"].append(mfe)

                #Skip the line containing the dot/parens description of the secondary structure
                line = handle.readline()

                #Read in the lines containing the base pairing description of the secondary structure
                #Continue reading until a % comment
                bp_x = []
                bp_y = []

                line = handle.readline()
                words = line.split()
                while not line == "\n" and not words[0] == "%":
                    bp_x.append(int(words[0]))
                    bp_y.append(int(words[1]))
                    words = handle.readline().split()

                self["ordered_basepairing_x"].append(bp_x)
                self["ordered_basepairing_y"].append(bp_y)

            line = handle.readline()
        handle.close()

    def _read_output_mfe(self):
    #Read the prefix.mfe output text file generated by NuPACK and write its data to instanced attributes
    #Output: total sequence length and minimum free energy
    #Output: list of base pairings describing the secondary structure

        handle = open(self.prefix+".mfe", "rU")

        #Skip the comments of the text file
        line = handle.readline()
        while line[0]=="%":
            line = handle.readline()

        self["mfe_basepairing_x"] = []
        self["mfe_basepairing_y"] = []
        self["mfe_energy"] = []
        self["totalnt"]=[]

        counter = 0
        while line:
            words=line.split()

            if not line == "\n" and not words[0] == "%" and not words[0] == "":

                #Read the line containing the number of total nucleotides in the complex
                totalnt = words[0]

                self["totalnt"].append(totalnt)
                counter += 1

                #Read the line containing the mfe
                words = handle.readline().split()
                mfe = float(words[0])

                self["mfe_energy"].append(mfe)

                #Skip the line containing the dot/parens description of the secondary structure
                line = handle.readline()

                #Read in the lines containing the base pairing description of the secondary structure
                #Continue reading until a % comment
                bp_x = []
                bp_y = []

                line = handle.readline()
                words = line.split()
                while not line == "\n" and not words[0] == "%":
                    bp_x.append(int(words[0]))
                    bp_y.append(int(words[1]))
                    words = handle.readline().split()

                self["mfe_basepairing_x"].append(bp_x)
                self["mfe_basepairing_y"].append(bp_y)

            line = handle.readline()
        handle.close()

        self["mfe_NumStructs"] = counter

    def _read_output_subopt(self):
    #Read the prefix.subopt output text file generated by NuPACK and write its data to instanced attributes
    #Output: total sequence length and minimum free energy
    #Output: list of base pairings describing the secondary structure

        handle = open(self.prefix+".subopt", "rU")

        #Skip the comments of the text file
        line = handle.readline()
        while line[0]=="%":
            line = handle.readline()

        self["subopt_basepairing_x"] = []
        self["subopt_basepairing_y"] = []
        self["subopt_energy"] = []
        self["totalnt"]=[]

        counter=0

        while line:
            words=line.split()

            if not line == "\n" and not words[0] == "%" and not words[0] == "":

                #Read the line containing the number of total nucleotides in the complex
                totalnt = words[0]

                self["totalnt"].append(totalnt)
                counter += 1

                #Read the line containing the mfe
                words = handle.readline().split()
                mfe = float(words[0])

                self["subopt_energy"].append(mfe)

                #Skip the line containing the dot/parens description of the secondary structure
                line = handle.readline()

                #Read in the lines containing the base pairing description of the secondary structure
                #Continue reading until a % comment
                bp_x = []
                bp_y = []

                line = handle.readline()
                words = line.split()
                while not line == "\n" and not words[0] == "%":
                    bp_x.append(int(words[0]))
                    bp_y.append(int(words[1]))
                    words = handle.readline().split()

                self["subopt_basepairing_x"].append(bp_x)
                self["subopt_basepairing_y"].append(bp_y)

            line = handle.readline()
        handle.close()

        self["subopt_NumStructs"] = counter

    def _cleanup(self,suffix):

        if os.path.exists(self.prefix+"."+suffix): os.remove(self.prefix+"."+suffix)

        return

    def export_PDF(self, complex_ID, name = "", filename = "temp.pdf", program = None):
        """Uses Zuker's sir_graph_ng and ps2pdf.exe to convert a secondary structure described in .ct format
        to a PDF of the RNA"""

        if program is None:
            program = self["program"]

        inputfile = "temp.ct"
        self.Convert_to_ct(complex_ID,name,inputfile,program)


        cmd = "sir_graph_ng" #Assumes it's on the path
        args = "-p" #to PostScript file
        output = popen2.Popen3(cmd + " " + args + " " + inputfile,"r")
        output.wait()
        if debug == 1: print output.fromchild.read()

        inputfile = inputfile[0:len(inputfile)-2] + "ps"

        cmd = "ps2pdf" #Assumes it's on the path
        output = popen2.Popen3(cmd + " " + inputfile,"r")
        output.wait()
        if debug == 1: print output.fromchild.read()

        outputfile = inputfile[0:len(inputfile)-2] + "pdf"

        #Remove the temporary file "temp.ct" if it exists
        if os.path.exists("temp.ct"): os.remove("temp.ct")

        #Remove the temporary Postscript file if it exists
        if os.path.exists(inputfile): os.remove(inputfile)

        #Rename the output file to the desired filename.
        if os.path.exists(outputfile): os.rename(outputfile,filename)
        #Done!

    def Convert_to_ct(self,complex_ID,name,filename = "temp.ct",program = "ordered"):
        """Converts the secondary structure of a single complex into the .ct file format, which is used
        with sir_graph_ng (or other programs) to create an image of the secondary structure."""

        #hacksy way of reading from data produced by 'complex', by 'mfe', or by 'subopt'
        data_x = program + "_basepairing_x"
        data_y = program + "_basepairing_y"
        mfe_name = program + "_energy"
        composition_name = program + "_composition"

        #Format of .ct file

        #Header: <Total # nt> \t dG = <# mfe> kcal/mol \t <name of sequence>
        #The Rest:
        #<nt num> \t <bp letter> \t <3' neighbor> \t <5' neighbor> \t <# of bp'ing, 0 if none> \t ...
        #<strand-specific nt num> \t <3' neighbor if connected by helix> \t <5' neighbor if connected by helix>

        #Extract the data for the desired complex using complex_ID
        bp_x = self[data_x][complex_ID]
        bp_y = self[data_y][complex_ID]
        mfe = self[mfe_name][complex_ID]

        if program == "mfe" or program == "subopt" or program == "energy":
            composition = self[composition_name]
        elif program == "ordered" or program == "unordered":
            composition = self[composition_name][complex_ID]


        #Determine concatenated sequence of all strands, their beginnings, and ends
        allseq = ""
        strand_begins = []
        strand_ends = []

        #Seemingly, the format of the composition is different for the program complex vs. mfe/subopt
        #for mfe/subopt, the composition is the list of strand ids
        #for complex, it is the number of each strand (in strand id order) in the complex
        #for mfe/subopt, '1 2 2 3' refers to 1 strand of 1, 2 strands of 2, and 1 strand of 3.
        #for complex, '1 2 2 3' refers to 1 strand of 1, 2 strands of 2, 2 strands of 3, and 3 strands of 4'.
        #what a mess.

        if program == "mfe" or program == "subopt" or program == "energy":
            for strand_id in composition:
                strand_begins.append(len(allseq) + 1)
                allseq = allseq + self["sequences"][strand_id-1]
                strand_ends.append(len(allseq))

        else:
            for (num_strands,strand_id) in zip(composition,range(len(composition))):
                for j in range(num_strands):
                    strand_begins.append(len(allseq) + 1)
                    allseq = allseq + self["sequences"][strand_id]
                    strand_ends.append(len(allseq))

        seq_len = len(allseq)

        #print "Seq Len = ", seq_len, "  Composition = ", composition
        #print "Sequence = ", allseq
        #print "Base pairing (x) = ", bp_x
        #print "Base pairing (y) = ", bp_y


        #Create the header
        header = str(seq_len) + "\t" + "dG = " + str(mfe) + " kcal/mol" + "\t" + name + "\n"

        #Open the file
        handle = open(filename,"w")

        #Write the header
        handle.write(header)

        #Write a line for each nt in the secondary structure
        for i in range(1,seq_len+1):


            for (nt,pos) in zip(strand_begins,range(len(strand_begins))):
                if i >= nt:
                    strand_id = pos


            #Determine 3' and 5' neighbor
            #If this is the beginning of a strand, then the 3' neighbor is 0
            #If this is the end of a strand, then the 5' neighbor is 0

            if i in strand_begins:
                nb_5p = 0
            else:
                nb_5p = i - 1

            if i in strand_ends:
                nb_3p = 0
            else:
                nb_3p = i + 1


            if i in bp_x or i in bp_y:
                if i in bp_x: nt_bp = bp_y[bp_x.index(i)]
                if i in bp_y: nt_bp = bp_x[bp_y.index(i)]
            else:
                nt_bp = 0

            #Determine strand-specific counter
            strand_counter = i - strand_begins[strand_id] + 1

            #Determine the 3' and 5' neighbor helical connectivity
            #If the ith nt is connected to its 3', 5' neighbor by a helix, then include it
            #Otherwise, 0
            #Helix connectivity conditions:
            #The 5' or 3' neighbor is connected via a helix iff:
            #a) helix start: i not bp'd, i+1 bp'd, bp_id(i+1) - 1 is bp'd, bp_id(i+1) + 1 is not bp'd
            #b) helix end: i not bp'd, i-1 bp'd, bp_id(i-1) - 1 is not bp'd, bp_id(i-1) + 1 is bp'd
            #c) helix continued: i and bp_id(i)+1 is bp'd, 5' helix connection is bp_id(bp_id(i)+1)
            #d) helix continued: i and bp_id(i)-1 is bp'd, 3' helix connection is bp_id(bp_id(i)-1)
            #Otherwise, zero.

            #Init
            hc_5p = 0
            hc_3p = 0

            if i in bp_x or i in bp_y:  #helix continued condition (c,d)
                if i in bp_x: bp_i = bp_y[bp_x.index(i)]
                if i in bp_y: bp_i = bp_x[bp_y.index(i)]

                if bp_i+1 in bp_x or bp_i+1 in bp_y: #helix condition c
                    if bp_i+1 in bp_x: hc_3p = bp_y[bp_x.index(bp_i+1)]
                    if bp_i+1 in bp_y: hc_3p = bp_x[bp_y.index(bp_i+1)]

                if bp_i-1 in bp_x or bp_i-1 in bp_y: #helix condition d
                    if bp_i-1 in bp_x: hc_5p = bp_y[bp_x.index(bp_i-1)]
                    if bp_i-1 in bp_y: hc_5p = bp_x[bp_y.index(bp_i-1)]

            else: #helix start or end (a,b)

                if i+1 in bp_x or i+1 in bp_y: #Start, condition a
                    if i+1 in bp_x: bp_3p = bp_y[bp_x.index(i+1)]
                    if i+1 in bp_y: bp_3p = bp_x[bp_y.index(i+1)]

                    if bp_3p + 1 not in bp_x and bp_3p + 1 not in bp_y:
                        hc_3p = i + 1

                if i-1 in bp_x or i-1 in bp_y: #End, condition b
                    if i-1 in bp_x: bp_5p = bp_y[bp_x.index(i-1)]

                    if i-1 in bp_y: bp_5p = bp_x[bp_y.index(i-1)]

                    if bp_5p - 1 not in bp_x and bp_5p - 1 not in bp_y:
                        hc_5p = i - 1


            line = str(i) + "\t" + allseq[i-1] + "\t" + str(nb_5p) + "\t" + str(nb_3p) + "\t" + str(nt_bp) + "\t" + str(strand_counter) + "\t" + str(hc_5p) + "\t" + str(hc_3p) + "\n"

            handle.write(line)

        #Close the file. Done.
        handle.close()

if __name__ == "__main__":

    import re

    #sequences = ["AAGATTAACTTAAAAGGAAGGCCCCCCATGCGATCAGCATCAGCACTACGACTACGCGA","acctcctta","ACGTTGGCCTTCC"]
    sequences = ["AAGATTAACTTAAAAGGAAGGCCCCCCATGCGATCAGCATCAGCACTACGACTACGCGA"]

    #Complexes
    #Input: Max number of strands in a complex. Considers all possible combinations of strands, up to max #.
    #'mfe': calculate mfe? 'ordered': consider ordered or unordered complexes?
    #Other options available (see function)

    AddComplexes = []
    test = NuPACK(sequences,"rna1999")
    test.complexes(3,mfe = 1, ordered=1)

    print test

    strand_compositions = test["ordered_composition"]
    num_complexes = len(strand_compositions)
    num_strands = len(sequences)

    for counter in range(num_complexes):
        output = "Complex #" + str(counter+1) + " composition: ("
        for strand_id in strand_compositions[counter][0:num_strands-1]:
            output = output + str(strand_id) + ", "
        output = output + str(strand_compositions[counter][num_strands-1]) + ")"

        output = output + "  dG (RT ln Q): " + str(test["ordered_energy"][counter]) + " kcal/mol"
        output = output + "  # Permutations: " + str(test["ordered_permutations"][counter])
        print output
        test.export_PDF(counter, name = "Complex #" + str(counter+1), filename = "Complex_" + str(counter) + ".pdf", program = "ordered")

    #Mfe
    #Input: Number of each strand in complex.
    #Options include RNA/DNA model, temperature, dangles, etc. (See function).
    #Example: If there are 3 unique strands (1, 2, 3), then [1, 2, 3] is one of each strand and [1, 1, 2, 2, 3, 3] is two of each strand.

    #test.mfe([1, 2], dangles = "all")
    #num_complexes = test["mfe_NumStructs"]  #Number of degenerate complexes (same energy)
    #dG_mfe = test["mfe_energy"]
    #print "There are ", num_complexes, " configuration(s) with a minimum free energy of ", dG_mfe, " kcal/mol."




